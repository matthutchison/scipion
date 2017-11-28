# **************************************************************************
# *
# * Authors:  Ruben Sanchez (rsanchez@cnb.csic.es), April 2017
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from __future__ import print_function

from six.moves import range
import numpy as np
import sys, os

import numpy as np
import scipy
import random
from math import ceil

import time
from os.path import expanduser
from pyworkflow.utils import Environ
from sklearn.metrics import roc_auc_score
import xmipp
import pyworkflow.em.metadata as md

DEBUG=False
if DEBUG: print("Debug MODE")

def updateEnviron(gpuNum=None):
  """ Create the needed environment for TensorFlow programs. """
  environ = Environ(os.environ)
  if not gpuNum is None:
    environ.update({'LD_LIBRARY_PATH': os.environ['CUDA_LIB']}, position=Environ.BEGIN)
    environ.update({'LD_LIBRARY_PATH': os.environ['CUDA_HOME']+"/extras/CUPTI/lib64"}, position=Environ.BEGIN)

    os.environ['CUDA_VISIBLE_DEVICES']=str(gpuNum)  #THIS IS FOR USING JUST GPU:# must be changed to select desired gpu

import tensorflow as tf
import tflearn
from pyworkflow.em.packages.xmipp3.networkDef import main_network
tf_intarnalError= tf.errors.InternalError

BATCH_SIZE= 128

EVALUATE_AT= 10
CHECK_POINT_AT= 100



class DeepTFSupervised(object):
  def __init__(self, rootPath, modelNum=0):
    '''
      @param rootPath: str. Root directory where neural net data will be saved.
                            Generally "extra/nnetData/"
                                                      tfchkpoints/
                                                      tflogs/
      @param learningRate: float. Learning rate for net training

    '''
    self.lRate= None
    self.rootPath= rootPath

    self.checkPointsNames= os.path.join(rootPath,"tfchkpoints_%d"%modelNum)
    self.checkPointsNames= os.path.join(self.checkPointsNames,"screening")

    self.logsSaveName= os.path.join(rootPath,"tflogs_%d"%modelNum)

    self.num_labels=2
    self.num_channels=None
    self.image_size=None

    # Tensorflow objects
    self.X= None
    self.Y= None
    self.saver = None
    self.session = None
    self.train_writer = None
    self.test_writer = None
    self.global_step= None
    self.loss= None
    self.optimizer= None

  def createNet(self, xdim, ydim, num_chan, num_extra_feat, nData=2**12):
    '''
      @param xdim: int. height of images
      @param ydim: int. width of images
      @param num_chan: int. number of channels of images
      @param num_extra_feat: int. number of extra features computed for each image
      @param nData: number of positive cases expected in data. Not needed
    '''
    print ("Creating net.")
    ##############################################################
    # INTIALIZE INPUT PLACEHOLDERS
    ##############################################################
    self.num_channels= num_chan
    self.image_size= (xdim, ydim)
    num_labels= self.num_labels
    tflearn.config.init_training_mode()
    self.X= tf.placeholder( tf.float32, shape=(None, self.image_size[0], self.image_size[1], num_chan), name="X")
    self.X_feats= tf.placeholder( tf.float32, shape=(None, num_extra_feat), name="X_feats")
    self.Y= tf.placeholder(tf.float32, shape=(None, num_labels), name="Y")
    self.lRate= tf.placeholder(tf.float32, name="learningRate")
    ######################
    # NEURAL NET CREATION
    ######################
    self.global_step = tf.Variable(initial_value=0, name='global_step', trainable=False)
    (self.y_pred, self.merged_summaries, self.optimizer,
     self.loss, self.accuracy)= main_network(self.X, self.X_feats, self.Y, num_labels, self.lRate, globalStep= self.global_step, nData= nData)

  def startSessionAndInitialize(self, numberOfThreads=8):
    '''
      @param numberOfThreads. Number of threads to use in cpu mode. If numberOfThreads==None
                              then default behaviour is expected (use GPU if available otherwise
                              one thread per cpu)
    '''
    print("Initializing tf session")

    save_dir, prefixName = os.path.split(self.checkPointsNames)
    if not os.path.exists(save_dir):
      os.makedirs(save_dir)

    self.saver = tf.train.Saver()
    print("numberOfThreads",numberOfThreads)
    if numberOfThreads is None:
      self.session = tf.Session()
    else:
      self.session= tf.Session(config=tf.ConfigProto(intra_op_parallelism_threads=numberOfThreads))

    try:
      print("Trying to restore last checkpoint from "+ os.path.abspath(save_dir))
      # Use TensorFlow to find the latest checkpoint - if any.
      last_chk_path = tf.train.latest_checkpoint(checkpoint_dir= save_dir)

      # Try and load the data in the checkpoint.

      self.saver.restore(self.session, save_path=last_chk_path)
      # If we get to this point, the checkpoint was successfully loaded.
      print("Restored checkpoint from:", last_chk_path)
    except Exception as e:
      print("last chk_point path:",last_chk_path)
      print(e)
      # If the above failed for some reason, simply
      # initialize all the variables for the TensorFlow graph.
      print("Failed to restore checkpoint. Initializing variables instead.")
      if( os.path.isdir(self.logsSaveName)):
        tf.gfile.DeleteRecursively(self.logsSaveName)
      os.makedirs(os.path.join(self.logsSaveName,"train"))
      os.makedirs(os.path.join(self.logsSaveName,"test"))
      self.session.run( tf.global_variables_initializer() )

    self.train_writer = tf.summary.FileWriter(os.path.join(self.logsSaveName,"train"), self.session.graph)
    self.test_writer = tf.summary.FileWriter(os.path.join(self.logsSaveName,"test"), self.session.graph)
    return self.session

  def close(self, saveModel= False):
    '''
      Closes a tensorflow connection and related objects.
      @param. saveModel: boolean. If True, model will be saved prior closing model

    '''
    if saveModel:
      self.saver.save(self.session, save_path= self.checkPointsNames, global_step= self.global_step)
      print("\nSaved checkpoint.")
    self.train_writer.close()
    self.test_writer.close()
    tf.reset_default_graph()
    self.session.close()


  def reset(self):
    '''
      Resets a tensorflow connection and related objects.
      Needed if 2 independent trains want to be done in the same program.

    '''
    self.train_writer.close()
    self.test_writer.close()
    self.session.close()
    tf.reset_default_graph()
    self.createNet()
    self.startSessionAndInitialize()

  def trainNet(self, numberOfBatches, dataManagerTrain, learningRate, dataManagerTest=None, auto_stop=False):
    '''
      @param numberOfBatches: int. The number of batches that will be used for training
      @param dataManagerTrain: DataManager. Object that will provide training batches (Xs and labels)
      @param dataManagerTest:  DataManager. Optional Object that will provide testing batches (Xs and labels)
                                            If not provided, no testing will be done
    '''

    ########################
    # TENSOR FLOW RUN
    ########################
    numberOfBatches= max( numberOfBatches, CHECK_POINT_AT+1)
    print("Learning rate %.1e"% learningRate)

    learningRate_At_Convergency= 0.01* learningRate
    print("auto_stop:", auto_stop)
    batchsPerEpoch= dataManagerTrain.getEpochSize()// dataManagerTrain.getBatchSize() +1
    saving_checkpoint_at= min(CHECK_POINT_AT, batchsPerEpoch)
    hasImproved=False
    numEpochsNoImprov=0
    epochImprovement= 0
    numEpochsNoImprov_Limit= 2
    if batchsPerEpoch < 5:
      numEpochsNoImprov_Limit= 10
    elif batchsPerEpoch< 10:
      numEpochsNoImprov_Limit= 5
    elif  batchsPerEpoch< 20:
      numEpochsNoImprov_Limit= 4

    trainDataBatch= dataManagerTrain.getRandomBatch()
    x_batchTrain, x_featsTrain, labels_batchTrain, md_ids = trainDataBatch
    feed_dict_train= {self.X : x_batchTrain, self.X_feats: x_featsTrain, self.Y: labels_batchTrain,
                        self.lRate: learningRate}
    tflearn.is_training(False, session=self.session)
    i_global,stepLoss= self.session.run( [self.global_step, self.loss], feed_dict= feed_dict_train)
    numberOfRemainingBatches= max(0, numberOfBatches- i_global)
    bestStepLoss= 2^30
    currentLoss=[]
    modelWasSaved=False
    if numberOfBatches >0:
      print("Training net for %d batches of size %d"%(numberOfRemainingBatches, dataManagerTrain.getBatchSize()))
      print("Initial loss %f"%stepLoss)
      self.testPerformance(i_global, trainDataBatch,dataManagerTest)
    else:
      return

    time0 = time.time()
    tflearn.is_training(True, session=self.session)
    for iterNum in range(numberOfRemainingBatches):
      trainDataBatch= dataManagerTrain.getRandomBatch()
      x_batchTrain, x_featsTrain, labels_batchTrain, md_ids = trainDataBatch

      feed_dict_train= {self.X : x_batchTrain, self.X_feats: x_featsTrain, self.Y: labels_batchTrain, self.lRate: learningRate }
      i_global, __, stepLoss, y_pred= self.session.run( [self.global_step, self.optimizer, self.loss,
                                                  self.y_pred],
                                                  feed_dict=feed_dict_train )

      currentLoss.append( stepLoss)

      print("iterNum %d/%d trainLoss: %3.4f"%((i_global), numberOfBatches, stepLoss))
      sys.stdout.flush()
      if dataManagerTest and i_global % EVALUATE_AT ==0:
        timeBatches= time.time() -time0
        timeTest=time.time()
        self.testPerformance(i_global,trainDataBatch,dataManagerTest)
        print("%d batches time: %f s. time for test %f s"%(EVALUATE_AT, timeBatches, time.time() -timeTest))
        time0 = time.time()
        sys.stdout.flush()

      if (i_global + 1 ) % saving_checkpoint_at==0:
        stepLossMean= np.min(currentLoss)
        currentLoss=[]
        improvement= bestStepLoss-stepLossMean
        print("Training improvement since last checkpoint %.5f"%(improvement))
        if improvement>0:
          epochImprovement= improvement
          modelWasSaved=True
          hasImproved= True
          bestStepLoss= stepLossMean
          self.saver.save(self.session, save_path= self.checkPointsNames, global_step= self.global_step)
          print("\nSaved checkpoint.")

      if auto_stop and ((i_global+1)% batchsPerEpoch)==0:
        print("Epoch %d finished. Learning rate %.1e. Epoch improvement: %f"%(i_global//batchsPerEpoch, learningRate, epochImprovement))
        epochImprovement=0
        if not hasImproved:
          numEpochsNoImprov+=1
          if numEpochsNoImprov== numEpochsNoImprov_Limit:
            learningRate*= 0.1
            print("reducing learning rate to %.1e"%learningRate)
            numEpochsNoImprov=0
          if learningRate < learningRate_At_Convergency:
            print("CONVERGENCY AUTO-DETECTED")
            break
        else:
          numEpochsNoImprov=0
          hasImproved= False

    self.testPerformance(i_global,trainDataBatch, dataManagerTest)
    if not modelWasSaved:
      self.saver.save(self.session, save_path= self.checkPointsNames, global_step= self.global_step)
      print("\nSaved checkpoint.")

  def accuracy_score(self, labels, predictions ):
    return (100.0 * np.sum(np.argmax(predictions, 1) == np.argmax(labels, 1))
            / predictions.shape[0])

  def testPerformance(self, stepNum, trainDataBatch, testDataManager=None):
    tflearn.is_training(False, session=self.session)

    batch_x, feats_x, batch_y, md_ids= trainDataBatch
    feed_dict= {self.X : batch_x, self.X_feats: feats_x, self.Y: batch_y}
    c_e_train, y_pred_train, merged = self.session.run( [self.loss, self.y_pred, self.merged_summaries],
                                                 feed_dict = feed_dict )
    train_accuracy=  self.accuracy_score(batch_y, y_pred_train)
    train_auc= roc_auc_score(batch_y, y_pred_train)
    self.train_writer.add_summary(merged, stepNum)
    if not testDataManager is None:
      y_pred_list=[]
      labels_list=[]
      for images, feats, labels in testDataManager.getIteratorTestBatch(20):
        feed_dict= {self.X : images, self.X_feats: feats, self.Y : labels }
        y_pred, merged= self.session.run( [self.y_pred,self.merged_summaries], feed_dict=feed_dict )
        y_pred_list.append(y_pred)
        labels_list.append(labels)

      y_pred= np.concatenate(y_pred_list)
      labels= np.concatenate(labels_list)
      test_accuracy= self.accuracy_score(labels, y_pred)
      test_auc= roc_auc_score(labels, y_pred)

      self.test_writer.add_summary(merged, stepNum)
      print("iterNum %d accuracy(train/test) %f / %f   auc  %f / %f"%(stepNum,train_accuracy,test_accuracy,
                                                                      train_auc, test_auc))
    tflearn.is_training(True, session=self.session)


  def predictNet(self, dataManger):
    tflearn.is_training(False, session=self.session)
    y_pred_list=[]
    labels_list=[]
    metadataId_list=[]
    for images, feats, labels, metadataIdTuple in dataManger.getIteratorPredictBatch():
      feed_dict= {self.X : images, self.X_feats: feats}
      y_pred= self.session.run( self.y_pred, feed_dict=feed_dict )
      y_pred_list.append(y_pred)
      labels_list.append(labels)
      metadataId_list+= metadataIdTuple
    y_pred= np.concatenate(y_pred_list)
    labels= np.concatenate(labels_list)
    y_pred_oneCol= y_pred[:,1]
    return y_pred_oneCol, labels, metadataId_list

class DataManager(object):
  COMP_FEAT_TYPES=[md.MDL_SCORE_BY_SCREENING, md.MDL_SCORE_BY_ZERNIKE, md.MDL_SCORE_BY_ENTROPY,
                   md.MDL_SCORE_BY_LBP, md.MDL_SCORE_BY_VARIANCE, md.MDL_SCORE_BY_RAMP ]

#  COMP_FEAT_TYPES=[md.MDL_SCORE_BY_SCREENING ]

  def __init__(self, posSetDict, negSetDict=None):
    '''
        posSetDict, negSetDict: { fname: [(SetOfParticles, weight:int)]
    '''
    self.mdListFalse=None
    self.nFalse=0 #Number of negative particles in dataManager

    (self.mdListTrue, self.fnListOfListTrue, self.featsListTrue, self.weightListTrue, self.nTrue,
         self.shape, self.nCompFeats)= self.colectMetadata(posSetDict)
    self.weightListTrue= [ elem/float(sum(self.weightListTrue)) for elem in  self.weightListTrue]
    self.batchSize= min(BATCH_SIZE, self.nTrue)
    self.splitPoint= self.batchSize//2

    self.batchStack = np.zeros((self.batchSize,)+self.shape)
    self.batchFeats = np.zeros((self.batchSize,)+(self.nCompFeats,))

    if not negSetDict is None:
      (self.mdListFalse, self.fnListOfListFalse, self.featsListFalse, self.weightListFalse, self.nFalse,
           shapeFalse, nCompFeats)=  self.colectMetadata(negSetDict)
      self.weightListFalse= [ elem/float(sum(self.weightListFalse)) for elem in  self.weightListFalse]
      assert shapeFalse== self.shape, "Negative images and positive images have differnt shape"
      assert nCompFeats== self.nCompFeats, "Negative images and positive images have differnt feats vect shape"

      self.getRandomBatch= self.getRandomBatchWorker
    else:
      self.getRandomBatch= self.NOgetRandomBatchWorker

    if DEBUG:
      self.nTrue=  2**9
      self.nFalse= 2**9

  def colectMetadata(self, dictData):
    mdList=[]
    fnamesList=[]
    featsList=[]
    weightsList= []
    nParticles=0
    shapeParticles=(None, None, 1)

    for fnameXMDF in sorted(dictData):
      setOfParticlesObject, weight= dictData[fnameXMDF]
      mdObject  = md.MetaData(fnameXMDF)
      featsList.append( self.getAllCompFeats(mdObject))
      mdList+= [mdObject]
      fnamesList+= [mdObject.getColumnValues(md.MDL_IMAGE)]
      xdim, ydim, _   = setOfParticlesObject.getDim()
      tmpShape= (xdim,ydim,1)
      tmpNumParticles= setOfParticlesObject.getSize()
      if shapeParticles!= (None, None, 1):
        assert tmpShape== shapeParticles, "Error, particles of different shapes mixed"
      else:
        shapeParticles= tmpShape
      weightsList+= [ weight ]
      nParticles+= tmpNumParticles
    nFeats= featsList[0].shape[-1]
#    print([elem.shape for elem in featsList])
#    print(nFeats)
    return mdList, fnamesList, featsList, weightsList, nParticles, shapeParticles, nFeats

  def getAllCompFeats(self, mdObject):
      imFeats=[]
      for featId in DataManager.COMP_FEAT_TYPES:
        feat= np.array(mdObject.getColumnValues(featId))

        if len(feat.shape)==1 and not np.any(feat==np.array(None)) :
          feat= np.expand_dims(feat, axis=1)
        elif len(feat.shape)==0:
          continue
        imFeats.append(feat)
      allCompFeats= np.concatenate(imFeats, axis=1)
      featsMin= np.min(allCompFeats, axis=0)
      featsMax= np.min(allCompFeats, axis=0)
      return allCompFeats

  def getCompFeatsOneId(self, mdObject, oneId):
      imFeats=[]
      for featId in DataManager.COMP_FEAT_TYPES:
        feat= mdObject.getValue(featId, oneId)
        if feat is None:
          raise ValueError("Some external featue was not computed")
        else:
          feat=np.array(feat)
        if len(feat.shape)==0 :
          feat= np.array([feat])
        imFeats.append(feat)
      return  np.concatenate(imFeats)

  def getMetadata(self, setNumber=None) :

    if setNumber is None:
      return [mdTrue for mdTrue in self.mdListTrue], [mdFalse for mdFalse in self.mdListFalse] if self.mdListFalse else None
    else:
      mdTrue= self.mdListTrue[setNumber]
      mdFalse= self.mdListFalse[setNumber]
      return  mdTrue, mdFalse

  def getNBatches(self, nEpochs):
    return  int(ceil(2*self.nTrue*nEpochs/self.batchSize))

  def getEpochSize(self):
    return 2*self.nTrue

  def getBatchSize(self):
    return self.batchSize

  def _random_flip_leftright(self, batch):
    for i in range(len(batch)):
      if bool(random.getrandbits(1)):
        batch[i] = np.fliplr(batch[i])
    return batch

  def _random_flip_updown(self, batch):
    for i in range(len(batch)):
      if bool(random.getrandbits(1)):
        batch[i] = np.flipud(batch[i])
    return batch

  def _random_90degrees_rotation(self, batch, rotations=[0, 1, 2, 3]):
    for i in range(len(batch)):
      num_rotations = random.choice(rotations)
      batch[i] = np.rot90(batch[i], num_rotations)
    return batch

  def _random_rotation(self, batch, max_angle):
    for i in range(len(batch)):
      if bool(random.getrandbits(1)):
        # Random angle
        angle = random.uniform(-max_angle, max_angle)
        batch[i] = scipy.ndimage.interpolation.rotate(batch[i], angle,reshape=False, mode="reflect")
    return batch

  def _random_blur(self, batch, sigma_max):
    for i in range(len(batch)):
      if bool(random.getrandbits(1)):
        # Random sigma
        sigma = random.uniform(0., sigma_max)
        batch[i] =scipy.ndimage.filters.gaussian_filter(batch[i], sigma)
    return batch

  def augmentBatch(self, batch):
    if bool(random.getrandbits(1)):
      batch= self._random_flip_leftright(batch)
      batch= self._random_flip_updown(batch)
    if bool(random.getrandbits(1)):
      batch= self._random_90degrees_rotation(batch)
    if bool(random.getrandbits(1)):
      batch= self._random_rotation(batch, 10.0)
    return batch

  def NOgetRandomBatchWorker(self):
    raise ValueError("needs positive and negative images to compute random batch. Just pos provided")

  def getRandomBatchWorker(self):

    batchSize = self.batchSize
    splitPoint = self.splitPoint
    batchStack = self.batchStack
    batchFeats = self.batchFeats

    batchLabels  = np.zeros((batchSize, 2))

    dataSetNumTrue=  np.random.choice(len(self.mdListTrue), p= self.weightListTrue)
    dataSetNumFalse= np.random.choice(len(self.mdListFalse), p=self.weightListFalse)
    print("DataSets in batch %d/%d"%(dataSetNumTrue, dataSetNumFalse))

#    metadataTrue= self.mdListTrue[dataSetNumTrue]
#    metadataFalse= self.mdListFalse[dataSetNumFalse]

    fnamesTrue= self.fnListOfListTrue[dataSetNumTrue]
    fnamesFalse= self.fnListOfListFalse[dataSetNumFalse]

    featsTrue= self.featsListTrue[dataSetNumTrue]
    featsFalse= self.featsListFalse[dataSetNumFalse]

    nTrue= len(fnamesTrue)
    nFalse= len(fnamesFalse)

    idxListTrue =  np.random.choice( nTrue,  splitPoint, False)
    idxListFalse = np.random.choice( nFalse, splitPoint, False)

    I = xmipp.Image()
    n = 0
    finalIds= []
    for idx in idxListTrue:
      fnImage = fnamesTrue[idx]
      I.read(fnImage)
      batchStack[n,...]= np.expand_dims(I.getData(),-1)
      batchFeats[n,...]= featsTrue[idx, ...]
      batchLabels[n, 1]= 1
      finalIds.append(idx)
      n+=1
      if n>=splitPoint:
          break
    for idx in idxListFalse:
      fnImage = fnamesFalse[idx]
      I.read(fnImage)
      batchStack[n,...]= np.expand_dims(I.getData(),-1)
      batchFeats[n,...]= featsFalse[idx, ...]
      batchLabels[n, 0]= 1
      finalIds.append(idx)
      n+=1
      if n>=batchSize:
          break

    shuffInd= np.random.choice(n,n, replace=False)
    batchStack= batchStack[shuffInd, ...]
    batchLabels= batchLabels[shuffInd, ...]
    finalIds= [finalIds[i] for i in shuffInd]
    return self.augmentBatch(batchStack), batchFeats, batchLabels, finalIds

  def getDataAsNp(self):
    allData= self.getIteratorPredictBatch()
    x, labels, __, __= zip(* allData)
    x= np.concatenate(x)
    y= np.concatenate(labels)
    print(x.shape, y.shape)
    print(np.min(x), np.mean(x), np.max(x))
    return x,y

  def getIteratorPredictBatch(self):
    batchSize = self.batchSize
    xdim,ydim,nChann= self.shape
    batchStack = np.zeros((self.batchSize, xdim,ydim,nChann))
    batchFeats = np.zeros((self.batchSize, self.nCompFeats))
    batchLabels  = np.zeros((batchSize, 2))
    I = xmipp.Image()
    n = 0
    idAndType=[]
    for dataSetNum in range(len(self.mdListTrue)):
      mdTrue= self.mdListTrue[dataSetNum]
      for objId in mdTrue:
#        print(objId)
        fnImage = mdTrue.getValue(md.MDL_IMAGE, objId)
        I.read(fnImage)
        batchStack[n,...]= np.expand_dims(I.getData(),-1)
        batchFeats[n,...]= self.getCompFeatsOneId(mdTrue, objId)
        batchLabels[n, 1]= 1
        idAndType.append((True,objId, dataSetNum))
        n+=1
        if n>=batchSize:
          yield batchStack, batchFeats, batchLabels, idAndType
          n=0
          batchLabels  = np.zeros((batchSize, 2))
          idAndType= []
    if not self.mdListFalse is None:
      for dataSetNum in range(len(self.mdListFalse)):
        mdFalse= self.mdListFalse[dataSetNum]
        for objId in mdFalse:
          fnImage = mdFalse.getValue(md.MDL_IMAGE, objId)
          I.read(fnImage)
          batchStack[n,...]= np.expand_dims(I.getData(),-1)
          batchFeats[n,...]= self.getCompFeatsOneId(mdFalse, objId)
          batchLabels[n, 0]= 1
          idAndType.append((False,objId,dataSetNum))
          n+=1
          if n>=batchSize:
            yield batchStack, batchFeats, batchLabels,  idAndType
            n=0
            idAndType= []
            batchLabels  = np.zeros((batchSize, 2))
    if n>0:
      yield batchStack[:n,...], batchFeats[:n,...], batchLabels[:n,...], idAndType[:n]

  def getIteratorTestBatch(self, nBatches= 10):
    batchSize = self.batchSize
    xdim,ydim,nChann= self.shape
    batchStack = np.zeros((self.batchSize, xdim,ydim,nChann))
    batchFeats = np.zeros((self.batchSize, self.nCompFeats))
    batchLabels  = np.zeros((batchSize, 2))


    I = xmipp.Image()
    n = 0

    currNBatches=0
    for iterNum in range(nBatches):
      dataSetNumTrue=  np.random.choice(len(self.mdListTrue), p= self.weightListTrue)
      dataSetNumFalse= np.random.choice(len(self.mdListFalse), p=self.weightListFalse)
      featsTrue= self.featsListTrue[dataSetNumTrue]
      featsFalse= self.featsListFalse[dataSetNumFalse]
      fnListTrue, fnListFalse= self.fnListOfListTrue[dataSetNumTrue], self.fnListOfListFalse[dataSetNumFalse]

      for idx, (fnImageTrue, fnImageFalse) in enumerate(zip(fnListTrue, fnListFalse)):
        I.read(fnImageTrue)
        batchStack[n,...]= np.expand_dims(I.getData(),-1)
        batchFeats[n,...]= featsTrue[idx, ...]
        batchLabels[n, 1]= 1
        n+=1
        if n>=batchSize:
          yield batchStack, batchFeats, batchLabels
          n=0
          batchLabels  = np.zeros((batchSize, 2))
          currNBatches+=1
          if currNBatches>=nBatches:
            break
        I.read(fnImageFalse)
        batchStack[n,...]= np.expand_dims(I.getData(),-1)
        batchFeats[n,...]= featsFalse[idx, ...]
        batchLabels[n, 0]= 1
        n+=1
        if n>=batchSize:
          yield batchStack, batchFeats, batchLabels
          n=0
          batchLabels  = np.zeros((batchSize, 2))
          currNBatches+=1
          if currNBatches>=nBatches:
            break
      if n>0:
        yield batchStack[:n,...], batchFeats[:n,...], batchLabels[:n,...]
