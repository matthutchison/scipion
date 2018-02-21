# **************************************************************************
# *
# * Authors:     Tomas Majtner (tmajtner@cnb.csic.es)
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

import os
from glob import glob
import pyworkflow.em as em
import pyworkflow.protocol.params as params
from pyworkflow.utils.path import cleanPath, cleanPattern

SIMPLE_HOME = 'SIMPLE_HOME'
SIMPLE_PRIME2D = "simple_distr_exec prg=prime2D_stream"  # streaming version


class ProtPrime2D(em.ProtClassify2D):
    """ Executes prime2d algorithm for 2d clustering of particles """
    _label = 'prime2d'

    #--------------------------- DEFINE param functions ---------------------
    
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputParts', params.PointerParam,
                      label="Input particles", important=True,
                      pointerClass='SetOfParticles',
                      help='Select the input set of particles from '
                           'the project.')
        form.addParam('msk', params.IntParam, default=36,
                      label='Mask radius', help="Mask radius in pixels")
        form.addParam('ncls_start', params.IntParam,
                      label='Number of clusters', default=5,
                      help="# of clusters required to start prime2D streaming")
        form.addParam('nparts', params.IntParam, default=2,
                      label='Number of partitions',
                      help="# of partitions in distributed exection")
        form.addParam('nptcls_per_cls', params.IntParam, default=10,
                      label='Number of images pre class',
                      help="# of images per class for 2D streaming")
        form.addParam('smpd', params.FloatParam, default=2.43,
                      label='Sampling rate',
                      help="Sampling distance, same as EMANs apix(in A)")

        form.addParam('doautoscale', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Automatic down-scaling?',
                      help='Automatic down-scaling (yes|no)')
        form.addParam('autoscale', params.BooleanParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Automatic down-scaling', condition='doautoscale',
                      help="Automatic down-scaling (yes|no)")
        form.addParam('cenlp', params.FloatParam, default=30.,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Low-pass limit for binarisation',
                      help="Low-pass limit for binarisation in "
                           "centering (in A)")
        form.addParam('center', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Center particles',
                      help='Center image(s)/class average(s)/volume(s)')
        form.addParam('dofilwidth', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Width of filament?',
                      help='Width of filament (in Angstroms).')
        form.addParam('filwidth', params.FloatParam, label='Width of filament',
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='dofilwidth',
                      help="Width of filament (in Angstroms)")
        form.addParam('dohp', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='High-pass limit?',
                      help='High-pass limit (in A).')
        form.addParam('hp', params.FloatParam, label='High-pass limit',
                      expertLevel=params.LEVEL_ADVANCED, condition='dohp',
                      help="High-pass limit (in A)")
        form.addParam('doinner', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Inner mask radius?',
                      help='Inner mask radius (in pixels).')
        form.addParam('inner', params.IntParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Inner mask radius', condition='doinner',
                      help="Inner mask radius (in pixels)")
        form.addParam('dolp', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Low-pass limit?', help='Low-pass limit (in A).')
        form.addParam('lp', params.FloatParam,
                      expertLevel=params.LEVEL_ADVANCED, condition='dolp',
                      label='Low-pass limit', help="Low-pass limit (in A)")
        form.addParam('match_filt', params.BooleanParam, default=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Matched filter on', help='Matched filter on')
        form.addParam('nthr', params.IntParam, default=1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='# of OpenMP threads',
                      help="Number of OpenMP threads")
        form.addParam('objfun', params.EnumParam, choices=['cc', 'ccres'],
                      expertLevel=params.LEVEL_ADVANCED, default=0,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label="Objective function",
                      help="Number of OpenMP threads")
        form.addParam('opt', params.EnumParam, choices=['bfgs', 'simplex'],
                      expertLevel=params.LEVEL_ADVANCED,
                      default=0, display=params.EnumParam.DISPLAY_HLIST,
                      label="Optimiser", help="Optimiser")
        form.addParam('phaseplate', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Images obtained with Volta phaseplate',
                      help='Images obtained with Volta phaseplate')
        form.addParam('dotrs', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Maximum halfwidth shift?',
                      help='Maximum halfwidth shift in pixels.')
        form.addParam('trs', params.FloatParam, condition='dotrs',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Maximum halfwidth shift',
                      help="Maximum halfwidth shift in pixels")
        form.addParam('weights2D', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Spectral weighting in 2D',
                      help='Spectral weighting in 2D')
        form.addParam('width', params.IntParam, default=10,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Falloff of inner mask or filter',
                      help="Falloff of inner mask or filter (in pixels)")


        form.addParallelSection(threads=4, mpi=0)
    
    #--------------------------- INSERT steps functions ----------------------
    
    def _insertAllSteps(self):
        """ Mainly prepare the command line for calling simple prime program"""
        
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runPrime2D')
        self._insertFunctionStep('createOutputStep')        

    #--------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        self.inputParts.get().writeStack(self._getExtraPath("ptcls_from_1.mrc"))

        line = ""
        for p in self.inputParts.get():
            line = line + "kv=300 fraca=0.1 cs=2.7 dfx=%f dfy=%f angast=%f\n" % \
                          (#p.!acceleration voltage,
                           #p.!fraction of amplitude contrast,
                           #p.!spherical abberation,
                           p._ctfModel._defocusU,
                           p._ctfModel._defocusV,
                           p._ctfModel._defocusAngle)
        params = open(self._getExtraPath("params.txt"), "w")
        params.write(line)
        params.close()

    def runPrime2D(self):
        args = "ctf=no dir_ptcls=extra deftab=extra/params.txt msk=%f " \
               "ncls_start=%d nparts=%d nptcls_per_cls=%d smpd=%f" \
               % (self.msk, self.ncls_start, self.nparts, self.nptcls_per_cls,
                  self.smpd)
        
        if self.doautoscale:
            if self.autoscale:
                args += " autoscale=yes"
            else:
                args += " autoscale=no"
        args += " cenlp=%f" % self.cenlp
        if self.center:
            args += " center=yes"
        else:
            args += " center=no"
        if self.dofilwidth:
            args += " filwidth=%f" % self.filwidth
        if self.dohp:
            args += " hp=%f" % self.hp
        if self.doinner:
            args += " inner=%f" % self.inner
        if self.dolp:
            args += " lp=%f" % self.lp
        if self.match_filt:
            args += " match_filt=yes"
        else:
            args += " match_filt=no"
        args += " nthr=%d" % self.nthr
        if self.objfun == 0:
            args += " objfun=cc"
        else:
            args += " objfun=ccres"
        if self.opt == 0:
            args += " opt=bfgs"
        else:
            args += " opt=simplex"
        if self.phaseplate:
            args += " phaseplate=yes"
        else:
            args += " phaseplate=no"
        if self.dotrs:
            args += " trs=%f" % self.trs
        if self.weights2D:
            args += " weights2D=yes"
        else:
            args += " weights2D=no"
        args += " width=%f" % self.width
        
        self.runJob(SIMPLE_PRIME2D, args, cwd=self._getPath())

    def getLastIteration(self):
        lastIter = 1
        pattern = self._getExtraPath("recvol_state1_iter%d.spi")
        while os.path.exists(pattern % lastIter):
            lastIter += 1
        return lastIter - 1

    
    def createOutputStep(self):
        lastIter = self.getLastIteration()
        
        if lastIter <= 1:
            return
        
        if self.Nvolumes == 1:
            vol = em.Volume()
            vol.setLocation(self._getExtraPath('recvol_state1_iter%d.spi' % lastIter))
            vol.setSamplingRate(self.inputClasses.get().getSamplingRate())
            self._defineOutputs(outputVol=vol)
        else:
            vol = self._createSetOfVolumes()
            vol.setSamplingRate(self.inputClasses.get().getSamplingRate())
            fnVolumes=glob(self._getExtraPath('recvol_state*_iter%d.spi') % lastIter)
            fnVolumes.sort()
            for fnVolume in fnVolumes:
                aux=em.Volume()
                aux.setLocation(fnVolume)
                vol.append(aux)
            self._defineOutputs(outputVolumes=vol)

        self._defineSourceRelation(self.inputClasses, vol)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        summary.append("Input particles: %s" % self.getObjectTag('inputParts'))
        return summary
    
    def _citations(self):
        return ['Elmlund2013']
    
    def _methods(self):
        if self.inputParts.get() is not None:
            retval="We used *simple_distr_exec prg=prime2D_stream* program " \
                   "[Elmlund2013] to produce set of classes from " \
                   "input set of particles %s."
            return [retval % self.getObjectTag('inputParts')]
        else:
            return []
