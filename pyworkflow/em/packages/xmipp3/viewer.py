# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO, CommandView
from pyworkflow.em.data import *
from pyworkflow.em.protocol import *

from xmipp3 import getXmippPath, getEnviron
from pyworkflow.em.data_tiltpairs import MicrographsTiltPair, ParticlesTiltPair, CoordinatesTiltPair
from convert import *
from os.path import dirname, join
from pyworkflow.utils import makePath, runJob, copyTree, cleanPath
import pyworkflow as pw
import xmipp
import pyworkflow.gui.dialog as dialog

from protocol_cl2d_align import XmippProtCL2DAlign
from protocol_cl2d import XmippProtCL2D
from protocol_compare_reprojections import XmippProtCompareReprojections
from protocol_compare_angles import XmippProtCompareAngles
from protocol_ctf_discrepancy import XmippProtCTFDiscrepancy
from protocol_extract_particles import XmippProtExtractParticles
from protocol_extract_particles_pairs import XmippProtExtractParticlesPairs
from protocol_helical_parameters import XmippProtHelicalParameters
from protocol_kerdensom import XmippProtKerdensom
from protocol_particle_pick_automatic import XmippParticlePickingAutomatic
from protocol_particle_pick import XmippProtParticlePicking
from protocol_particle_pick_pairs import XmippProtParticlePickingPairs
from protocol_preprocess import XmippProtPreprocessVolumes
from protocol_preprocess_micrographs import XmippProtPreprocessMicrographs
from protocol_rotational_spectra import XmippProtRotSpectra
from protocol_screen_particles import XmippProtScreenParticles
from protocol_ctf_micrographs import XmippProtCTFMicrographs
from pyworkflow.em.showj import *
from protocol_validate_nontilt import XmippProtValidateNonTilt
from protocol_multireference_alignability import XmippProtMultiRefAlignability
from protocol_assignment_tilt_pair import XmippProtAssignmentTiltPair
from protocol_movie_gain import XmippProtMovieGain


class XmippViewer(Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [
                XmippProtCompareReprojections,
                XmippProtCompareAngles,
                XmippParticlePickingAutomatic,
                XmippProtExtractParticles,
                XmippProtExtractParticlesPairs,
                XmippProtKerdensom,
                XmippProtParticlePickingPairs,
                XmippProtRotSpectra,
                XmippProtScreenParticles,
                XmippProtCTFMicrographs,
                XmippProtValidateNonTilt,
                XmippProtAssignmentTiltPair,
                XmippProtMultiRefAlignability,
                XmippProtMovieGain
                ]

    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)
    
    def __createTemporaryCtfs(self, obj, setOfMics):
        pwutils.cleanPath(obj._getPath("ctfs_temporary.sqlite"))
        self.protocol._createFilenameTemplates()
        ctfSet = self.protocol._createSetOfCTF("_temporary")

        for mic in setOfMics:
            micDir = obj._getExtraPath(removeBaseExt(mic.getFileName()))
            ctfparam = self.protocol._getFileName('ctfparam', micDir=micDir)

            if exists(ctfparam) or exists('xmipp_default_ctf.ctfparam'):
                if not os.path.exists(ctfparam):
                    ctfparam = 'xmipp_default_ctf.ctfparam'
                ctfModel = readCTFModel(ctfparam, mic)
                self.protocol._setPsdFiles(ctfModel, micDir)
                ctfSet.append(ctfModel)

        if not ctfSet.isEmpty():
            ctfSet.write()
            ctfSet.close()

        return ctfSet

    def _visualize(self, obj, **kwargs):
        cls = type(obj)

        if issubclass(cls, XmippProtCTFMicrographs):
            if obj.hasAttribute('outputCTF'):
                ctfSet = obj.outputCTF
            else:
                mics = obj.inputMicrographs.get()
                ctfSet = self.__createTemporaryCtfs(obj, mics)

            if ctfSet.isEmpty():
                self._views.append(self.infoMessage("No CTF estimation has finished yet"))
            else:
                self._views.append(CtfView(self._project, ctfSet))

        elif (issubclass(cls, XmippProtExtractParticles) or
              issubclass(cls, XmippProtScreenParticles)):
            particles = obj.outputParticles
            self._visualize(particles)

            fn = obj._getPath('images.xmd')
            if os.path.exists(fn): # it doesnt unless cls is Xmipp
                md = xmipp.MetaData(fn)
                # If Zscore on output images plot Zscore particle sorting
                if md.containsLabel(xmipp.MDL_ZSCORE):
                    from plotter import XmippPlotter
                    xplotter = XmippPlotter(windowTitle="Zscore particles sorting")
                    xplotter.createSubPlot("Particle sorting", "Particle number", "Zscore")
                    xplotter.plotMd(md, False, mdLabelY=xmipp.MDL_ZSCORE)
                    self._views.append(xplotter)

        elif issubclass(cls, XmippProtMovieGain):
            self._visualize(obj.outputMovies)
            movieGainMonitor = MonitorMovieGain(obj,
                                                workingDir=obj.workingDir.get(),
                                                samplingInterval=60,
                                                monitorTime=300,
                                                stddevValue=0.04,
                                                ratio1Value=1.15,
                                                ratio2Value=4.5)
            self._views.append(MovieGainMonitorPlotter(movieGainMonitor))

        elif issubclass(cls, XmippProtRotSpectra):
            self._visualize(obj.outputClasses,
                            viewParams={'columns': obj.SomXdim.get(),
                                        RENDER: ' spectraPlot._filename average._filename',
                                        ZOOM: 30,
                                        VISIBLE:  'enabled id _size average._filename spectraPlot._filename',
                                        'labels': 'id _size',
                                        SORT_BY: 'id'})

        elif issubclass(cls, XmippProtKerdensom):
            self._visualize(obj.outputClasses,
                            viewParams={'columns': obj.SomXdim.get(),
                                       'render': 'average._filename _representative._filename',
                                       'labels': '_size',
                                       'sortby': 'id'})



        elif issubclass(cls, XmippProtCompareReprojections):
                fn = obj.outputParticles.getFileName()
                labels = 'id enabled _index _xmipp_image._filename _xmipp_imageRef._filename _xmipp_imageResidual._filename _xmipp_imageCovariance._filename _xmipp_cost _xmipp_zScoreResCov _xmipp_zScoreResMean _xmipp_zScoreResVar _xmipp_continuousA _xmipp_continuousB _xmipp_continuousX _xmipp_continuousY'
                labelRender = "_xmipp_image._filename _xmipp_imageRef._filename _xmipp_imageResidual._filename _xmipp_imageCovariance._filename"
                self._views.append(ObjectView(self._project, obj.outputParticles.strId(), fn,
                                              viewParams={ORDER: labels,
                                                      VISIBLE: labels,
                                                      SORT_BY: '_xmipp_cost asc', RENDER:labelRender,
                                                      MODE: MODE_MD}))

        elif issubclass(cls, XmippProtCompareAngles):
                fn = obj.outputParticles.getFileName()
                labels = 'id enabled _index _filename _xmipp_shiftDiff _xmipp_angleDiff'
                labelRender = "_filename"
                self._views.append(ObjectView(self._project, obj.outputParticles.strId(), fn,
                                              viewParams={ORDER: labels,
                                                      VISIBLE: labels,
                                                      SORT_BY: '_xmipp_angleDiff asc', RENDER:labelRender,
                                                      MODE: MODE_MD}))

        elif issubclass(cls, XmippParticlePickingAutomatic):
            micSet = obj.getInputMicrographs()
            mdFn = getattr(micSet, '_xmippMd', None)
            inTmpFolder = False
            if mdFn:
                micsfn = mdFn.get()
            else:  # happens if protocol is not an xmipp one
                micsfn = self._getTmpPath(micSet.getName() + '_micrographs.xmd')
                writeSetOfMicrographs(micSet, micsfn)
                inTmpFolder = True

            posDir = obj._getExtraPath()
            memory = '%dg'%obj.memory.get(),
            launchSupervisedPickerGUI(micsfn, posDir, obj, mode='review', memory=memory, inTmpFolder=inTmpFolder)

         # We need this case to happens before the ProtParticlePicking one
        elif issubclass(cls, XmippProtAssignmentTiltPair):
            if obj.getOutputsSize() >= 1:
                coordsSet = obj.getCoordsTiltPair()
                self._visualize(coordsSet)

        elif issubclass(cls, XmippProtValidateNonTilt):
            outputVols = obj.outputVolumes
            labels = 'id enabled comment _filename weight'
            self._views.append(ObjectView(self._project, outputVols.strId(), outputVols.getFileName(),
                                          viewParams={MODE: MODE_MD, VISIBLE:labels, ORDER: labels,
                                                      SORT_BY: 'weight desc', RENDER: '_filename'}))

        elif issubclass(cls, XmippProtMultiRefAlignability):
            outputVols = obj.outputVolumes
            labels = 'id enabled comment _filename weightAlignabilityPrecision weightAlignabilityAccuracy'
            self._views.append(ObjectView(self._project, outputVols.strId(), outputVols.getFileName(),
                                          viewParams={MODE: MODE_MD, VISIBLE:labels, ORDER: labels,
                                                      SORT_BY: 'weightAlignabilityAccuracy desc', RENDER: '_filename'}))

            fn = obj.outputParticles.getFileName()
            labels = 'id enabled _index _filename _xmipp_scoreAlignabilityAccuracy _xmipp_scoreAlignabilityPrecision'
            labelRender = "_filename"
            self._views.append(ObjectView(self._project, obj.outputParticles.strId(), fn,
                                            viewParams={ORDER: labels,
                                                      VISIBLE: labels,
                                                      SORT_BY: '_xmipp_scoreAlignabilityAccuracy desc', RENDER:labelRender,
                                                      MODE: MODE_MD}))

            fn = obj._getExtraPath('vol001_pruned_particles_alignability.xmd')
            md = xmipp.MetaData(fn)
            from plotter import XmippPlotter
            from pyworkflow.em.plotter import EmPlotter
            plotter = XmippPlotter()
            plotter.createSubPlot('Soft-alignment validation plot','Angular Precision', 'Angular Accuracy')
            plotter.plotMdFile(md, xmipp.MDL_SCORE_BY_ALIGNABILITY_PRECISION, xmipp.MDL_SCORE_BY_ALIGNABILITY_ACCURACY,
                               marker='.', markersize=.55, color='red', linestyle='')
            self._views.append(plotter)

        elif issubclass(cls, XmippProtExtractParticlesPairs):
            self._visualize(obj.outputParticlesTiltPair)

        return self._views




