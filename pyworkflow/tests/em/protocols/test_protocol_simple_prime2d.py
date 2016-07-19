# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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


import pyworkflow.em.packages.simple as simple
from pyworkflow.em import ImageHandler
from pyworkflow.em.protocol import ProtImportParticles

from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.tests.em.workflows.test_workflow import TestWorkflow

  
   
class TestSimplePrime(TestWorkflow):
    @classmethod
    def setUpClass(cls):    
        # Create a new project
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('mda')
        cls.particlesFn = cls.dataset.getFile('particles')
        protImport = cls.newProtocol(ProtImportParticles,
                                      filesPath=cls.particlesFn,
                                      samplingRate=3.5)

        cls.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles '
                            'is None.' % cls.particlesFn)

        cls.protImport = protImport


    def test_prime2D(self):
        n = 4
        protClassify = self.newProtocol(simple.ProtPrime2D,
                                        generateReferences=True,
                                        numberOfClasses=n,
                                        maskRadius=45)
        protClassify.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protClassify)
        classes = getattr(protClassify, 'outputClasses', None)
        self.assertIsNotNone(classes, "There was a problem with the "
                                      "simple - prime 2D outputClasses")

        # Do some basic validations
        self.assertEqual(classes.getSize(), n)
        firstClass = classes.getFirstItem()
        ih = ImageHandler()
        self.assertTrue(ih.existsLocation(firstClass.getRepresentative()))

        return protClassify

    def test_prime3D_initialFromParticles(self):
        protInitial = self.newProtocol(simple.ProtPrime3DInitial,
                                       objLabel='from particles',
                                       maskRadius=45)
        protInitial.inputSet.set(self.protImport.outputParticles)
        self.launchProtocol(protInitial)

    def test_prime3D_initialFromClasses(self):
        protClassify = self.test_prime2D()

        protInitial = self.newProtocol(simple.ProtPrime3DInitial,
                                       objLabel='from classes',
                                       maskRadius=45)
        protInitial.inputSet.set(protClassify.outputClasses)
        self.launchProtocol(protInitial)



