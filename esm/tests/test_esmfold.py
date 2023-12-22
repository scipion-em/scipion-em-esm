# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportSequence
from pwem.tests.protocols.test_protocols_import_sequence import TestImportBase

from ..protocols import ProtESMFoldPrediction

class TestESMFold(TestImportBase):
    NAME = 'USER_SEQ'
    DESCRIPTION = 'User description'
    AMINOACIDSSEQ1 = 'LARKJLAKPABXZJUO********VAVAVALK'

    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)
        cls._runImportSequence()
        cls._waitOutput(cls.protImportSequence, 'outputSequence', sleepTime=5)

    @classmethod
    def _runImportSequence(cls):
        args = {'inputSequenceName': cls.NAME,
                'inputSequenceDescription': cls.DESCRIPTION,
                'inputRawSequence': cls.AMINOACIDSSEQ1
                }
        cls.protImportSequence = cls.newProtocol(ProtImportSequence, **args)
        cls.protImportSequence.setObjLabel('Import aminoacid seq')
        cls.proj.launchProtocol(cls.protImportSequence, wait=False)

    def _runESMFold(self):
        protESMFold = self.newProtocol(
            ProtESMFoldPrediction,
            inputSequence=self.protImportSequence.outputSequence)

        self.launchProtocol(protESMFold)
        return protESMFold

    def test(self):
        protESMFold = self._runESMFold()

        self._waitOutput(protESMFold, 'outputStructure', sleepTime=10)
        self.assertIsNotNone(getattr(protESMFold, 'outputStructure', None))


