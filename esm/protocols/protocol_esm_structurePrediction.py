# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import os, glob, shutil, subprocess

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct

from .. import Plugin as esmPlugin
from ..constants import ESM_DIC

scriptName = 'runESMFold.py'

class ProtESMFoldPrediction(EMProtocol):
  """Run a structural prediction using a ESMFold model over a protein sequence"""
  _label = 'esm fold structure prediction'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)
    self.stepsExecutionMode = params.STEPS_PARALLEL

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSequence', params.PointerParam, pointerClass="Sequence",
                    label='Input protein sequence: ',
                    help="Protein sequences to perform the structur prediction on")

    mGroup = form.addGroup('Model')
    mGroup.addParam('modelName', params.EnumParam, choices=['esmfold_v1'],
                    label='Model to use: ', default=0,
                    help='Choose a model for structure prediction. \nCurrently, only v1 is available')

  def _insertAllSteps(self):
    self._insertFunctionStep(self.predictStep)
    self._insertFunctionStep(self.createOutputStep)

  def predictStep(self):
    sequence = self.getInputSequence()
    model = self.getEnumText('modelName')
    seqName = self.getInputName()
    cwd = os.path.join(esmPlugin.getVar(ESM_DIC['home']), 'esm')

    args = f'-i {sequence} -m {model} -o {seqName} -od {os.path.abspath(self._getPath())} '
    esmPlugin.runScript(self, scriptName, args, envDict=ESM_DIC, cwd=cwd)

  def createOutputStep(self):
    fnOut = self._getPath(f'{self.getInputName()}.pdb')
    if os.path.exists(fnOut):
      target = AtomStruct(filename=fnOut)
      self._defineOutputs(outputStructure=target)

  def getInputSequence(self):
    return self.inputSequence.get().getSequence()

  def getInputName(self):
    return self.inputSequence.get().getSeqName()