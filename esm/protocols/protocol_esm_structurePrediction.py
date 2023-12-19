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
from pwem.convert.atom_struct import toCIF, AtomicStructHandler, addScipionAttribute

from .. import Plugin as esmPlugin
from ..constants import ESM_DIC

scriptName = 'runESMFold.py'

class ProtESMFoldPrediction(EMProtocol):
  """Run a structural prediction using a ESMFold model over a protein sequence"""
  _label = 'ESMFold structure prediction'
  _ATTRNAME = 'ESMFoldScore'
  _OUTNAME = 'outputStructure'
  _possibleOutputs = {_OUTNAME: AtomStruct}

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)
    self.stepsExecutionMode = params.STEPS_PARALLEL

  def _defineParams(self, form):
    form.addHidden(params.GPU_LIST, params.StringParam, default='0', label="Choose GPU IDs",
                   help="Add a list of GPU device that can be used")
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSequence', params.PointerParam, pointerClass="Sequence",
                    label='Input protein sequence: ',
                    help="Protein sequences to perform the structur prediction on")

    mGroup = form.addGroup('Model')
    mGroup.addParam('modelName', params.EnumParam, choices=['esmfold_v0', 'esmfold_v1'],
                    label='Model to use: ', default=1,
                    help='Choose a model for structure prediction. \nCurrently, only v1 is available')

    mGroup.addParam('nRecycles', params.IntParam, label='Number of recycles: ', default=4,
                    help='Number of recycles to run. Defaults to number used in training (4)')
    mGroup.addParam('maxTokens', params.IntParam, label='Maximum number of tokens per GPU forward-pass: ', default=4,
                    expertLevel=params.LEVEL_ADVANCED,
                    help='Maximum number of tokens per gpu forward-pass. This will group shorter sequences together '
                         'for batched prediction. Lowering this can help with o')
    mGroup.addParam('chunkSize', params.IntParam, label='Chunk size: ', default=64, expertLevel=params.LEVEL_ADVANCED,
                    help='Chunks axial attention computation to reduce memory usage from O(L^2) to O(L). '
                         'Equivalent to running a for loop over chunks of of each dimension. '
                         'Lower values will result in lower memory usage at the cost of speed.')


  def _insertAllSteps(self):
    self._insertFunctionStep(self.predictStep)
    self._insertFunctionStep(self.createOutputStep)

  def predictStep(self):
    sequence = self.getInputSequence()
    model = self.getEnumText('modelName')
    seqName = self.getInputName()
    cwd = os.path.join(esmPlugin.getVar(ESM_DIC['home']), 'esm')

    args = f' -i {sequence} -m {model} -o {seqName} -od {os.path.abspath(self._getPath())}' \
           f' -g {self.gpuList.get().split(",")[0]}' \
           f' -cs {self.chunkSize.get()} -nr {self.chunkSize.get()} -mt {self.maxTokens.get()}'
    esmPlugin.runScript(self, scriptName, args, envDict=ESM_DIC, cwd=cwd)

  def createOutputStep(self):
    fnOut = self._getPath(f'{self.getInputName()}.pdb')

    if os.path.exists(fnOut):
      outStructFileName = self._getPath('outputStructureESMFold.cif')
      # Write conservation in a section of the output cif file
      ASH = AtomicStructHandler()

      esmDic = self.getESMFoldScoreDic()
      inpAS = toCIF(fnOut, self._getTmpPath('inputStruct.cif'))
      cifDic = ASH.readLowLevel(inpAS)
      cifDic = addScipionAttribute(cifDic, esmDic, self._ATTRNAME)
      ASH._writeLowLevel(outStructFileName, cifDic)

      outAS = AtomStruct(filename=outStructFileName)
      self._defineOutputs(outputStructure=outAS)


  def getESMFoldScoreDic(self):
    fnOut = self._getPath(f'{self.getInputName()}.pdb')
    ASH = AtomicStructHandler()
    ASH.read(fnOut)

    esmDic = {}
    for model in ASH.structure:
      for residue in model.get_residues():
        fId = residue.get_resname()
        chainName, resNumber = fId[2], fId[3][1]
        resId = '{}:{}'.format(chainName, resNumber)
        fscq_atom = residue.get_bfactor()
        esmDic[resId] = fscq_atom

    return esmDic

  def getInputSequence(self):
    return self.inputSequence.get().getSequence()

  def getInputName(self):
    return self.inputSequence.get().getSeqName()