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
  """
  The ProtESMFoldPrediction class implements a structural prediction
  protocol based on the ESMFold deep learning model for protein structure
  inference from amino acid sequences. The protocol is designed to predict
  three-dimensional atomic models directly from an input protein sequence,
  providing an automated workflow for structural modeling within Scipion
  frameworks.

  The protocol operates by receiving a protein sequence as input together
  with a selected ESMFold model configuration. Users can define execution
  parameters such as the number of recycles used during prediction,
  chunk size for memory optimization, and GPU selection for accelerated
  computation. The implementation supports parallel execution and is
  optimized for handling large protein sequences efficiently through axial
  attention chunking strategies that reduce memory consumption.

  During execution, the protocol constructs the required command-line
  arguments and launches the ESMFold prediction script within the
  configured environment. The prediction step generates a structural model
  in PDB format corresponding to the inferred three-dimensional
  conformation of the input sequence.

  After prediction, the protocol converts the generated structure into a
  CIF representation and enriches it with confidence-related information
  extracted from the model output. The implementation reads atom-level
  B-factor values from the predicted structure and interprets them as
  ESMFold confidence scores. These scores are stored as additional
  Scipion-compatible attributes associated with the atomic coordinates,
  allowing downstream protocols and visualization tools to access residue
  confidence information directly.

  The protocol finally produces an output atomic structure object
  containing the predicted model and its associated confidence metadata.
  This output can subsequently be used for structural visualization,
  comparative analysis, fitting into cryo-EM maps, or additional
  computational modeling workflows.

  Overall, ProtESMFoldPrediction provides an integrated framework for
  protein structure prediction using ESMFold models, combining automated
  sequence-based inference, GPU-accelerated execution, memory-efficient
  processing, and confidence score integration into a unified structural
  modeling workflow.
  """

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
           f' -cs {self.chunkSize.get()} -nr {self.nRecycles.get()}'
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
      cifDic = addScipionAttribute(cifDic, esmDic, self._ATTRNAME, recipient='atoms')
      ASH._writeLowLevel(outStructFileName, cifDic)

      outAS = AtomStruct(filename=outStructFileName)
      self._defineOutputs(outputStructure=outAS)


  def getESMFoldScoreDic(self):
    fnOut = self._getPath(f'{self.getInputName()}.pdb')
    ASH = AtomicStructHandler()
    ASH.read(fnOut)

    esmDic = {}
    for model in ASH.structure:
      for atom in model.get_atoms():
        fId = atom.get_full_id()
        chainName, resNumber, atomName = fId[2], fId[3][1], fId[4][0]
        atomId = '{}:{}@{}'.format(chainName, resNumber, atomName)
        esmScore = atom.get_bfactor()
        esmDic[atomId] = esmScore

    return esmDic

  def getInputSequence(self):
    return self.inputSequence.get().getSequence()

  def getInputName(self):
    return self.inputSequence.get().getSeqName()