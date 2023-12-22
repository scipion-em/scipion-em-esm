# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
# *  e-mail address 'ddelhoyo@cnb.csic.es'
# *
# **************************************************************************

import os, argparse
import torch, esm

if __name__ == "__main__":
    '''Use: python <scriptName> -i/--inputFilename <sysFile> -o/--outputName <outputName> 
    -t/--inputTraj [<trjFile>]
    If only an MD system file is provided, then the script will convert the system to the required format,
    based on outputName.
    If a trajectory file is provided too, then the script will convert the trajectory to the required format,
    based on outputName, using the system file to load it.
    '''
    parser = argparse.ArgumentParser(description='Handles the IO for molecule files using openbabel')
    parser.add_argument('-i', '--inputSequence', type=str, help='Input sequence')
    parser.add_argument('-o', '--outputName', type=str, help='Output name')
    parser.add_argument('-od', '--outputDir', type=str, help='Output directory')
    parser.add_argument('-m', '--ESMModel', type=str, default='esmfold_v1', help='ESMFold model to use')

    parser.add_argument('-nr', '--numberRecycles', type=int, default=4, help='Number of recycles')
    parser.add_argument('-cs', '--chunkSize', type=int, default=128, help='Chunk size to use the model')
    parser.add_argument('-g', '--gpuId', type=int, default=0, help='GPU index to use')

    args = parser.parse_args()
    sequence, modelName = args.inputSequence, args.ESMModel
    outFile = os.path.join(args.outputDir, args.outputName) + '.pdb'

    model = getattr(esm.pretrained, modelName)()
    model = model.to(torch.device(f'cuda:{args.gpuId}'))
    model.set_chunk_size(args.chunkSize)
    model.eval()

    with torch.no_grad():
      output = model.infer_pdb(sequence, num_recycles=args.numberRecycles)

    with open(outFile, "w") as f:
      f.write(output)


