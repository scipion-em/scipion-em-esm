# **************************************************************************
# *
# * Authors:	Carlos Oscar Sorzano (coss@cnb.csic.es)
# *			 	Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *			 	Martín Salinas Antón (martin.salinas@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains protocols for creating and using ConPLex models for virtual screening
"""

# General imports
import os, subprocess, json

# Scipion em imports
import pwem
from scipion.install.funcs import InstallHelper

# Plugin imports
from .bibtex import _bibtexStr
from .constants import *

# Pluging variables
_logo = 'mit_logo.png'

class Plugin(pwem.Plugin):
	"""
	"""
	_dfdHome = os.path.join(pwem.Config.EM_ROOT, ESM_DIC['name'] + '-' + ESM_DIC['version'])

	@classmethod
	def _defineVariables(cls):
		cls._defineEmVar(ESM_DIC['home'], cls._dfdHome)

	@classmethod
	def defineBinaries(cls, env):
		"""
        This function defines the binaries for each package.
        """
		cls.addEsmPackage(env)

	@classmethod
	def addEsmPackage(cls, env, default=True):
		""" This function provides the neccessary commands for installing AutoDock. """
		# Instantiating the install helper
		installer = InstallHelper(ESM_DIC['name'], packageHome=cls.getVar(ESM_DIC['home']),
															packageVersion=ESM_DIC['version'])

		# Installing package
		installer.getCloneCommand('https://github.com/facebookresearch/esm.git', targeName='ESM_CLONED') \
			.addCommand(f'cd esm && conda env create -f environment.yml && conda rename -n esmfold {cls.getEnvName(ESM_DIC)} && '
									f'{cls.getEnvActivationCommand(ESM_DIC)} && conda install -y einops=0.6.1 && '
									f'python setup.py build && python setup.py install && '
									f'pip install "fair-esm[esmfold]" && pip install "dllogger @ git+https://github.com/NVIDIA/dllogger.git" && '
									f'pip install "openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307" ',
									'ESMFOLD_INSTALLED')\
			.addPackage(env, ['git', 'conda', 'pip'], default=default)


	# ---------------------------------- Protocol functions-----------------------
	@classmethod
	def getPluginHome(cls, path=""):
		import esm
		fnDir = os.path.split(esm.__file__)[0]
		return os.path.join(fnDir, path)

	@classmethod
	def getEnvName(cls, packageDictionary):
		""" This function returns the name of the conda enviroment for a given package. """
		return '{}-{}'.format(packageDictionary['name'], packageDictionary['version'])

	@classmethod
	def getEnvActivationCommand(cls, packageDictionary, condaHook=True):
		""" This function returns the conda enviroment activation command for a given package. """
		return '{}conda activate {}'.format(cls.getCondaActivationCmd() if condaHook else '',
																				cls.getEnvName(packageDictionary))

	@classmethod
	def getScriptsDir(cls, scriptName):
		return cls.getPluginHome('scripts/%s' % scriptName)

	@classmethod
	def runScript(cls, protocol, scriptName, args, envDict, cwd=None, isSubprocess=False):
		""" Run ESM script from a given protocol. """
		scriptName = cls.getScriptsDir(scriptName)
		fullProgram = '%s && %s %s' % (cls.getEnvActivationCommand(envDict), 'python', scriptName)
		if not isSubprocess:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)
