#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = pDynamo_Scripts.py

#--------------------------------------------------------------
import os, glob, sys

from commonFunctions import *
from SimulationSystem import SimulationSystem
from Simulation import *
from Analysis import *

#==============================================================
class Scripts:
	'''
	Set scripts and control the library of the presetted simulations
	'''

	def __init__(self, _projectFolder=None):
		'''
		Initialize project folder and paths.
		'''
		self.activeSystem    = None 
		self.system_historic = []
		self.projectFolder   = os.getcwd()

		if not _projectFolder: self.projectFolder = os.getcwd()
		else: self.projectFolder = os.path.join( os.getcwd(), _projectFolder )		
		if not os.path.exists(self.projectFolder): os.makedirs(self.projectFolder)
		

	#-----------------------------------------
	@classmethod
	def From_Input(selfClass,_inputFile):
		'''
		Reads input and automate system setting, simulation and analysis.
		'''
		
		RUN_SIMULATIONS = 0
		RUN_ANALYSIS    = 0 
		SET_CRD_NMB     = 0

		Spherical_Pruning = False
		Fixed_atoms = False
		reactions_crds = []
		_parameters = {}

		inpFile = open(_inputFile,"r")
		for line in lines:
			lines = line.split()
			if   lines[0] == "#INPUT_TYPE":
				_parameters["Input_Type"] = lines[1]
				if lines[1] == "geometry":
					_parameters["crd_file"] = lines[2]
				elif lines[1] == "amber" or lines[1] == "gromacs":
					_parameters["top_file"] = lines[2]
					_parameters["crd_file"] = lines[3]
				elif lines[1] == "pkl_file":
					_parameters["pkl_file"] = lines[2]
				elif lines[1] == "protein":
					_parameters["pdb_file"] = lines[2]
			elif lines[0] == "#SPHERICAL_PRUNNING":
				_parameters["spherical_prune"] = lines[1] # put the patterns to center the selection
				_parameters["spherical_prune_radius"] = float(lines[2])
			elif lines[0] == "#FIXED_ATOMS":
				_parameters["set_fixed_atoms"] = lines[1] # put the patterns to centert the selection
				_parameters["free_atoms_radius"] = float(lines[2])
			elif lines[0] == "#SET_REACTION_CRD":
				SET_CRD_NMB+=1
				_parameters["type_rc"+str(SET_CRD_NMB)]     = lines[1]
				_parameters["mass_ctr_rc"+str(SET_CRD_NMB)] = lines[2]
				for i in range(0, int(lines[3])):
					_parameters["atoms_rc"+str(SET_CRD_NMB)] = lines[i+4]
			elif lines[0] == "#SET_INITIAL_CRD":
				_parameters["set_initial_crd"] = lines[1]
			elif lines[0] == "#SET_QMMM_REGION":
				_parameters["set_qc_region"] = "yes"
				if lines[1] == "#residues_paterns":
					list_res = []
					for i in range( 0,int(lines[2]) ):
						list_res.append(lines[i+3]) 
					_parameters["residue_patterns"] = list_res
				elif lines[1] == "#center_atom":
					_parameters["center_atom"] = lines[2]
					_parameters["radius"] = float(lines[3])
			elif lines[0] == "#SET_ENERGY_MODEL_QM":
				_parameters["method_class"] = lines[1]
				


		_parameters["set_reaction_crd"] = SET_CRD_NMB
			
		self.Set_System(_parameters)

	#-----------------------------------------
	def Set_System(self,_parameters):
		'''
		Sets up the simulation system based on the provided parameters.
		'''
		_system4load =  SimulationSystem()

		mass_constraints = []

		if "mass_constraints" in _parameters:
			for mc in _parameters["mass_constraints"]:
				if mc == "yes": mass_constraints.append(True)
				elif mc == "no": mass_constraints.append(False)

		_debug_ok = False
		if "DEBUG" in _parameters: _debug_ok = True 
		if "Input_Type" not in _parameters:	raise KeyError("Missing required parameter: Input_Type")

		#Load system based on Input_Type using a dictionary for clarity
		input_methods = {
			"geometry": SimulationSystem.From_Coordinates,
			"amber": SimulationSystem.From_AMBER,
			"gromacs": SimulationSystem.From_Gromacs,
			"pkl": SimulationSystem.From_PKL,
			"protein":SimulationSystem.Protein_From_Coordinates,
		}

		try:
			input_type = _parameters["Input_Type"]
			load_function = input_methods[input_type]
			if input_type == "geometry":
				_system4load = load_function(_parameters["crd_file"])
			elif input_type == "amber":
				_system4load = load_function(_parameters["top_file"], _parameters["crd_file"])
				print(_parameters["top_file"],_parameters["crd_file"])
			elif input_type == "gromacs":
				_system4load = load_function(_parameters["top_file"], _parameters["crd_file"])
			elif input_type == "pkl":
				_system4load = load_function(_parameters["pkl_file"])
			elif input_type == "protein":
				_system4load = load_function(_parameters["pdb_file"])
			else: raise ValueError(f"Unsupported Input_Type: {input_type}")

		except (KeyError, TypeError) as e:
		# Handle potential issues with missing keys or invalid values
			raise ValueError(f"Error loading system: {e}") from e  # Chain the original exception
	
		self.activeSystem = _system4load
		#------------------------------------
		if "spherical_prune" in _parameters:			
			self.activeSystem.Spherical_Pruning(_parameters["spherical_prune"],float(_parameters["spherical_prune_radius"]))
		if "set_fixed_atoms" in _parameters:
			self.activeSystem.Setting_Free_Atoms(_parameters["set_fixed_atoms"],float(_parameters["free_atoms_radius"]))
		if "set_reaction_crd" in _parameters:
			for rc in range(0,_parameters["set_reaction_crd"]):
				self.activeSystem.Set_Reaction_crd( _parameters["atoms_rc"+str(rc+1)],_parameters["type_"+str(rc+1)],mass_constraints[rc])
		if "set_initial_crd" in _parameters:
			if ( _parameters["set_initial_crd"][-4:] ) == ".pkl":
				try:				
					self.activeSystem.system.coordinates3 = Unpickle(_parameters["set_initial_crd"])[0]
				except:
					self.activeSystem.system.coordinates3 = Unpickle(_parameters["set_initial_crd"])
			else: self.activeSystem.system.coordinates3 = ImportCoordinates3(_parameters["set_initial_crd"])
		if "set_qc_region" in _parameters:
			_residue_list = []
			_centerAtom = None
			_radius = None
			if "residue_patterns" in _parameters: _residue_list = _parameters["residue_patterns"]
			if "center_atom"      in _parameters: _centerAtom   = _parameters["center_atom"]
			if "radius"           in _parameters: _radius       = _parameters["radius"]
			if _parameters["set_qc_region"] == "yes":				
				self.activeSystem.Set_QCMM_Region(_residue_list,_centerAtom,_radius)
		if "set_energy_model" in _parameters:
			if _parameters["set_energy_model"] == "QM":
				self.activeSystem.Set_QC_Method(_parameters,_DEBUG=_debug_ok)

		#-----------------------------------
		self.system_historic.append(_system4load)
		self.activeSystem.Check()


	#-----------------------------------------
	def Run_Simulation(self,_parameters):
		'''
		'''
		_parameters["active_system"] = self.activeSystem
		_parameters["project_folder"]= self.projectFolder
		self.system_historic.append(self.activeSystem) 
		#----------------------------------------------- 
		_Run_Simulation   = Simulation(_parameters)
		self.activeSystem.system = _Run_Simulation.Execute()

	#-----------------------------------------
	def Run_Analysis(self,_parameters):
		'''
		'''
		_parameters["active_system"] = self.activeSystem
		self.system_historic.append(self.activeSystem) 
		#--------------------------------------------
		_Analysis = Analysis(_parameters)
		_Analysis.Execute()
		
	#-----------------------------------------
	def PrintSystems(self):
		'''
		Method to print the summary of the loaded systems 
		'''
		print("There are {} loaded systems".format( len(self.system_historic) ) )
		ctn = input("Type any key to print the Summary of the Systems, or 'N' to cancel this")
		if not ctn == "N":
			if len(self.system_historic) > 0:
				for system in self.system_historic:
					system.system.Summary()
					print("***************************************************")
				print("Now, printing the current system Summary:")
				self.activeSystem.system.Summary()
		
			elif len(self.system_historic) == 1: print("There is only the current System loaded!\n Printing its information below!")
			else:                                print( "There are no loaded systems!")

	#.-------------------------------------------------------------------------
	def SaveProject(self):
		'''
		The complete version of this function intends to save in pkl and another coordinate format
		the systems and trajectories worked in this simulations
		Though, in the current state only will save the current system to a pkl file
		'''
		_baseName = self.activeSystem.basename
	
		savePathPkl = os.path.join(self.projectFolder,_baseName+".pkl")
		savePathPdb = os.path.join(self.projectFolder,_baseName+".pdb")
	
		for _System in self.system_historic:
			_label = _System.label
			savePathPkl = os.path.join(savePathPkl, _baseName+"_"+_label+".pkl")
			Pickle( savePathPkl, _System )        
    #.-------------------------------------------------------------------------
	def SaveSystem( self, _cname=None):
		'''
		'''
	
		_baseName = self.activeSystem.baseName
		if _cname:
			savePathPkl = os.path.join(self.projectFolder,_cname+".pkl")
			savePathPdb = os.path.join(self.projectFolder,_cname+".pdb")
		else:
			savePathPkl = os.path.join(self.projectFolder,_baseName+".pkl")
			savePathPdb = os.path.join(self.projectFolder,_baseName+".pdb")
	
			i = 0;
			while os.path.exists(savePathPdb):
				i += 1
				savePathPdb = os.path.join(self.projectFolder,_baseName+".pdb")
				savePathPdb = savePathPdb[:-4] + "_#{}.pdb".format(i)
			while os.path.exists(savePathPkl):
				i += 1
				savePathPkl = os.path.join(self.projectFolder,_baseName+".pkl")
				savePathPkl = savePathPkl[:-4] + "_#{}.pkl".format(i)
		#----------------------------------------------------------------
		self.activeSystem.system.Summary()
		ExportSystem( savePathPkl,self.activeSystem.system )
		ExportSystem( savePathPdb,self.activeSystem.system )
#==============================================================

if __name__=="__main__":
	'''
	'''
	pass
