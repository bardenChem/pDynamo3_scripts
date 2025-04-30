#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os
#===================================
def Run_Test():
	'''
	Test the setting of all quantum methods classes available.
	Internal semi-empirical methods in pDynamo
	DFTB+ 
	Orca 
	Mopac
	pySCF 
	'''
	_parameters = {
		"Input_Type":"geometry",
		"crd_file":os.path.join("data","cyclohexane_single_frame.xyz"),
		"set_energy_model":"QM",
		"functional":"HF",
		"method_class":"ORCA",
		"basis":"6-31G*",
		"save_frequency":5,
		"scratch":"test_02_orca_cyclohex",
		"NmaxThreads":1
	}
	
	test_02 = Scripts("test_02_orca_cyclohex")
	test_02.Set_System(_parameters)	
	_parameters["simulation_type"] = "Geometry_Optimization"
	#test_02.Run_Simulation(_parameters)
	test_02.SaveSystem()
	#-------------------------------
	#test QC/MM from gromacs
	_parameters["Input_Type"]      ="pkl"
	_parameters["pkl_file"]        ="test_01/1atp_peptide.pkl"
	_parameters["set_qc_region"]   ="yes"
	_parameters["residue_patterns"]=["*:ARG.19:*"]
	_parameters["QCcharge"]        = 1
	_parameters["scratch"]		   = "test_02_orca_qmmm_1atp"

	test_03 = Scripts("test_02_orca_qmmm_1atp")
	test_03.Set_System(_parameters)
	test_03.SaveSystem()

	#test QC/MM from AMBER
	_parameters["residue_patterns"]= ["*:LIG.248:*","*:GLU.164:*","*:HIE.94:*"]
	_parameters["pkl_file"]        = "test_01/7tim.pkl"
	_parameters["scratch"]		   = "test_02_orca_qmmm_tim"
	_parameters["QCcharge"]        = -3

	test_04 = Scripts("test_02_orca_qmmm_tim")
	test_04.Set_System(_parameters)
	test_04.SaveSystem()





	
#===================================
if __name__ == '__main__': Run_Test()