#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os
#===================================
def Run_Test():
	'''
	'''
	_parameters = {
		"Input_Type":"geometry",
		"crd_file":os.path.join("data","cyclohexane_single_frame.xyz"),		
	}
	#test load xyz
	test_01 = Scripts()
	test_01.Set_System(_parameters)
	#test load gromacs topology and coordinate files 
	_parameters["Input_Type"] = "gromacs"
	_parameters["crd_file"] = os.path.join("data","crambin.gro")
	_parameters["top_file"] = os.path.join("data","crambin.top")	
	test_02 = Scripts()
	test_02.Set_System(_parameters)
	#test load amber force field topology and coordinate files 
	_parameters["Input_Type"] = "amber"
	_parameters["crd_file"] = os.path.join("data","7tim.crd")
	_parameters["top_file"] = os.path.join("data","7tim.top")
	test_03 = Scripts()
	test_03.Set_System(_parameters)
	#test load pkl and test spherical pruning and fixed atoms
	_parameters["Input_Type"] = "pkl"
	_parameters["pkl_file"] = "7tim.pkl"
	_parameters["spherical_prune"] = "*:LIG.248:C02"
	_parameters["spherical_prune_radius"] = 25.0
	_parameters["set_fixed_atoms"] = "*:LIG.248:C02"
	_parameters["free_atoms_radius"] = 20.0
	test_04 = Scripts()
	test_04.Set_System(_parameters)
#===================================
if __name__ == '__main__': Run_Test()