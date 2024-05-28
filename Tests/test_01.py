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
	#test load force field topology and coordinate files 
	




#===================================
if __name__ == '__main__': Run_Test()