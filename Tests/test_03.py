#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os, sys
#===================================
def Run_Test():
	'''
	Test geometry optimization algorithms
	'''

	algs = ["ConjugatedGradient",
			"LFBGS"             ,
			"SteepestDescent"   ,
			"FIRE"              ]
	
	if not os.path.exists( os.path.join("test_01","7tim.pkl") ):
		try: os.system("python3 test_01.py")
		except: 
			print("There is no input file for this example! Run example #01!")
			return(False)

	_parameters = {
		"Input_Type":"pkl",
		"pkl_file":"test_01/7tim.pkl",
		"simulation_type":"Geometry_Optimization",
		"save_format":".dcd",
		"save_frequency":20,
		"rmsGradient":0.1,
		"maxIterations":1200,
		"folder":"test_03"
	}
	#------------------------------------
	test_01 = Scripts("test_03")
	for alg in algs:
		test_01.Set_System(_parameters)
		_parameters["optmizer"]=alg
		_parameters["trajectory_name"]="7timMMopt_"+alg+".ptGeo"
		test_01.Run_Simulation(_parameters)
		test_01.SaveSystem("7tim_opt"+alg)
	#-----------------------------------
	
	
#===================================
if __name__ == '__main__': Run_Test()