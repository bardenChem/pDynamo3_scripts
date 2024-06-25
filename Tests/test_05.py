#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os, sys
#===================================
def Run_Test():
	'''
	Test 
	'''

	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_03","7tim_optLFBGS.pkl"),
		"set_qc_region":"yes",
		"residue_patterns":["*:LIG.248:*","*:GLU.164:*","*:HIE.94:*"],
		"set_energy_model":"QM",
		"Hamiltonian":"am1",
		"method_class":"SMO",
		"QCcharge":1,
		"set_reaction_crd":2,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		#"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2","*:HIE.94:NE2"],
		"mass_constraint":"true",
		"type":"distance"	
	}

	optmization_parameters = {
		"simulation_type":"Geometry_Optimization",
		"save_format":".dcd",
		"save_frequency":20,
		"rmsGradient":0.1,
		"maxIterations":1200,
		"folder":"test_05"
	}
	
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_RC1":-0.1,
		"nsteps_RC!":12
	}

	if not os.path.exists( os.path.join("test_03","7tim_optLFBGS.pkl") ):
		try: os.system("python3 test_03.py")
		except: 
			print("There is no input file for this example! Run example #03!")
			return(False)

	test_01 = Scripts("test_05")
	test_01.Set_System(system_parameters)
	test_01.SaveSystem("7timAM1_qcMM.pkl")
	test_01.Run_Simulation(optmization_parameters)	
	test_01.Run_Simulation(scan1_parameters)





	
	
#===================================
if __name__ == '__main__': Run_Test()