#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os, sys
#===================================
def Run_Test():
	'''
	Test umbrella sampling 
	'''

	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_optam1","7tim_am1_opt_PF.pkl"),		
		"set_reaction_crd":1,	
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		#"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2","*:HIE.94:NE2"],
		"type":"Distance",
		"mass_constraint":"True",
	}

	_path   = "test_05/Multiple_Distance_rm1/ScanTraj.ptGeo"
	
	US_parameters = {
				   "ndim": 1 					        ,
				   "sampling_production":100	   		,
				   "equilibration_nsteps":2000			,
				   "production_nsteps":20000   			,
				   "source_folder":_path 		       	,
				   "MD_method":"LeapFrog"		      	,
				   "MC_RC1":True				        ,
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":10		  				,
				   "xnbins":10                          ,
				   "xwindows":20                        ,  
				   "trajectory_name":"US_test"          ,
				   }
	#------------------------------------
	test_01 = Scripts("test_11")
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(US_parameters)
	test_01.SaveSystem()
	#-----------------------------------	
	
	
#===================================
if __name__ == '__main__': Run_Test()