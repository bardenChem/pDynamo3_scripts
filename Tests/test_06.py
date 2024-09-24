#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os, sys
#==============================================

#-------------------------------------------------
def Simple_Distance2D(_hamiltonian):
	'''
	'''
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_opt"+_hamiltonian,"7tim_"+_hamiltonian+"_opt_PF.pkl"),		
		"set_reaction_crd":2,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02"],
		"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2"],
		"type":"Distance",
		"mass_constraint":"True",
	}		
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_rc1":0.1,
		"nsteps_rc1":16,
		"dincre_rc2":0.1,
		"nsteps_rc2":16,
		"maxIterations":2200,
		"optmizer":"SteepestDescent",
		"NmaxThreads":8,
		"force_constants":[4000.0,4000.0]
	}
	#test simple distance
	test_01 = Scripts("test_06/Simple_Distance2D_"+_hamiltonian)
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)
	test_01.SaveSystem("Simple_DistanceScan2D")


#-------------------------------------------------
def Mixed_Distance2D(_hamiltonian):
	'''
	'''
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_opt"+_hamiltonian,"7tim_"+_hamiltonian+"_opt_PF.pkl"),		
		"set_reaction_crd":2,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2"],
		"type":"Distance",
		"maxIterations":2200,
		"mass_constraint":"True",
	}		
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_rc1":0.15,
		"dincre_rc2":0.15,
		"optmizer":"SteepestDescent",
		"maxIterations":2200,
		"nsteps_rc1":16,
		"nsteps_rc2":16,
		"NmaxThreads":8,
		"force_constants":[4000.0,4000.0]
	}
	#test simple distance
	test_01 = Scripts("test_06/Mixed_Distance_"+_hamiltonian)
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)
	test_01.SaveSystem("Mixed_DistanceScan")


#-----------------------------------------------
def Multiple_Distance2D(_hamiltonian):
	'''
	'''
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_opt"+_hamiltonian,"7tim_"+_hamiltonian+"_opt_PF.pkl"),		
		"set_reaction_crd":2,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:HIE.94:NE2","*:HIE.94:HE2","*:LIG.*:O06"],
		"type":"multipleDistance",
		"maxIterations":2200,
		"mass_constraint":"True",
	}		
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_rc1":0.1,
		"dincre_rc2":0.1,
		"optmizer":"ConjugatedGradient",
		"maxIterations":2200,
		"nsteps_rc1":16,
		"nsteps_rc2":16,
		"NmaxThreads":16,
		"force_constants":[2000.0,2000.0]

	}
	#test simple distance
	test_01 = Scripts("test_06/Multiple_Distance_"+_hamiltonian)
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)
	test_01.SaveSystem("Multiple_DistanceScan")



#-----------------------------------------------
def Run_Test():
	'''
	Test 
	'''

	if not os.path.exists( os.path.join("test_05","7tim.pkl") ):
		Prepare_MM_System()
	
	if not os.path.exists( os.path.join("test_05","7tim_optMM.pkl") ):
		Prepare_Prune_System()	

	#Simple_Distance2D("")	
	#Mixed_Distance2D("pm6")
	_methods = ["am1dphot","pddgpm3"]

	for smo in _methods:
		Multiple_Distance2D(smo)
		
	
#===================================
if __name__ == '__main__': Run_Test()