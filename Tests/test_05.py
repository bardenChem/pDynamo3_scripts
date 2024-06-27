#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os, sys
#==============================================
def Prepare_MM_System():
	'''
	'''
	_opt_pars = {
		"Input_Type":"amber",
		"simulation_type":"Geometry_Optimization",
		"save_format":".dcd",
		"save_frequency":20,
		"rmsGradient":1,
		"maxIterations":2200,
		"folder":"test_03",
		"trajectory_name":"opt_full_tim.ptGeo"
	}

	#preparation of 7tim structure	
	#test load amber force field topology and coordinate files 
	_opt_pars["crd_file"] = os.path.join("data","7tim.crd")
	_opt_pars["top_file"] = os.path.join("data","7tim.top")
	test_03 = Scripts("test_05")
	test_03.Set_System(_opt_pars)
	test_03.Run_Simulation(_opt_pars)
	test_03.SaveSystem()

#-----------------------------------------------
def Prepare_Prune_System():
	'''
	'''
	_prune_pars= {
		"Input_Type":"pkl",
		"pkl_file":"test_05/7tim.pkl",
		"simulation_type":"Geometry_Optimization",
		"spherical_prune":"*:LIG.248:C02",
		"spherical_prune_radius":25.0,
		"set_fixed_atoms":"*:LIG.248:C02",
		"free_atoms_radius":20.0,
		"save_format":".dcd",
		"save_frequency":20,
		"rmsGradient":0.1,
		"maxIterations":2200,
		"trajectory_name":"opt_pruned_tim.ptGeo"		
	}

	test_03_b = Scripts("test_05")
	test_03_b.Set_System(_prune_pars)
	test_03_b.Run_Simulation(_prune_pars)
	test_03_b.SaveSystem("7tim_optMM")


#-----------------------------------------------
def Set_QC_MM():
	'''
	'''
	_qc_mmpars = {
		"Input_Type":"pkl",
		"pkl_file":"test_05/7tim_optMM.pkl",
		"set_energy_model":"QM",
		"Hamiltonian":"am1",
		"method_class":"SMO",
		"set_qc_region":"yes",
		"residue_patterns":["*:LIG.248:*","*:GLU.164:*","*:HIE.94:*"],
		"QCcharge":1,
		"save_format":".dcd",
		"save_frequency":20,
		"rmsGradient":0.1,
		"maxIterations":2200,
		"trajectory_name":"opt_qcmm_tim.ptGeo",
		"simulation_type":"Geometry_Optimization"
	}

	test_03_c = Scripts("test_05")
	test_03_c.Set_System(_qc_mmpars)
	test_03_c.Run_Simulation(_qc_mmpars)
	test_03_c.SaveSystem("7tim_qcmm_opt_PF")

#-----------------------------------------------
def Run_Test():
	'''
	Test 
	'''

	if not os.path.exists( os.path.join("test_05","7tim.pkl") ):
		Prepare_MM_System()
	
	if not os.path.exists( os.path.join("test_05","7tim_optMM.pkl") ):
		Prepare_Prune_System()

	if not os.path.exists( os.path.join("test_05","7tim_qcmm_opt_PF.pkl") ):
		Set_QC_MM()	

	
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_03","7tim_qcmm_opt_PF.pkl"),		
		"set_reaction_crd":1,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02"],
		#"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		#"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2","*:HIE.94:NE2"],
		"type":"distance",
		"mass_constraint":"true",
	}		
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_RC1":0.1,
		"nsteps_RC1":12
	}
	
	#test simple distance
	test_01 = Scripts("test_05")
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)	
	
#===================================
if __name__ == '__main__': Run_Test()