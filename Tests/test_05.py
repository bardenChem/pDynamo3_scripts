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

	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_03","7tim_optLFBGS.pkl"),
		"set_qc_region":"yes",
		"residue_patterns":["*:LIG.248:*","*:GLU.164:*","*:HIE.94:*"],
		"set_energy_model":"QM",
		"Hamiltonian":"am1",
		"method_class":"SMO",
		"QCcharge":1,
	}


	if not os.path.exists( os.path.join("test_03","7tim_optLFBGS.pkl") ):
		try: os.system("python3 test_03.py")
		except: 
			print("There is no input file for this example! Run example #03!")
			return(False)

	test_01 = Scripts("test_05")
	#test_01.Set_System(system_parameters)
	#test_01.SaveSystem("7timAM1_qcMM.pkl")	


	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_03","7tim_optLFBGS.pkl"),
		"set_qc_region":"yes",
		"residue_patterns":["*:LIG.248:*","*:GLU.164:*","*:HIE.94:*"],
		"set_energy_model":"QM",
		"Hamiltonian":"b3lyp",
		"basis":"sto-3g",
		"method_class":"pySCF",
		"QCcharge":1,
	}

	test_02 = Scripts("test_05")
	test_02.Set_System(system_parameters)
	test_02.SaveSystem("7tim_PySCF_qcMM.pkl")	
	
#===================================
if __name__ == '__main__': Run_Test()