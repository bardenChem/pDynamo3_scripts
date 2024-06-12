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
	SMOmodels = ["am1","am1dphot","pddgpm3","pm3","rm1","pm6"]		
	_parameters = {
		"Input_Type":"geometry",
		"crd_file":os.path.join("data","cyclohexane_single_frame.xyz"),
		"set_energy_model":"QM",
		"Hamiltonian":"am1",
		"method_class":"SMO"
	}
	for smo in SMOmodels:
		_parameters["Hamiltonian"] = smo
		test_01 = Scripts()
		test_01.Set_System(_parameters)
	#------------------------------------
	'''
	_parameters["method_class"]="DFTB"
	_parameters["Hamiltonian"] ="DFTB"
	test_02 = Scripts()
	test_02.Set_System(_parameters)
	#-----------------------------------
	_parameters["method_class"]="orca"
	_parameters["functional"]  ="HF"
	_parameters["basis"]       ="6-31G*"
	test_03 = Scripts()
	test_03.Set_System(_parameters)	
	#---------------------------------
	'''
	_parameters["method_class"]="pySCF"
	_parameters["functional"]  ="b3lyp"
	_parameters["basis"]       ="6-31G*"
	test_04 = Scripts()
	test_04.Set_System(_parameters)
	#----------------------------------
	_parameters["method_class"]="abinitio"
	_parameters["basis"]       ="dgauss-dzvp"
	test_05 = Scripts()
	test_05.Set_System(_parameters)
	#----------------------------------



	
#===================================
if __name__ == '__main__': Run_Test()