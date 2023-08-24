 #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = Analysis.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#--------------------------------------------------------------
import os, glob, sys
import numpy as np
#--------------------
from pMolecule import *
#==============================================================
class QuantumMethods:
	'''
	Classe to set up quantum chemical method to the system.
	'''	
	def __init__(self):
		'''
		Default constructor
		'''
		self.methodClass 	= "SMO" #  SMO, HF, DFT, ORCA, DFTB, MOPAC and PYSCF
		self.selectionType  = "Hybrid" # "Hybrid", "Whole" 
		self.selection      = None 
		self.systemBase     = None 
		self.QuantumSystem  = None 
		self.converger      = "standard" # veryLoose, loose, standard, tight, veryTight 

	#-------------------------------------
	@classmethod
	def From_Parameters(selfClass,_parameters):
		'''
		Create object setting the quantum chemical method.
		Returns the system with quantum chemical enabled.
		'''
		if _parameters == ""


		pass 

	#-------------------------------------
	def Set_Converger(self):
		'''
		'''
		pass 

	#-------------------------------------
	def Change_ConvergerPars(self,_parameters):
		'''
		'''
		pass 
	#------------------------------------- 
	def Change_QC_Region(self,_parameters):
		'''
		'''
		pass 
#=================================================
