#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = Analysis.py

#--------------------------------------------------------------
import os, glob, sys

from commonFunctions import *
from Simulation import Simulation
from Analysis import *

#==============================================================
class Interface:
	'''
	'''

	def __init__(self,_projectFolder = None):
		'''
		'''
		self.activeSystem    = None 
		self.system_historic = []
		self.projectFolder   = os.getcwd()

		if _projectFolder: self.projectFolder = os.getcwd()		
        if not os.path.exists(self.projectFolder): os.makedirs(self.projectFolder)


	#-----------------------------------------
	@classmethod
	def From_Input(selfClass,_inputFile):
		'''
		'''
		pass

	#-----------------------------------------
	def Set_System(self,_parameters):
		'''
		'''
		pass

	#-----------------------------------------
	def Run_Simulation(self,_parameters):
		'''
		'''
		_parameters["active_system"] = self.activeSystem
		self.system_historic.append(self.activeSystem) 
		#----------------------------------------------- 
		_Run_Simulation   = Simulation(_parameters)
        self.activeSystem = _Run_Simulation.Execute()		 

	#-----------------------------------------
	def Run_Individual_Analysis(self,_parameters):
		'''
		'''
		_parameters["active_system"] = self.activeSystem
		self.system_historic.append(self.activeSystem) 
		#--------------------------------------------
		_Analysis = Analysis(_parameters)
		self.activeSystem = _Analysis.Execute()
		


#==============================================================

if __name__=="__main__":
	'''
	'''
	pass