#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = Analysis.py

#--------------------------------------------------------------
import os, glob, sys

from commonFunctions import *
from Simulation import Simulation
from Analysis import *

#==============================================================
class Scripts:
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
		
	#-----------------------------------------
	def PrintSystems(self):
        '''
        Method to print the summary of the loaded systems 
        '''
        print("There are {} loaded systems".format( len(self.system_historic) ) )
        ctn = input("Type any key to print the Summary of the Systems, or 'N' to cancel this")
        if not ctn == "N":
            if len(self.system_historic) > 0:
                for system in self.system_historic:
                    system.system.Summary()
                    print("***************************************************")
                print("Now, printing the current system Summary:")
                self.activeSystem.system.Summary()

            elif len(self.system_historic) == 1: print("There is only the current System loaded!\n Printing its information below!")
            else:                                print( "There are no loaded systems!")

	#.-------------------------------------------------------------------------
    def SaveProject(self):
        '''
        The complete version of this function intends to save in pkl and another coordinate format
        the systems and trajectories worked in this simulations
        Though, in the current state only will save the current system to a pkl file
        '''
        _baseName = self.activeSystem.basename

        savePathPkl = os.path.join(self.projectFolder,_baseName+".pkl")
        savePathPdb = os.path.join(self.projectFolder,_baseName+".pdb")

        for _System in self.system_historic:
        	_label = _System.label
            savePathPkl = os.path.join(savePathPkl, _baseName+"_"+_label+".pkl")
            Pickle( savePathPkl, _System )        
    #.-------------------------------------------------------------------------
    def SaveSystem( self, _cname=None):
        '''
        '''

        _baseName = self.activeSystem.basename
        if _cname:
            savePathPkl = os.path.join(self.projectFolder,_cname+".pkl")
            savePathPdb = os.path.join(self.projectFolder,_cname+".pdb")
        else:
            savePathPkl = os.path.join(self.projectFolder,_baseName+".pkl")
            savePathPdb = os.path.join(self.projectFolder,_baseName+".pdb")

            i = 0;
            while os.path.exists(savePathPdb):
                i += 1
                savePathPdb = os.path.join(self.projectFolder,_baseName+".pdb")
                savePathPdb = savePathPdb[:-4] + "_#{}.pdb".format(i)
            while os.path.exists(savePathPkl):
                i += 1
                savePathPkl = os.path.join(self.projectFolder,_baseName+".pkl")
                savePathPkl = savePathPkl[:-4] + "_#{}.pkl".format(i)
        #----------------------------------------------------------------
        Pickle( savePathPkl,self.activeSystem.system )
        ExportSystem( savePathPdb,self.activeSyste.system )
#==============================================================

if __name__=="__main__":
	'''
	'''
	pass