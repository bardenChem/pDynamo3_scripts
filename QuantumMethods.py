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
		self.Hybrid  		= False  
		self.selection      = None 
		self.systemBase     = None 
		self.QuantumSystem  = None 
		self.convergerLevel = "standard" # veryLoose, loose, standard, tight, veryTight 
		self.converger      = None
		self.qcSystem       = None 
		self.system         = None
		self.qcModel        = None
	#==================================================
	@classmethod
	def From_Parameters(selfClass,_parameters):
		'''
		Create object setting the quantum chemical method.
		Returns the system with quantum chemical enabled.
		'''
		self.system = _parameters["active_system"]
		self.methodClass = _parameters["method_class"]
		if "region" in _parameters: self.selectionType = "Hybrid"

		#---------------------------------------------
		if self.Hybrid:
			atomlist = []
        	for sel in _parameters["region"]:
            	if type(sel) == int:
                	self.selection.append(sel)
            	elif type(_parameters["region"]) == list:
                	for i in range( len(sel) ):
                    	self.selection.append( sel[i] ) 
                self.selection = Selection.FromIterable(self.selection) 
                
       		NBmodel = self.system.nbModel
			self.system.nbModel = None
        	oldSystem = Clone(self.system)
       	#--------------------

    	self.Set_Converger()
        self.system.electronicState = ElectronicState.WithOptions( charge = _parameters["QCcharge"],multiplicity = _parameters["multiplicity"] )
        if self.methodClass == "SMO":
        	self.qcModel = QCModelMNDO.WithOptions( hamiltonian = _parameters["Hamiltonian"], converger=self.converger )
        elif self.methodClass == "abinitio":
        	self.qcModel = QCModelDFT.WithOptions(converger   = self.converger 			 ,
        										  functional  = _parameters["functional"],
        										  orbitalBasis= _parameters["basis"] 	 ,
        										  fitBasis    = _parameters["fit_basi"]  )

        #------------------------------------------------------------------------
        if self.hybrid: self.system.DefineQCModel( self.qcModel, qcSelection=_QCRegion )
    	else: self.system.DefineQCModel(self.qcModel)
        self.system.DefineNBModel( NBmodel )
        #------------------------------------------------------------------------
        energy = self.system.Energy()
        
	#-------------------------------------
	def Export_QC_System(self):
		'''
		'''		
		ExportSystem(self.baseName+"/qcSystem.pdb",qcSystem)
        ExportSystem(self.baseName+"/qcSystemEntire.pdb",self.system)

	def Set_Converger(self,_parameters=None):
		'''
		'''
		EnergyTolerance  = 3.0e-4
		DensityTolerance = 1.0e-8
		MaxIterations    = 500 

		if self.convergerLevel == "Very_Loose":
			EnergyTolerance  = 3.0e-4
			DensityTolerance = 1.0e-8
			MaxIterations    = 500
		elif  self.convergerLevel == "Loose":
			EnergyTolerance  = 3.0e-4
			DensityTolerance = 1.0e-8
			MaxIterations    = 500
		elif self.convergerLevel == "standard":
			EnergyTolerance  = 3.0e-4
			DensityTolerance = 1.0e-8
			MaxIterations    = 500
		elif self.convergerLevel == "Tight":
			EnergyTolerance  = 3.0e-4
			DensityTolerance = 1.0e-8
			MaxIterations    = 500
		elif self.convergerLevel == "Very_Tight":
			EnergyTolerance  = 3.0e-4
			DensityTolerance = 1.0e-8
			MaxIterations    = 500

		converger = DIISSCFConverger.WithOptions( energyTolerance   = EnergyTolerance,
                                                  densityTolerance  = DensityTolerance,
                                                  maximumIterations = MaxIterations )
	
	#------------------------------------- 
	def Change_QC_Region(self,_parameters):
		'''
		'''
		pass 
#=================================================
