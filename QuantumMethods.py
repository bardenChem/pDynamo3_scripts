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
from pBabel                    import *                                     
from pCore                     import *                                     
from pSimulation               import PruneByAtom
from pMolecule.MMModel         import *
from pMolecule.NBModel         import *                                     
from pMolecule.QCModel         import *
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
		self.selection      = [] 
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
		self = selfClass()
		self.system = _parameters["active_system"]
		self.methodClass = _parameters["method_class"]
		if "region" in _parameters: self.Hybrid = True

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
       	#--------------------
		self.Set_Converger()
		self.system.electronicState = ElectronicState.WithOptions( charge = _parameters["QCcharge"],multiplicity = _parameters["multiplicity"] )
		#----------------------------
		if self.methodClass == "SMO":
			self.qcModel = QCModelMNDO.WithOptions( hamiltonian = _parameters["Hamiltonian"], converger=self.converger )
		#...................................
		elif self.methodClass == "abinitio":
			self.qcModel = QCModelDFT.WithOptions(converger   = self.converger 			 ,
												  functional  = _parameters["functional"],
												  orbitalBasis= _parameters["basis"] 	 ,
												  fitBasis    = _parameters["fit_basi"]  )
		#...................................
		elif self.methodClass == "ORCA":
			options =  "\n% output\n"
        	options += "print [ p_mos ] 1\n"
        	options += "print [ p_overlap ] 5\n"
        	options += "end # output"
			self.qcModel = QCModelORCA.WithOptions( keywords = [ _parameters["functional"],_parameters["basis"],options  ], 
                                            deleteJobFiles  = False                      							  ,
                                            scratch         =_parameters["scratch"]                                   )
			NBmodel  = NBModelORCA.WithDefaults()
		#...................................
		elif self.methodClass == "DFTB":
			self.qcModel = QCModelDFTB.WithOptions( deleteJobFiles = False   ,
			                                		randomScratch  = True    ,
			                                 		scratch        = _parameters["scratch"],
			                                 		skfPath        = _parameters["skfPath"],
			                                 		useSCC         = True    )
			NBmodel = NBModelDFTB.WithDefaults()
		#...................................
		elif self.methodClass == "pySCF":
			pass
		#------------------------------------------------------------------------
		if self.Hybrid: self.system.DefineQCModel( self.qcModel, qcSelection=self.selection )
		else: self.system.DefineQCModel(self.qcModel)
		self.system.DefineNBModel( NBmodel )
        #------------------------------------------------------------------------
		energy = self.system.Energy()
		return(self)
        
	#-----------------------------------------
	def Export_QC_System(self,baseName = None):
		'''
		'''		
		if baseName == None: baseName = os.getcwd()
		self.qcSystem = PruneByAtom(self.system,self.selection)
		ExportSystem( os.path.join( baseName,"qcSystem.pdb"), self.qcSystem ) 
		ExportSystem( os.path.join( baseName,"qcSystem.pkl"), self.qcSystem )
	#-----------------------------------------
	def Set_Converger(self,_parameters=None):
		'''
		'''
		EnergyTolerance  = 3.0e-4
		DensityTolerance = 1.0e-8
		MaxIterations    = 500 

		if self.convergerLevel == "Very_Loose":
			EnergyTolerance  = 1.0e-4
			DensityTolerance = 1.0e-6
			MaxIterations    = 3500
		elif  self.convergerLevel == "Loose":
			EnergyTolerance  = 3.0e-4
			DensityTolerance = 1.0e-7
			MaxIterations    = 2500
		elif self.convergerLevel == "standard":
			EnergyTolerance  = 3.0e-4
			DensityTolerance = 1.0e-8
			MaxIterations    = 1000
		elif self.convergerLevel == "Tight":
			EnergyTolerance  = 1.0e-5
			DensityTolerance = 1.0e-8
			MaxIterations    = 800
		elif self.convergerLevel == "Very_Tight":
			EnergyTolerance  = 3.0e-6
			DensityTolerance = 1.0e-8
			MaxIterations    = 500

		self.converger = DIISSCFConverger.WithOptions( energyTolerance   = EnergyTolerance,
												  densityTolerance  = DensityTolerance,
												  maximumIterations = MaxIterations )
	
	#------------------------------------- 
	def Change_QC_Region(self,_parameters):
		'''
		'''
		pass 
	#-------------------------------------

#*********************************************
def QuantumMethods_UnitTest(_parameters):
	'''
	'''					
	qs =  QuantumMethods.From_Parameters(_parameters)
	qs.Export_QC_System()
	return(qs.system)
	
#=================================================
