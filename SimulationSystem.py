#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = CoreInterface.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

import os, glob, sys
#--------------------------------------------------------------
#Loading own libraries
from commonFunctions import *
from Simulation import Simulation
from Analysis import *
from QuantumMethods import *
#--------------------------------------------------------------
#loading pDynamo Libraries
from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                              
from pMolecule.MMModel         import *
from pMolecule.NBModel         import *                                     
from pMolecule.QCModel         import *
from pSimulation               import PruneByAtom

from pBabel            import ExportSystem                    , \
                              GromacsDefinitionsFileReader    , \
                              GromacsParameterFileReader      , \
                              ImportCoordinates3              , \
                              ImportSystem
#==========================================================================

#**************************************************************************
class SimulationSystem:
    '''
    Class to Wrapper information around Sysyem class from pDynamo3 Libraries
    '''  
    #.-------------------------------------------------------------------------
    def __init__(self,_label="No specified"):
        '''
        Class constructor
        Parameters:
        '''        
        self.baseName            = None
        self.label               = _label
        self.system              = None # Instance of the System pDnamo class
        self.NBmodel             = None 
        self.QCmodel             = None
        self.MMmodel             = None
        self.Hybrid              = None 
        self.quantumRegion       = []
        self.protein             = False 
        self.ReactionCoordinates = [] 
        self.refEnergy           = 0.0 

    #===================================================================================
    @classmethod
    def From_PKL(selfClass,_pklPath,_FolderName=None):
        '''
        Initialize project object from a pdynamo PKL file.
        '''
        self                = selfClass()
        self.system         = ImportSystem(_pklPath)
        _name               = os.path.basename(_pklPath)
        self.baseName       = _name[:-4]
        # test if is quantum, if hass mmModel and NbModel
        return(self)
    #=================================================================================== 
    @classmethod
    def From_Force_Field(selfClass,_parameters):
        '''
        Initialize project from force field topology and coordinate files.
        '''
        self = selfClass()
        self.NBmodel = NBModelCutOff.WithDefaults()      
        self.system               = ImportSystem(_topologyFile)
        self.system.DefineNBModel = self.NBmodel
        self.system.coordinates3  = ImportCoordinates3(_coordinateFile)               
        _name                     = os.path.basename(_topologyFile)
        self.baseName             = _name[:-4]        
        return(self)     
    #===================================================================================
    @classmethod
    def From_Gromacs(selfClass,_parameters):
        '''
        '''
        self = selfClass()
        fileName     = os.path.join                                ( dataPath, label + "_" + ff )
        parameters   = GromacsParameterFileReader.PathToParameters ( fileName + ".top" )
        self.system  = GromacsDefinitionsFileReader.PathToSystem   ( fileName + ".top", parameters = parameters )
        self.system.coordinates3 = ImportCoordinates3              ( fileName + ".gro" )
        return(self)        
    #===================================================================================
    @classmethod
    def From_Coordinates(selfClass,_coordinateFile):
        '''
        Initialize project from coordinates 
        '''
        self = selfClass()        
        self.system   = ImportSystem(_coordinateFile)
        _name         = os.path.basename(_coordinateFile)
        self.baseName = _name[:-4]
        return(self) 
    #===================================================================================
    @classmethod
    def Protein_From_Coordinates(selfClass,_coordinateFile,_modelNumber=1):
        '''
        Initialize project from coordinate file with OPLS general force field
        '''
        self = selfClass()
        self.system      = ImportSystem(_coordinateFile, modelNumber = _modelNumber, useComponentLibrary = True)
        self.MMmodel = MMModelOPLS.WithParameterSet("protein")
        self.NBmodel = NBModelCutOff.WithDefaults()
        self.system.DefineNBModel(self.NBmodel)                     
        _name        = os.path.basenem(_coordinateFile)
        self.baseName= _name[:-4]
        self.protein = True
        return(self) 
    #====================================================================================
    @property
    def Check(self):
        '''
        Calculates single point energy.
        '''
        self.system.Summary()
        self.refEnergy = self.system.Energy(doGradients = True)
        return self.refEnergy    
    #====================================================================================
    def Spherical_Pruning(self,_centerAtom,_radius):
        '''
        Perform a spherical pruning from a certain atom center coordinates.
        Parameters:
            _centerAtom:
            _radius    :
        '''
        #---------------------------------------------------
        oldSystem = copySystem(self.system)
        #---------------------------------------------------
        atomref      = AtomSelection.FromAtomPattern( oldSystem, _centerAtom )
        core         = AtomSelection.Within(oldSystem,atomref,_radius)
        core2        = AtomSelection.ByComponent(oldSystem,core)
        #---------------------------------------------------
        newLabel     = self.label + "_pruned"
        self.system  = PruneByAtom( oldSystem,Selection(core2) )
        self.system.DefineNBModel( self.NBmodel )        
    #======================================================================================
    def Setting_Free_Atoms(self,_centerAtom,_radius,_DEBUG=False):
        '''
        Set the list of atoms to keep with the positions fixed through the next simulations
        Parameters:
            _centerAtom:
            _radius    :
        '''
        newSystem = copySystem(self.system)
        self.systemCoutCurr += 1
        #-----------------------------------------------------
        atomref = AtomSelection.FromAtomPattern( newSystem, _centerAtom )
        core    = AtomSelection.Within(newSystem,atomref,_radius)
        mobile  = AtomSelection.ByComponent(newSystem,core)        
        #-----------------------------------------------------
        newLabel= self.system.label + "_fixed"
        if _DEBUG:
            MobileSys = PruneByAtom( newSystem, Selection(mobile) )
            ExportSystem("MobileSystemCheck.pdb",MobileSys)
        #------------------------------------------------------
        self.systems = newSystem
        self.systems.freeAtoms = mobile       
        self.system.DefineNBModel( self.NBmodel )
        self.system = self.systems[newLabel]
        if self.DEBUG: self.Energy
    
    #=========================================================================
    def Set_QC_Method(self,_parameters,_DEBUG=false):
        '''
        '''        
        _parameters["active_system"] = self.system 
        qs =  QuantumMethods.From_Parameters(_parameters)
        if _DEBUG: qs.Export_QC_System()
        newLabel = self.system.label + "QC_system_"
        if "Hamiltonian" in _parameters: newLabel += _parameters["Hamiltonian"] 
        if "functional" in _parameters: newLabel  += _parameters["functional"] 
        self.system = qs.system

    #=========================================================================
    def Set_Reaction_Coordinates(self,_parameters):
        '''
        '''
        pass

#================================================================================

