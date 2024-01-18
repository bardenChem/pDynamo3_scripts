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
#==========================================================================

#**************************************************************************
class SimulationProject:
    '''
    Class to setup pDynamo simulations from a remote framework, i.e. without using VisMol GUI 
    '''  
    #.-------------------------------------------------------------------------
    def __init__(self,_projectFolder=None,DEBUG=False):
        '''
        Class constructor
        Parameters:
            _projectFolder: name of the project and the main folder where all results will be produced. Please pass the full path
            DEBUG       : if this paramters it set True, some extra steps in the system setting will be performe to check the things up
        '''        
        self.baseName       = None
        self.folderName     = _projectFolder
        self.systems        = {} #System Dict loaded in this project 
        self.simulation     = [] #List given the simulations executed in this project
        self.systemKeys     = []
        self.systemCoutCurr = 0 
        self.system         = None
        self.DEBUG          = DEBUG
        #Some status variable for current Sytem instance ( For now without any great utility )
        self.NBmodel        = None 
        self.QCmodel        = {} #QC models Dict loaded in this project
        self.MMmodel        = {} #MM models Dict loaded in this project
        self.refEnergy      = 0.0 

        if not self.folderName:
            self.folderName = os.getcwd()
        if not os.path.exists(self.folderName): os.makedirs(self.folderName)

    #===================================================================================
    @classmethod
    def From_PKL(selfClass,_pklPath,_FolderName=None):
        '''
        Initialize project object from a pdynamo PKL file.
        '''
        self                = selfClass(_projectFolder=_FolderName)
        _system             = ImportSystem(_pklPath)
        _name               = os.path.basename(_pklPath)
        self.baseName       = _name[:-4]
        self.systems[_name] = _system
        self.system         = _system
        self.systemCoutCurr += 1 
        self.systemKeys.append(_name)
        #-----------------------------
        return(self)
    #=================================================================================== 
    @classmethod
    def From_Force_Field(selfClass,_topologyFile,_coordinateFile,_FolderName=None):
        '''
        Initialize project from force field topology and coordinate files.
        '''
        self = selfClass(_projectFolder=_FolderName)
        #...........................................
        self.NBmodel = NBModelCutOff.WithDefaults()      
        #...........................................
        _system              = ImportSystem(_topologyFile)
        _system.DefineNBModel(self.NBmodel) 
        _system.coordinates3 = ImportCoordinates3(_coordinateFile)               
        _name                = os.path.basename(_topologyFile)
        self.baseName        = _name[:-4]
        self.systems[_name]  = _system
        self.MMmodel[_name]  = _system.mmModel
        self.systemKeys.append(_name)
        self.system          = _system
        self.systemCoutCurr +=1        
        #-----------------------------
        return(self)        
    #===================================================================================
    @classmethod
    def From_Coordinates(selfClass,_coordinateFile,_FolderName=None):
        '''
        Initialize project from coordinates 
        '''
        self = selfClass(_projectFolder=_FolderName)
        #...........................................
        _system      = ImportSystem(_coordinateFile)
        _system.coordinates3 = ImportCoordinates3(_coordinateFile)               
        _name                = os.path.basename(_coordinateFile)
        self.baseName        = _name[:-4]
        self.systems[_name]  = _system
        #self.MMmodels[_name] = _system.mmModel
        self.system          = _system
        self.systemKeys.append(_name)
        #-----------------------------
        return(self) 
    #===================================================================================
    @classmethod
    def Protein_From_Coordinates(selfClass,_coordinateFile,_FolderName=None):
        '''
        Initialize project from coordinate file with OPLS general force field
        '''
        self = selfClass(_projectFolder=_FolderName)
        #...........................................
        _system      = ImportSystem(_coordinateFile)
        self.NBmodel = NBModelCutOff.WithDefaults()
        self.MMmodel = MMModelOPLS.WithParameterSet("protein") 
        _system.DefineNBModel(self.NBmodel) 
        #...........................................
        _system              = ImportSystem(_coordinateFile)
        _system.coordinates3 = ImportCoordinates3(_coordinateFile)               
        _name                = os.path.basenem(_coordinateFile)
        self.baseName        = _name[:-4]
        self.systems[_name]  = _system
        self.MMmodels[_name] = _system.mmModel
        self.system          = _system
        self.systemKeys.append(_name)  
        #-----------------------------
        return(self) 
    #====================================================================================
    @property
    def Energy(self):
        '''
        Calculates single point energy.
        '''
        self.refEnergy = self.system.Energy()
        return self.refEnergy    
    #====================================================================================
    def Spherical_Pruning(self,_centerAtom, _radius):
        '''
        Perform a spherical pruning from a certain atom center coordinates.
        Parameters:
            _centerAtom:
            _radius    :
        '''
        #---------------------------------------------------
        oldSystem = copySystem(self.system)
        self.systemCoutCurr += 1
        #---------------------------------------------------
        atomref      = AtomSelection.FromAtomPattern( oldSystem, _centerAtom )
        core         = AtomSelection.Within(oldSystem,atomref,_radius)
        core2        = AtomSelection.ByComponent(oldSystem,core)
        #---------------------------------------------------
        newLabel     = self.systemKeys[-1] + "_pruned"
        self.systems[newLabel] = PruneByAtom( oldSystem,Selection(core2) )
        self.systems[newLabel].DefineNBModel( self.NBmodel )
        self.system = self.systems[newLabel]
        #---------------------------------------------------
        if self.DEBUG: self.Energy
    #======================================================================================
    def Setting_Free_Atoms(self,_centerAtom,_radius):
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
        newLabel= self.systemKeys[-1] + "_fixed"
        if self.DEBUG:
            MobileSys = PruneByAtom( newSystem, Selection(mobile) )
            ExportSystem("MobileSystemCheck.pdb",MobileSys)
        #------------------------------------------------------
        self.systems[newLabel] = newSystem
        self.systems[newLabel].freeAtoms = mobile       
        self.systems[newLabel].DefineNBModel( self.NBmodel )
        self.system = self.systems[newLabel]
        if self.DEBUG: self.Energy
    
    #=========================================================================
    def Set_QC_Method(self,_parameters):
        '''
        '''        
        _parameters["active_system"] = self.system 
        qs =  QuantumMethods.From_Parameters(_parameters)
        if self.DEBUG: qs.Export_QC_System()
        newLabel = "QC_system_"
        print(_parameters)
        if "Hamiltonian" in _parameters: newLabel += _parameters["Hamiltonian"] 
        if "functional" in _parameters: newLabel  += _parameters["functional"] 
        self.systems[newLabel]                     = qs.system
        self.system                                = self.systems[newLabel]
    #=========================================================================
    def Set_QCMM_SMO(self,_parameters):
        '''
        Define hybrid system initial configutation. 
        '''
        atomlist = []
        for sel in _parameters["region"]:
            if type(sel) == int:
                atomlist.append(sel)
            elif type(_parameters["region"]) == list:
                for i in range( len(sel) ):
                    atomlist.append( sel[i] ) 
        #--------------------------------------------
        if self.DEBUG:
            print(atomlist)           
        #---------------------------------------------
        #define QC atoms selection
        converger = DIISSCFConverger.WithOptions( energyTolerance   = 3.0e-4,
                                                  densityTolerance  = 1.0e-8,
                                                  maximumIterations = 2500  )
        _QCRegion = Selection.FromIterable(atomlist)    
        #---------------------------------------------
        #Appending sytem
        self.system = self.systems[self.systemKeys[-1]]
        self.NBmodel = self.system.nbModel
        self.system.nbModel = None
        oldSystem = Clone(self.system)       
        self.system.label = self.baseName + "#{} {} Hamiltonian and QC region Set".format(self.systemCoutCurr,_parameters["Hamiltonian"])
        self.systemCoutCurr += 1
        #Setting QC model 
        self.system.electronicState = ElectronicState.WithOptions ( charge = _parameters["QCcharge"] , multiplicity = _parameters["multiplicity"] )
        _QCmodel = QCModelMNDO.WithOptions( hamiltonian = _parameters["Hamiltonian"], converger=converger )        
        #------------------------------------------------------------------------
        self.system.DefineQCModel( _QCmodel, qcSelection=_QCRegion )
        self.system.DefineNBModel( self.NBmodel )
        self.systems["QCMM_system"] = self.system
        #------------------------------------------------------------------------
        if self.DEBUG:
            qcSystem = PruneByAtom(self.system,_QCRegion)
            ExportSystem(self.baseName+"/qcSystem.pdb",qcSystem)
            ExportSystem(self.baseName+"/qcSystemEntire.pdb",self.system)
        energy = self.system.Energy()

    #=========================================================================
    def Run_Simulation(self,_parameters):
        '''
        Execute a preset simulation for the current system. 
        Parameters:
           _parameters: Dict with all paramters for setting the simulations.
        '''
        #----------------------------------------------------------------------
        #---------------------------------------------------------------------
        _parameters["active_system"] = self.system
        if not "folder" in _parameters:_parameters["folder"] = self.folderName        
        process = Simulation(_parameters)
        process.Execute()
    #=========================================================================
    def Run_Analysis(self,_parameters):
        '''
        Execute analysis for a given done simulation
        Parameters:
            _parameters: Dict with all possible parameters for analysis of the given simulation 
        '''
        _parameters["active_system"] = self.system
        if not "folder" in _parameters:_parameters["folder"] = self.folderName        
        process = Analysis(_parameters)
        process.Execute()
        pass 
    #========================================================================================
    def PrintSystems(self):
        '''
        Method to print the summary of the loaded systems 
        '''
        print("There are {} loaded systems".format( self.systemCoutCurr) )
        ctn = input("Type any key to print the Summary of the Systems, or 'N' to cancel this")
        if not ctn == "N":
            if len(self.SystemStates) > 0:
                for system in self.SystemStates:
                    system.Summary()
                    print("***************************************************")
                print("Now, printing the current system Summary:")
                self.System.Summary()

            elif self.systemCoutCurr == 1: print("There is only the current System loaded!\n Printing its information below!")
            else:                          print( "There are no loaded systems!")
    #==========================================================================================
    def SaveProject(self):
        '''
        The complete version of this function intends to save in pkl and another coordinate format
        the systems and trajectories worked in this simulations
        Though, in the current state only will save the current system to a pkl file
        '''
        savePathPkl = os.path.join(self.folderName,self.baseName+".pkl")
        savePathPdb = os.path.join(self.folderName,self.baseName+".pdb")

        for key in self.systems:
            savePathPkl = os.path.join(savePathPkl, self.baseName+"_"+key+".pkl")
            Pickle( savePathPkl, self.systems[key] )        
    #.-------------------------------------------------------------------------
    def SaveSystem( self, _cname=None):
        '''
        '''

        if _cname:
            savePathPkl = os.path.join(self.folderName,_cname+".pkl")
            savePathPdb = os.path.join(self.folderName,_cname+".pdb")
        else:
            savePathPkl = os.path.join(self.folderName,self.baseName+".pkl")
            savePathPdb = os.path.join(self.folderName,self.baseName+".pdb")

            i = 0;
            while os.path.exists(savePathPdb):
                i += 1
                savePathPdb = os.path.join(self.folderName,self.baseName+".pdb")
                savePathPdb = savePathPdb[:-4] + "_#{}.pdb".format(i)
            while os.path.exists(savePathPkl):
                i += 1
                savePathPkl = os.path.join(self.folderName,self.baseName+".pkl")
                savePathPkl = savePathPkl[:-4] + "_#{}.pkl".format(i)
        #----------------------------------------------------------------
        Pickle( savePathPkl,self.system )
        ExportSystem( savePathPdb,self.system )
        
    #.-------------------------------------------------------------------------
    def FinishRun(self):
        '''
        Finalize the run.
        '''
        
#==============================================================================
if __name__ == "__main__":
    #----------------------------------------------------
    '''
    Fast test
    scratch_path = os.path.join("TestsScratch")
    ex_path      = "examples"
    timTop       = os.path.join(ex_path,"TIM","7tim.top")
    timCrd       = os.path.join(ex_path,"TIM","7tim.crd")
    balapkl      = os.path.join(ex_path,"bala","bAla.pkl")
    meth         = os.path.join(ex_path,"pdb","methane.pdb")
    #----------------------------------------------------
    test1 = SimulationProject() #default initialization
    test2 = SimulationProject.From_PKL(balapkl,_FolderName="bala_test")
    print(test2.Energy)
    '''
