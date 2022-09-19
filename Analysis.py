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
#--------------------------------------------------------------
VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')
if not VISMOL_HOME == None: sys.path.append(os.path.join(VISMOL_HOME,"easyhybrid/pDynamoMethods") ) 
else:                       sys.path.append(os.path.join("/home/igorchem/easyhybrid/pDynamoMethods") ) 
#-----------------------------------------------------------------------------------------------------
#Loading own libraries
#-------------------------------------------------------------
from EnergyAnalysis     	import EnergyAnalysis
from TrajectoryAnalysis 	import TrajectoryAnalysis
#-------------------------------------------------------------
from GeometrySearcher 	    import GeometrySearcher
from RelaxedScan 			import SCAN
from MolecularDynamics  	import MD
from UmbrellaSampling  	    import US
from PotentialOfMeanForce   import PMF
from ReactionCoordinate 	import ReactionCoordinate
from EnergyRefinement	 	import EnergyRefinement
#--------------------------------------------------------------
#loading pDynamo Libraries
from pBabel                    import *                                     
from pCore                     import *
#---------------------------------------                                     
from pMolecule                 import *                              
from pMolecule.MMModel         import *
from pMolecule.NBModel         import *                                     
from pMolecule.QCModel         import *
#---------------------------------------
from pScientific               import *                                     
from pScientific.Arrays        import *                                     
from pScientific.Geometry3     import *                                     
from pScientific.RandomNumbers import *                                     
from pScientific.Statistics    import *
from pScientific.Symmetry      import *
#--------------------------------------                                    
from pSimulation               import *
#=======================================================================
class Analysis:
	'''
	'''
	pass
