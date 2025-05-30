#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = commonFunctions.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################
import pBabel
from  pCore import Clone
import pMolecule                            
from pMolecule.MMModel         import *
from pMolecule.NBModel         import *                                     
from pMolecule.QCModel         import *
import numpy as np
import os, sys, glob
import math
#===============================================================================
orcaScratchBase="/home/igorchem/CCDIR/scratch"
skfPath        ="/home/igorchem/CCDIR/3ob-3-1"
#==============================================================================
#Atom dictionary with relevant information.
atomic_dic = {#Symbol     name         number    Cov(r)     VdW(r)     Mass
                "H"  : ["Hydrogen"     , 1   ,  0.330000 , 1.200000,  1.007940   ],
                "He" : ["Helium"       , 2   ,  0.700000 , 1.400000,  4.002602   ],
                "Li" : ["Lithium"      , 3   ,  1.230000 , 1.820000,  6.941000   ],
                "Be" : ["Beryllium"    , 4   ,  0.900000 , 1.700000,  9.012182   ],
                "B"  : ["Boron"        , 5   ,  0.820000 , 2.080000,  10.811000  ],
                "C"  : ["Carbon"       , 6   ,  0.770000 , 1.950000,  12.010700  ],
                "N"  : ["Nitrogen"     , 7   ,  0.700000 , 1.850000,  14.006700  ],
                "O"  : ["Oxygen"       , 8   ,  0.660000 , 1.700000,  15.999400  ],
                "F"  : ["Fluorine"     , 9   ,  0.611000 , 1.730000,  18.998404  ],
                "Ne" : ["Neon"         , 10  ,  0.700000 , 1.540000,  20.179701  ],
                "Na" : ["Sodium"       , 11  ,  3.06     , 2.270000,  22.989771  ],
                "Mg" : ["Magnesium"    , 12  ,  1.360000 , 1.730000,  24.305000  ],
                "Al" : ["Aluminium"    , 13  ,  1.180000 , 2.050000,  26.981539  ],
                "Si" : ["Silicon"      , 14  ,  0.937000 , 2.100000,  28.085501  ],
                "P"  : ["Phosphorus"   , 15  ,  0.890000 , 2.080000,  30.973761  ],
                "S"  : ["Sulphur"      , 16  ,  1.040000 , 2.000000,  32.064999  ],
                "Cl" : ["Chlorine"     , 17  ,  0.997000 , 1.970000,  35.452999  ],
                "Ar" : ["Argon"        , 18  ,  1.740000 , 1.880000,  39.948002  ],
                "K"  : ["Potassium"    , 19  ,  2.030000 , 2.750000,  39.098301  ],
                "Ca" : ["Calcium"      , 20  ,  1.740000 , 1.973000,  40.077999  ],
                "Sc" : ["Scandium"     , 21  ,  1.440000 , 1.700000,  44.955910  ],
                "Ti" : ["Titanium"     , 22  ,  1.320000 , 1.700000,  47.867001  ],
                "V"  : ["Vanadium"     , 23  ,  1.220000 , 1.700000,  50.941502  ],
                "Cr" : ["Chromium"     , 24  ,  1.180000 , 1.700000,  51.996101  ],
                "Mn" : ["Manganese"    , 25  ,  1.170000 , 1.700000,  54.938049  ],
                "Fe" : ["Iron"         , 26  ,  1.170000 , 1.700000,  55.845001  ],
                "Co" : ["Cobalt"       , 27  ,  1.160000 , 1.700000,  58.933201  ],
                "Ni" : ["Nickel"       , 28  ,  1.150000 , 1.630000,  58.693401  ],
                "Cu" : ["Copper"       , 29  ,  1.170000 , 1.400000,  63.546001  ],
                "Zn" : ["Zinc"         , 30  ,  1.250000 , 1.390000,  65.408997  ],
                "Ga" : ["Gallium"      , 31  ,  1.260000 , 1.870000,  69.723000  ],
                "Ge" : ["Germanium"    , 32  ,  1.188000 , 1.700000,  72.639999  ],
                "As" : ["Arsenic"      , 33  ,  1.200000 , 1.850000,  74.921600  ],
                "Se" : ["Selenium"     , 34  ,  1.170000 , 1.900000,  78.959999  ],
                "Br" : ["Bromine"      , 35  ,  1.167000 , 2.100000,  79.903999  ],
                "Kr" : ["Krypton"      , 36  ,  1.910000 , 2.020000,  83.797997  ],
                "Rb" : ["Rubidium"     , 37  ,  2.160000 , 1.700000,  85.467796  ],
                "Sr" : ["Strontium"    , 38  ,  1.910000 , 1.700000,  87.620003  ],
                "Y"  : ["Yttrium"      , 39  ,  1.620000 , 1.700000,  88.905853  ],
                "Zr" : ["Zirconium"    , 40  ,  1.450000 , 1.700000,  91.223999  ],
                "Nb" : ["Niobium"      , 41  ,  1.340000 , 1.700000,  92.906380  ],
                "Mo" : ["Molybdenum"   , 42  ,  1.300000 , 1.700000,  95.940002  ],
                "Tc" : ["Technetium"   , 43  ,  1.270000 , 1.700000,  98.000000  ],
                "Ru" : ["Ruthenium"    , 44  ,  1.250000 , 1.700000,  101.070000 ],
                "Rh" : ["Rhodium"      , 45  ,  1.250000 , 1.700000,  102.905502 ],
                "Pd" : ["Palladium"    , 46  ,  1.280000 , 1.630000,  106.419998 ],
                "Ag" : ["Silver"       , 47  ,  1.340000 , 1.720000,  107.868202 ],
                "Cd" : ["Cadmium"      , 48  ,  1.480000 , 1.580000,  112.411003 ],
                "In" : ["Indium"       , 49  ,  1.440000 , 1.930000,  114.818001 ],
                "Sn" : ["Tin"          , 50  ,  1.385000 , 2.170000,  118.709999 ],
                "Sb" : ["Antimony"     , 51  ,  1.400000 , 2.200000,  121.760002 ],
                "Te" : ["Tellurium"    , 52  ,  1.378000 , 2.060000,  127.599998 ],
                "I"  : ["Iodine"       , 53  ,  1.387000 , 2.150000,  126.904472 ],
                "Xe" : ["Xenon"        , 54  ,  1.980000 , 2.160000,  131.292999 ],
                "Cs" : ["Cesium"       , 55  ,  2.350000 , 1.700000,  132.905457 ],
                "Ba" : ["Barium"       , 56  ,  1.980000 , 1.700000,  137.326996 ],
                "La" : ["Lanthanum"    , 57  ,  1.690000 , 1.700000,  138.905502 ],
                "Ce" : ["Cerium"       , 58  ,  1.830000 , 1.700000,  140.115997 ],
                "Pr" : ["Praseodymium" , 59  ,  1.820000 , 1.700000,  140.907654 ],
                "Nd" : ["Neodymium"    , 60  ,  1.810000 , 1.700000,  144.240005 ],
                "Pm" : ["Promethium"   , 61  ,  1.800000 , 1.700000,  145.000000 ],
                "Sm" : ["Samarium"     , 62  ,  1.800000 , 1.700000,  150.360001 ],
                "Eu" : ["Europium"     , 63  ,  1.990000 , 1.700000,  151.964005 ],
                "Gd" : ["Gadolinium"   , 64  ,  1.790000 , 1.700000,  157.250000 ],
                "Tb" : ["Terbium"      , 65  ,  1.760000 , 1.700000,  158.925339 ],
                "Dy" : ["Dysprosium"   , 66  ,  1.750000 , 1.700000,  162.500000 ],
                "Ho" : ["Holmium"      , 67  ,  1.740000 , 1.700000,  164.930313 ],
                "Er" : ["Erbium"       , 68  ,  1.730000 , 1.700000,  167.259003 ],
                "Tm" : ["Thulium"      , 69  ,  1.720000 , 1.700000,  168.934204 ],
                "Yb" : ["Ytterbium"    , 70  ,  1.940000 , 1.700000,  173.039993 ],
                "Lu" : ["Lutetium"     , 71  ,  1.720000 , 1.700000,  174.966995 ],
                "Hf" : ["Hafnium"      , 72  ,  1.440000 , 1.700000,  178.490005 ],
                "Ta" : ["Tantalum"     , 73  ,  1.340000 , 1.700000,  180.947906 ],
                "W"  : ["Tungsten"     , 74  ,  1.300000 , 1.700000,  183.839996 ],
                "Re" : ["Rhenium"      , 75  ,  1.280000 , 1.700000,  186.207001 ],
                "Os" : ["Osmium"       , 76  ,  1.260000 , 1.700000,  190.229996 ],
                "Ir" : ["Iridium"      , 77  ,  1.270000 , 1.700000,  192.216995 ],
                "Pt" : ["Platinum"     , 78  ,  1.300000 , 1.720000,  195.078003 ],
                "Au" : ["Gold"         , 79  ,  1.340000 , 1.660000,  196.966553 ],
                "Hg" : ["Mercury"      , 80  ,  1.490000 , 1.550000,  200.589996 ],
                "Tl" : ["Thallium"     , 81  ,  1.480000 , 1.960000,  204.383301 ],
                "Pb" : ["Lead"         , 82  ,  1.480000 , 2.020000,  207.199997 ],
                "Bi" : ["Bismuth"      , 83  ,  1.450000 , 1.700000,  208.980377 ],
                "Po" : ["Polonium"     , 84  ,  1.460000 , 1.700000,  209.000000 ],
                "At" : ["Astatine"     , 85  ,  1.450000 , 1.700000,  210.000000 ],
                "Rn" : ["Radon"        , 86  ,  2.400000 , 1.700000,  222.000000 ],
                "Fr" : ["Francium"     , 87  ,  2.000000 , 1.700000,  223.000000 ],
                "Ra" : ["Radium"       , 88  ,  1.900000 , 1.700000,  226.000000 ],
                "Ac" : ["Actinium"     , 89  ,  1.880000 , 1.700000,  227.000000 ],
                "Th" : ["Thorium"      , 90  ,  1.790000 , 1.700000,  232.038101 ],
                "Pa" : ["Protactinium" , 91  ,  1.610000 , 1.700000,  231.035873 ],
                "U"  : ["Uranium"      , 92  ,  1.580000 , 1.860000,  238.028915 ],
                "Np" : ["Neptunium"    , 93  ,  1.550000 , 1.700000,  237.000000 ],
                "Pu" : ["Plutionium"   , 94  ,  1.530000 , 1.700000,  244.000000 ],
                "Am" : ["Americium"    , 95  ,  1.070000 , 1.700000,  243.000000 ],
                "Cm" : ["Curium"       , 96  ,  0.000000 , 1.700000,  247.000000 ],
                "Bk" : ["Berkelium"    , 97  ,  0.000000 , 1.700000,  247.000000 ],
                "Cf" : ["Californium"  , 98  ,  0.000000 , 1.700000,  251.000000 ],
                "Es" : ["Einsteinium"  , 99  ,  0.000000 , 1.700000,  252.000000 ],
                "Fm" : ["Fermium"      , 100 ,  0.000000 , 1.700000,  257.000000 ],
                "Md" : ["Mendelevium"  , 101 ,  0.000000 , 1.700000,  258.000000 ],
                "No" : ["Nobelium"     , 102 ,  0.000000 , 1.700000,  259.000000 ],
                "Lr" : ["Lawrencium"   , 103 ,  0.000000 , 1.700000,  262.000000 ],
                "Rf" : ["Rutherfordiu" , 104 ,  0.000000 , 1.700000,  261.000000 ],
                "Db" : ["Dubnium"      , 105 ,  0.000000 , 1.700000,  262.000000 ],
                "Sg" : ["Seaborgium"   , 106 ,  0.000000 , 1.700000,  263.000000 ],
                "Bh" : ["Bohrium"      , 107 ,  0.000000 , 1.700000,  264.000000 ],
                "Hs" : ["Hassium"      , 108 ,  0.000000 , 1.700000,  265.000000 ],
                "Mt" : ["Meitnerium"   , 109 ,  0.000000 , 1.700000,  268.000000 ],
                "Xx" : ["Dummy"        , 0   ,  0.000000 , 0.000000,  0.000000   ],
                "X"  : ["Dummy"        , 0   ,  0.000000 , 0.000000,  0.000000   ]
              }
#==============================================================================
#possible SMO names
SMOnames = ["am1","am1dphot","mndostong","pddgmndo","pddgpm3","pm3","pm6","rm1"]
#==============================================================================
def VerifyMNDOKey( key ):
    for smo in SMOnames:
        if key == smo:
            return(True)
    print("Invalid name for MNDO Hamiltonin Defined on pDynamo3")
    return(False)
#==============================================================================
def GetTotalCharge(_system):
    '''    
    Calculate the total charge from the atoms of the
    passed System instance object.
    Parameter #1: System instanced object
    Return: Double holding the total charge
    '''
    totalCharge = 0
    Charge = _system.energyModel.mmAtoms.AtomicCharges()
    for i in range( len(Charge) ):
        totalCharge += Charge[i]
    return (totalCharge)
#==============================================================================
def ReescaleCharges(_system, tc):
    '''Function to reescale atomic charges from System instance object
    and make their sum a real integer. 
    Parameter #1: System instanced object
    Parameter #2: Total charge to be achieved 
    Return: System instanced object with the scaled charges.
    '''    
    scaled_system = Clone(_system)    
    Charges = scaled_system.energyModel.mmAtoms.AtomicCharges()
    pTC = GetTotalCharge(scaled_system)    
    print ("Old Charges Sum:",pTC)    
    nAtoms = len(scaled_system.atoms.items)
    #---------------------------------------------------------
    new_charges = []
    new_tc = 0
    frac = (tc - pTC)/nAtoms
    #---------------------------------------------------------
    for i in range( len(Charges) ):
        new_charges.append(Charges[i] + frac)    
    for i in range( len(new_charges) ):
        Charges[i] = new_charges[i]    
    #---------------------------------------------------------
    scaled_system.energyModel.mmAtoms.SetAtomicCharges(Charges)    
    new_tc = GetTotalCharge(scaled_system)
    print ("New Charges Sum:",new_tc)        
    return(scaled_system)        
#==============================================================================
def copySystem(system):
    '''
    Make a System deepy copy handling with the Non-bonded methods definitions 
    '''
    nbmodel_hold   = system.nbModel
    system.nbModel = None
    newSystem      = Clone(system)
    newSystem.DefineNBModel ( nbmodel_hold )
    system.nbModel = nbmodel_hold      
    return newSystem
#==============================================================================
def GetAtomicMass(atomN):
    '''
    Return the atomic mass float value from the atomic information dictionary
    '''
    ls = list( atomic_dic.values() )
    atomMass = ls[atomN-1][4]
    return (atomMass)
#==============================================================================
def GetAtomicSymbol(atomN):
    '''
    Return the atomic element symbol string from the atomic information dictionary
    '''
    ls = list( atomic_dic )
    _symbol = ls[atomN-1]
    return(_symbol)
#=========================================================================================
def GetFrameIndex(fname):
    '''
    Get the indices of the frame from string name
    Pass file name without extension
    '''
    idxs = []
    wkstr  = os.path.basename(fname)
    ssplit = wkstr.split("_")
    if len(ssplit) == 1:
        if wkstr[:5] == "frame":
            if "." in wkstr:
                wkstr = wkstr.split(".")
                idxs.append( int( wkstr[0][5:] ) )
            else: idxs.append( int( wkstr[5:] ) )
    elif len(ssplit) == 2: 
        if ssplit[0][:5] == "frame":
            idxs.append( int( ssplit[0][5:] ) )
            if "." in ssplit[1]:
                wkstr = ssplit[1].split(".")
                idxs.append( int(wkstr[0]) )   
            else: idxs.append( int(ssplit[1]) )
    return(idxs)
#=========================================================================================
def write_base_input():
    '''
    '''
    pass