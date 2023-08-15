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
#-----------------------------------------------------------------------------------------------------
#Loading own libraries
#-------------------------------------------------------------
from EnergyAnalysis     	import EnergyAnalysis
from TrajectoryAnalysis 	import TrajectoryAnalysis
#-------------------------------------------------------------
#=======================================================================
class Analysis:
	'''
	'''
	def __init__(self,_parameters):
		'''
		Default Constructor
		'''
		self.parameters = _parameters


	def Execute(self):
		'''
		'''
		Type = self.parameters["analysis_type"]
		if   Type == "Trajectory_Analysis": self.TrajectoryPlots() 		
		elif Type == "Energy_Plots":		self.EnergyPlots()
		elif Type == "PMF":                 self.PMF()

	#=========================================================================	
	def TrajectoryPlots(self) :
		'''
		Mandatory keys in self.parameters:
		Optional keys in self.parameters:
		'''
		RCs  = None
		if "show" in self.parameters: show = self.parameters["show"]
		t_time = self.parameters["nsteps"]*0.001
		DA = TrajectoryAnalysis(MDrun.trajectoryNameCurr,self.molecule,t_time)
		DA.CalculateRG_RMSD()
		DA.PlotRG_RMS(show)						
		if "calculate_distances" in self.parameters:
			if self.parameters["calculate_distances"] == True:
				rc1 = ReactionCoordinate(self.parameters["ATOMS_RC1"],False,0)
				rc1.GetRCLabel(self.molecule)
				RCs = [rc1]
				rc2 = None
				if "ATOMS_RC2" in self.parameters:
					rc2 = ReactionCoordinate(self.parameters["ATOMS_RC2"],False,0)						
					rc2.GetRCLabel(self.molecule)
					RCs.append(rc2)
				DA.DistancePlots(RCs,show)
	#=========================================================================
	def EnergyPlots(self):
		'''
		Produce Energy plots from previus simulations log files
		Mandatory keys in self.parameters:
		Optional keys in self.parameters:
		'''		
		multiPlot = False
		ndim      = 1 
		crd1_label= "Reaction Coordinate #1"
		crd2_label= "Reaction Coordinate #2"
		cnt_lines = 0 
		ysize     = 0
		if "ysize" in self.parameters: ysize = self.parameters["ysize"]
		xlim      = [ 0, self.parameters["xsize"] ]
		ylim 	  = [ 0, ysize ]
		show 	  = False
		#--------------------------------------------------------
		if "contour_lines" 	in self.parameters: cnt_lines  = self.parameters["contour_lines"]
		if "crd1_label" 	in self.parameters: crd1_label = self.parameters["crd1_label"]
		if "crd2_label" 	in self.parameters:	crd2_label = self.parameters["crd2_label"]
		if "xlim_list" 		in self.parameters: xlim  	   = self.parameters["xlim_list"]
		if "ylim_list" 		in self.parameters: ylim       = self.parameters["ylim_list"]
		if "show" 			in self.parameters: show       = self.parameters["show"]
		if "multiple_plot" 	in self.parameters: multiPlot  = self.parameters["multiple_plot"]		
		if ysize > 0: ndim = 2
		if "log_names" in self.parameters: multiPlot = True 
		#--------------------------------------------------------
		EA = EnergyAnalysis(self.parameters["xsize"],ysize,_type=self.parameters["type"] )
		if multiPlot:
			for log in self.parameters["log_names"]:
				EA.ReadLog( log )
				EA.MultPlot1D()
		else:	EA.ReadLog( self.parameters["log_name"] )
		#--------------------------------------------------------
		if 	 ndim == 1: EA.Plot1D(crd1_label,XLIM=xlim,SHOW=show)
		elif ndim == 2:	EA.Plot2D(cnt_lines,crd1_label,crd2_label,xlim,ylim,show)

	#=========================================================================
	def PMFAnalysis(self):
		'''
		Calculate potential of mean force and Free energy from restricted molecular dynamics
		Mandatory keys: 
			"source_folder"	:
			"xbins"			:
			"ybins"			:
			"temperature"	:
		Optinal keys        :
		plot keys           :				
		'''
		ynbins = 0 
		if "ynbins" in self.parameters: ynbins = self.parameters["ynbins"]
		potmean = PMF( self.molecule, self.parameters["source_folder"], self.baseFolder )
		potmean.CalculateWHAM(self.parameters["xnbins"],ynbins,self.parameters["temperature"])
		#================================================================
		#Set default plot parameters
		cnt_lines  = 12
		crd1_label = ""
		crd2_label = ""
		nRC2       = ynbins
		show       = False
		xwin       = 0
		ywin       = 0 
		#-----------------------------------------------------------------
		nDims = 1
		if ynbins > 0: nDims = 2
		xlims = [ 0,  self.parameters['xnbins'] ]
		ylims = [ 0,  ynbins ]
		OneDimPlot = False
		#-------------------------------------------------------------
		#check parameters for plot
		if "contour_lines" 	in self.parameters: cnt_lines  = self.parameters["contour_lines"]		
		if "xlim_list" 		in self.parameters: xlims 	   = self.parameters["xlim_list"]
		if "ylim_list" 		in self.parameters:	ylims 	   = self.parameters["ylim_list"]
		if "show" 			in self.parameters:	show 	   = self.parameters["show"]
		if "crd1_label" 	in self.parameters:	crd1_label = self.parameters["crd1_label"]
		if "crd2_label" 	in self.parameters:	crd2_label = self.parameters["crd2_label"]
		if "xwindows" 		in self.parameters: xwin 	   = self.parameters["xwindows"]
		if "ywindows" 		in self.parameters:	ywin 	   = self.parameters["ywindows"]
		if "oneDimPlot"     in self.parameters: OneDimPlot = self.parameters["oneDimPlot"]
		#------------------------------------------------------------
		if   nDims == 2: TYPE = "WHAM2D"
		elif nDims == 1: TYPE = "WHAM1D"	
		#------------------------------------------------------------
		# Plot PMF graphs
		EA = EnergyAnalysis(self.parameters['xnbins'],nRC2,_type=TYPE)
		EA.ReadLog( os.path.join(potmean.baseName,"PotentialOfMeanForce.dat") ) 
		#-------------------------------------------------------------
		if   nDims == 2: EA.Plot2D(cnt_lines,crd1_label,crd2_label,xlims,ylims,show)
		elif nDims == 1: EA.Plot1D(crd1_label,SHOW=show)
		#-------------------------------------------
		#Plot Free energy of the calculated windows
		if 	 OneDimPlot == True: TYPE = "FE1D"
		elif nDims 		== 2: 	 TYPE = "FE2D"
		elif nDims 		== 1: 	 TYPE = "FE1D"

		xlims = [ np.min(EA.RC1), np.max(EA.RC1) ]

		if nDims  == 2:  ylims = [ np.min(EA.RC2), np.max(EA.RC2) ]	
		#------------------------------------------
		EAfe = EnergyAnalysis(xwin,ywin,_type=TYPE)
		EAfe.ReadLog( os.path.join(potmean.baseName,"FreeEnergy.log") ) 
		#-------------------------------------------------------------
		if nDims == 2: 
			if OneDimPlot: EAfe.Plot1D_FreeEnergy(crd1_label,crd2_label,show)
			else 		 : EAfe.Plot2D(cnt_lines,crd1_label,crd2_label,xlims,ylims,show)
		elif nDims == 1: EAfe.Plot1D(crd1_label,XLIM=xlims,SHOW=show)
