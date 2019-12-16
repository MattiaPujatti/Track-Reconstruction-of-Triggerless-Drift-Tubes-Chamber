#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 15:26:24 2019

@author: mattia
"""




# grouping phisycal parameters

# cell's dimensions
XCELL = 42.   #mm
ZCELL = 13.   	#mm

# chamber's dimensions
LENGHT = 700.         #mm
WIDTH = 4*ZCELL    #mm


# Chambers 3D location
#LOCATION = {
#        
#        0: []
#            
#            }




# drift's parameters
TMAX = 390. # ns
Vd = XCELL/(2*TMAX) # mm/ns

# shift to global coordinates
X_SHIFT = [0,0,0,0]
Z_SHIFT = [0,62,782,844]

# error on x_points for local fit
XERR = 4.2

# expected muon tracks regions (local coordinates)
#REG = {0:[540,720],1:[540,720],2:[0,200],3:[0,200]}

# Chi-square limit
Chi_square_max = 999

# slope limits for local fit
Slope_limits = [1,999]







