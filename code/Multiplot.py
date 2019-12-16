#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 21:08:37 2019

@author: mattia
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random

from tracks_building import Tracks_Builder
from dataframe_builder import DataFrameBuilder
from eventplot import Event_Plot
from Raw_data_plot import Raw_data_plots
import constants as c


class Multi_event_plots:
    
    def __init__(self,datafile,events):
        
        self.file = datafile
        self.ev = events
        
        
    def Check_single_event(self,number):
    
        if number == 'random':
            num_lines = sum(1 for line in self.file)
            number = random.randint(1,num_lines)
    
        Raw_data_plots(self.ev).Hit_Matrix_single_plot(number)
        df = DataFrameBuilder(self.file).DataFrame_getn(number,stamp=True)
        
        TB = Tracks_Builder(df)
        Event_Plot(df,TB).Hits_Plot(detectors_zoom=True,bkgplot=True)
    
        return df
        
    
    def Random_events_plot(self,how_many):
        
        num_lines = sum(1 for line in self.file)
        selected_events = pd.DataFrame()
        rand_lst = []
        
        good_cnt = 0
        
        fig = plt.figure()
        ax = fig.add_axes([1,1,1,1])

        if how_many == 'all':
            rand_lst = list(range(1,num_lines))
            
            for i in range (1,num_lines):
                df = DataFrameBuilder(self.file).DataFrame_getn(i)
                TB = Tracks_Builder(df)
                if not TB.Checkgoodevent():
                    continue
                Figure_plot(ax,df)
        else:
            while good_cnt < how_many:
                index = random.randint(1,num_lines)
                rand_lst.append(index)
                df = DataFrameBuilder(self.file).DataFrame_getn(index)
                
                TB = Tracks_Builder(df)
                if not TB.Checkgoodevent():
                    continue
                
                Figure_plot(ax,df)
                good_cnt += 1
        
            
        selected_events = self.ev[self.ev['EVENT_NUMBER'].isin(rand_lst)]
        Raw_data_plots(selected_events).Hit_Matrix_plot()
            
        return
            
        











def Figure_plot(ax,df):
    # chamber's dimensions
    width = c.WIDTH
    height = c.HEIGHT


    # chambers global coordinates
    x0 = c.X_SHIFT[0]
    z0 = c.Z_SHIFT[0] - c.ZCELL/2
    
    x1 = c.X_SHIFT[1]
    z1 = c.Z_SHIFT[1] - c.ZCELL/2
    
    x2 = c.X_SHIFT[2]
    z2 = c.Z_SHIFT[2] - c.ZCELL/2
    
    x3 = c.X_SHIFT[3]
    z3 = c.Z_SHIFT[3] - c.ZCELL/2
    
    # function that take single hit's coordinates and define a patch for the 
    # corresponding cell in which signal has been found
    def BuildCell(hit):
        wire_x = (hit.xleft_g + hit.xright_g)/2
        xc = wire_x - c.XCELL/2
        zc = hit.z_g - c.ZCELL/2
        cell = patches.Rectangle((xc,zc), c.XCELL, c.ZCELL, fill = False, color = "lightgrey")
        return cell


    # get hit's cell and print it
    def StampCell(event,axe):
        cell = BuildCell(event)
        axe.add_patch(cell)
        return
    

    def fitting_func(x,par):
        return par[1]+x*par[0]
                
    
    det0 = patches.Rectangle((x0,z0), width, height, fill=False, color = "lightgrey")
    det1 = patches.Rectangle((x1,z1), width, height, fill=False, color = "lightgrey")
    det2 = patches.Rectangle((x2,z2), width, height, fill=False, color = "lightgrey")
    det3 = patches.Rectangle((x3,z3), width, height, fill=False, color = "lightgrey")        
        
    
    # plotting chamber positions
    ax.add_patch(det0)
    ax.add_patch(det1)
    ax.add_patch(det2)
    ax.add_patch(det3)
            
    # drawing options
    plt.xlim([-1000,1000])
    plt.ylim([-100,1000])
    plt.xlabel("x [mm]")
    plt.ylabel("z [mm]")
    plt.title('Many events')
    
    results = Tracks_Builder(df).Get_Results()
    
    tracks = []
    if Tracks_Builder(df).checksx:
        tracks.append('sx')
    if Tracks_Builder(df).checkdx:
        tracks.append('dx')
            
    for side in tracks:
        par = results[side][0].beta
        x=np.linspace(-1000,1000,100)
        plt.plot(x,fitting_func(x,par),color = "black")
    
    # drawing points
#    plt.scatter(df['xleft_g'],df['z_g'],c='blue',label='left_hits')
#    plt.scatter(df['xright_g'],df['z_g'],c='red',label='right_hits')
        
    # drawing cells
#    for event in df.itertuples():
#        StampCell(event,ax)
    
    return
    
    
    
    
    
    
    
    
