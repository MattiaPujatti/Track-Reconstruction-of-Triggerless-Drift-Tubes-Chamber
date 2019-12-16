#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 10:01:36 2019

@author: mattia
"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.ticker import MultipleLocator

class Raw_data_plots:
    
    def __init__(self,events):
        
        self.events = events  
        self.max_ev = events['EVENT_NUMBER'].max()
    
    
    def Drift_times_plot(self,nbins=60):
        
        # Select hits where position is different from 0
        events_drift = self.events[self.events['POSITION']!=0]
        # Get all drift times
        drift_times = events_drift['TIME_NS']-events_drift['t0']
        # Drift Times Frequency
        figure = plt.figure(figsize=(14,8))
        ax = figure.add_subplot(111)
        y, edges, bins = ax.hist(drift_times, bins = nbins, label='Drift Times Frequency', alpha=0.6)
        ax.set_ylabel("Number of samples / Bin")
        ax.set_xlabel("Drift time")
        ax.set_title("Drift Times Frequency", fontsize=20)
        mean_point = (edges[1:] + edges[:-1])/2
        ax.errorbar(mean_point, y, yerr = y**-0.5, marker = '.', drawstyle = 'steps-mid', label = 'error')
        
        plt.show()
        
        return 
    
    
    def Hits_per_event_plot(self):
        
        # Group the events by event number
        grouped_ev_num = self.events.groupby('EVENT_NUMBER')
        
        # For every event count the number of hits (i.e. return the lenght of the array containing the hits of a certain event)
        hits_per_event = []
        for i in grouped_ev_num.groups.keys():
                hits_per_event.append(len(np.array(grouped_ev_num.groups[i])))
        hits_per_event = np.array(hits_per_event)
        
        # Plot the distribution of Hits/Event
        fig = plt.figure(figsize=(10,6))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlim(0,20)
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.set_xlabel('Number of hits')
        ax.set_ylabel('Number of Events')
        
        ax.set_title("Hits Per Event Histogram and Distribution", fontsize=20, verticalalignment='bottom')
        
        sns.set()
        sns.distplot(hits_per_event, bins=100, ax=ax, color = 'blue')
        
        return
    
    
    def Hit_Matrix_plot(self):
        
        df = self.events
    
        # Hit Matrix initialization
        raw_mat = np.zeros((16,32))
        rows = np.array((df['CHAMBER']-1)*4+df['LAYER']-1)
        columns = np.array((df['CELL']-1)*2)
        
        # Count the hits associated to each cell
        for i in range(len(rows)):
            raw_mat[rows[i],columns[i]] = raw_mat[rows[i],columns[i]]+1
            raw_mat[rows[i],columns[i]+1] = raw_mat[rows[i],columns[i]+1]+1
            
        # Reshape the hit matrix
        final_mat = np.zeros((8,66))
    
        # Chamber 1,2
        for i in range(8):
            if (i%2 == 0):
                final_mat[7-i,1:33] = raw_mat[i,:32]
            else:
                final_mat[7-i,:32] = raw_mat[i,:32]
            
            # Chamber 3,4
            for i in range(8,16):
                if (i%2 == 0):
                    final_mat[7-(i-8),34:66] = raw_mat[i,:32]
                else:
                    final_mat[7-(i-8),33:65] = raw_mat[i,:32]
                        
        # Show the results
        plt.figure(figsize=(12,4))
        plt.imshow(final_mat, cmap='plasma')
        
        plt.annotate('CHAMBER 0', xy=(0.5,0.5), xytext=(12,13), fontsize=12, color= 'red')
        plt.annotate('CHAMBER 1', xy=(0.5,0.5), xytext=(12,-3), fontsize=12, color= 'red')
        plt.annotate('CHAMBER 2', xy=(0.5,0.5), xytext=(46,13), fontsize=12, color= 'red')
        plt.annotate('CHAMBER 3', xy=(0.5,0.5), xytext=(46,-3), fontsize=12, color= 'red')
               
        plt.axvline(x=32.5, color='white', linewidth=2)
        plt.axhline(y=3.5, color='white', linewidth=2)
        
        # Tick on the axis
        plt.xticks(np.concatenate((np.arange(1.5,33,2),np.arange(34.5,66,2))),
               ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16',
                '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'])
        plt.yticks([x for x in range(8)], ['4','3','2','1','4','3','2','1'])
    
        plt.xlabel('CELL',fontsize=10)
        plt.ylabel('LAYER',fontsize=10)
        plt.colorbar()
             
        return 
    
    
    
    def Hit_Matrix_single_plot(self,eN):
        
        self.events = self.events[self.events['EVENT_NUMBER']==eN]
        self.Hit_Matrix_plot() 
        plt.pause(0.1)
        
        return
  
    
    def Hits_per_chamber(self):
        
        df = self.events
        
        ch0 = np.empty([1,self.max_ev])
        ch1 = np.empty([1,self.max_ev])
        ch2 = np.empty([1,self.max_ev])
        ch3 = np.empty([1,self.max_ev])
        
        # Count the hits associated to each chamber for every event
        #for i in range(10):
            #print df[df['EVENT_NUMBER']==i]    
            #df1 = df[df['EVENT_NUMBER']==i]
        df1 = df.groupby(['EVENT_NUMBER','CHAMBER']).size()
        print df1['CHAMBER']
        
        return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
            
        