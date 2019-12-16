#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 10:01:36 2019

@author: mattia
"""



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns

from matplotlib.ticker import MultipleLocator


class Raw_data_plots:
    
    def __init__(self,events):
        
        self.events = events  
        self.max_ev = events['EVENT_NUMBER'].max()
    
    
    def Drift_times_plot(self,nbins=60,plot=False,save=False):
        
        # Select hits where position is different from 0
        events_drift = self.events[self.events['POSITION']!=0]
        # Get all drift times
        measured_times = np.array(events_drift['TIME_NS'])
        tzeros = np.array(events_drift['t0'])
        drift_times = np.array(events_drift['TIME_NS']-events_drift['t0'])
        
        Times = [measured_times,tzeros,drift_times]
        titles = ["Measured times Frequency","Time pedestal distribution","Drift Times Frequency"]
        labels = ["Measured times (s)","Time pedestal (ns)","Drift time (ns)"]
        
        for i in range(len(Times)):
            
            # Drift Times Frequency
            y, edges, bins = plt.hist(Times[i], bins = nbins, label='Drift Times Frequency', alpha=0.6,color = 'blue')
            plt.ylabel("Number of samples / Bin")
            plt.xlabel(labels[i])
            plt.title(titles[i], fontsize=20)

            if plot: 
                plt.show()
                plt.close()
                
            elif save: plt.savefig(titles[i] + '.png')
        
        return 
    
    
    def Hits_per_event_plot(self,plot=False,save=False):
        
        # Group the events by event number
        grouped_ev_num = self.events.groupby('EVENT_NUMBER')
        
        # For every event count the number of hits (i.e. return the lenght of the array containing the hits of a certain event)
        hits_per_event = []
        for i in list(grouped_ev_num.groups.keys()):
                hits_per_event.append(len(np.array(grouped_ev_num.groups[i])))
        hits_per_event = np.array(hits_per_event)

        # Plot the distribution of Hits/Event
        fig = plt.figure(figsize=(10,6))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlim(0,max(hits_per_event))
        ax.set_xticks = [i for i in range(max(hits_per_event))]
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.set_xlabel('Number of hits')
        ax.set_ylabel('Number of Events')
        
        ax.set_title("Hits Per Event Histogram", fontsize=20, verticalalignment='bottom')
        
        ax.hist(hits_per_event, bins=max(hits_per_event));
        
        if plot:
            plt.show()
            plt.close()
        elif save: plt.savefig('Hits_per_event.png')
        
        return
    
    
    def Hit_Matrix_plot(self,plot=False,save=False):
        
        df = self.events
    
        # Hit Matrix initialization
        matrix = np.zeros((8,32)) 
        final_matrix = np.zeros((8,66))
        
        grpd = df.groupby(['CHAMBER','LAYER','CELL']).size().to_frame('SIZE')
        grpd.reset_index(inplace=True)  
        
        grpd = pd.pivot_table(grpd,index=['CHAMBER','LAYER','CELL'],values='SIZE',fill_value = 0,dropna=False,aggfunc=np.sum)
        grpd.reset_index(inplace=True)  
        
        for i in range(4):
            
            mask_pattern_chamber = (grpd['CHAMBER'] == 1) | (grpd['CHAMBER'] == 3)
            mask_pattern_layer = (grpd['LAYER']==4-i)
            mask_pattern = mask_pattern_chamber & mask_pattern_layer
            matrix[i] = (grpd.loc[mask_pattern,'SIZE'])

        for i in range(4,8):
            mask_pattern_chamber = (grpd['CHAMBER'] == 0) | (grpd['CHAMBER'] == 2)
            mask_pattern_layer = (grpd['LAYER']==8-i)
            mask_pattern = mask_pattern_chamber & mask_pattern_layer
            matrix[i] = (grpd.loc[mask_pattern,'SIZE'])

                
        # Reshape the hit matrix
        for row in range(8):
            if (row%2)==0: final_index = 0
            elif (row%2)!=0: final_index = 1
            for index in range(16):
                final_matrix[row][final_index] = matrix[row][index]
                final_matrix[row][final_index+1] = matrix[row][index]
                final_index+=2
            final_index += 1
            for index in range(16,32):
                final_matrix[row][final_index] = matrix[row][index]
                final_matrix[row][final_index+1] = matrix[row][index]
                final_index+=2
                

                
        #Show the results
        plt.imshow(final_matrix, cmap='plasma')
        
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
        
        if plot:
            plt.show()
            plt.close()
        if save: plt.savefig('Hit_matrix.png')
             
        return 
    
    
    
    def Hit_Matrix_single_plot(self,eN,plot=False,save=False):
        
        self.events = self.events[self.events['EVENT_NUMBER']==eN]
        self.Hit_Matrix_plot() 
        plt.pause(0.1)
        
        return
  
    
    def Hits_per_chamber(self,plot=False,save=False):
        
        df = self.events
        
        ch0, ch1, ch2, ch3 = [], [], [], []
                
        ch_dict = {0:ch0,1:ch1,2:ch2,3:ch3}
                
        # Count the hits associated to each chamber for every event  
        grpd = df.groupby(['EVENT_NUMBER','CHAMBER']).size().to_frame('SIZE')
        grpd.reset_index(inplace=True)
        max_nhits = max(grpd['SIZE'])
        
        grpd['CHAMBER0'] = grpd.loc[grpd['CHAMBER']==0,'SIZE']  
        grpd['CHAMBER1'] = grpd.loc[grpd['CHAMBER']==1,'SIZE'] 
        grpd['CHAMBER2'] = grpd.loc[grpd['CHAMBER']==2,'SIZE'] 
        grpd['CHAMBER3'] = grpd.loc[grpd['CHAMBER']==3,'SIZE']
        grpd.fillna(0,inplace=True)
        grpd2 = grpd.groupby('EVENT_NUMBER').sum()
        
        ch_dict[0] = grpd2['CHAMBER0']
        ch_dict[1] = grpd2['CHAMBER1']
        ch_dict[2] = grpd2['CHAMBER2']
        ch_dict[3] = grpd2['CHAMBER3']
        
        
               
        # Plotting hits ch1 vs ch0
        plt.figure()         
        plt.hist2d(np.array(ch_dict[0]),np.array(ch_dict[1]),bins=np.arange(max_nhits+1)-0.5,norm=LogNorm())
        plt.xlabel('Hits ch0')
        plt.ylabel('Hits ch1')
        plt.title('Hits ch1 vs ch0')
        plt.colorbar()
        
        if plot:
            plt.show()
            plt.close()
        if save: plt.savefig('Hits_ch1_vs_ch0.png')
        
        # Plotting hits ch3 vs ch2
        plt.figure()         
        plt.hist2d(np.array(ch_dict[2]),np.array(ch_dict[3]),bins=np.arange(max_nhits+1)-0.5,norm=LogNorm())
        plt.xlabel('Hits ch2')
        plt.ylabel('Hits ch3')
        plt.title('Hits ch3 vs ch2')
        plt.colorbar()
        
        if plot:
            plt.show()
            plt.close()
        if save: plt.savefig('Hits_ch3_vs_ch2.png')
            
        return
    
    
    
    def Hits_per_FPGA(self,plot=False,save=False):
        
        df = self.events

        fpga0, fpga1 = [],[]
        
        dict_fpga = {0:fpga0,1:fpga1}
        
        # Count the hits associated to each chamber for every event  
        grpd = df.groupby(['EVENT_NUMBER','FPGA']).size().to_frame('SIZE')
        grpd.reset_index(inplace=True)
        
        max_nhits = max(grpd['SIZE'])
        
        grpd['FPGA0'] = grpd.loc[grpd['FPGA']==0,'SIZE']
        grpd['FPGA1'] = grpd.loc[grpd['FPGA']==1,'SIZE']
        grpd.fillna(0,inplace=True)
        grpd2 = grpd.groupby('EVENT_NUMBER').sum()
        dict_fpga[0] = grpd2['FPGA0']
        dict_fpga[1] = grpd2['FPGA1']
        
        
        print('Number of zero in FPGA0: ' + str(len(dict_fpga[0][dict_fpga[0] == 0])))
        print('Number of zero in FPGA1: ' + str(len(dict_fpga[0][dict_fpga[1] == 0])))
        
        # Plotting hits ch1 vs ch0
        plt.figure()         
        plt.hist2d(np.array(dict_fpga[1]),np.array(dict_fpga[0]),bins=np.arange(max_nhits+1)-0.5,norm=LogNorm(),range=[[0,max_nhits+1],[0,max_nhits+1]])
        plt.xlabel('Hits FPGA1')
        plt.ylabel('Hits FPGA0')
        plt.title('Hits FPGA0 vs FPGA1')
        plt.colorbar()
        
        if plot:
            plt.show()
            plt.close()
        if save: plt.savefig('Hits_FPGA0_vs_FPGA1.png')
    
        return
    



    
    
    
            
        