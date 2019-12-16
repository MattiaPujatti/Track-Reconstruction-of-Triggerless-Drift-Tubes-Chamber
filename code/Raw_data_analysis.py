#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 13:57:09 2019

@author: mattia
"""




import pandas as pd
import matplotlib.pyplot as plt
import time
from tqdm import tqdm

from Raw_data_plot import Raw_data_plots


def Good_Chamber(df,ch):
    
    # Good is a chamber with at least 3 hits in 3 different layers
    if ch not in df['CHAMBER'].values:
        return False
    if df[df['CHAMBER']==ch].shape[0] < 3:
        return False
    if len(df.loc[df['CHAMBER']==ch,'LAYER'].unique()) < 3:
        return False

    return True
        
        

class Raw_Data_Analysis:
    
    def __init__(self,data,events,stamp,save):
        
        self.data = data
        self.events = events
        self.stamp = stamp
        self.save = save

    def Trigger_efficiency(self,verbose=False):    

        print('Calculating selection efficiency')        
        
        total_events = len(set((self.events['EVENT_NUMBER'])))
        self.Write('\n')
        self.Write('Total events identified: ' + str(total_events))
        
        # Printing avaiable triggers for the analysis (the used one is known, i hope)
        trigger_used = set()
        trigger_used = self.data.loc[self.data['TDC_CHANNEL']>129,'TDC_CHANNEL'].unique()

        self.Write('External triggers avaiable: ' + str(list(trigger_used)))
        
        
        hits_selected = self.events.shape[0]
        hits_selected_percentage = (hits_selected/self.data[self.data['TDC_CHANNEL']<129].shape[0])*100
        
        # Visualize all the hit from orbits marked by trigger 
        self.Write('Total hits detected (signal + background): '+ str(self.data[self.data['TDC_CHANNEL']<129].shape[0]))
        self.Write('Number of hits belonging to orbit marked by selected trigger: ' + str(hits_selected) + ' ' + '(' + str(hits_selected_percentage) + '%)')        
        
        hits_with_a_position = self.events[self.events['POSITION']!=0].shape[0]
        hits_with_a_position_percentage = (hits_with_a_position/self.data[self.data['TDC_CHANNEL']<129].shape[0])*100
        
        # Visualize only hit with a legit calculated position
        self.Write('Number of hits which show a legit calculated position: ' + str(hits_with_a_position) + ' ' + '(' + str(hits_with_a_position_percentage) + '%)')
        
        
        # Computing number of 'showers'
        
        # Group the events by event number
        grouped_ev_num = self.events.groupby('EVENT_NUMBER').size().to_frame('SIZE')
        grouped_ev_num.reset_index(inplace=True)
        
        shower_twenty = grouped_ev_num[grouped_ev_num['SIZE']>20].shape[0]
        shower_twenty_percentage = (shower_twenty/grouped_ev_num.shape[0])*100
        shower_forty = grouped_ev_num[grouped_ev_num['SIZE']>40].shape[0]
        shower_forty_percentage = (shower_forty/grouped_ev_num.shape[0])*100
        
        self.Write('Number of events made by more than 20 hits: ' + str(shower_twenty) + ' ' + '(' + str(shower_twenty_percentage) + '%)')
        self.Write('Number of events made by more than 40 hits: ' + str(shower_forty) + ' ' + '(' + str(shower_forty_percentage) + '%)')
        
        
        if not verbose:
            return
        
        reduced_events = self.events[['EVENT_NUMBER','CHAMBER','LAYER']]

        ch_cnt = [0,0,0,0]
        conditional_cnt_02,conditional_cnt_20 = 0,0
        conditional_cnt_13,conditional_cnt_31 = 0,0
        
        print('Computing chamber efficiency')
        for EVENT_NUMBER, temp_df in tqdm(reduced_events.groupby('EVENT_NUMBER')):
            
            # Counting good events
            if Good_Chamber(temp_df,0):
                ch_cnt[0] += 1
                if Good_Chamber(temp_df,2):
                    conditional_cnt_20 += 1
                    
            if Good_Chamber(temp_df,2):
                ch_cnt[2] += 1
                if Good_Chamber(temp_df,0):
                    conditional_cnt_02 += 1
                    
            if Good_Chamber(temp_df,1):
                ch_cnt[1] += 1
                if Good_Chamber(temp_df,3):
                    conditional_cnt_31 += 1
            
            if Good_Chamber(temp_df,3):
                ch_cnt[3] += 1
                if Good_Chamber(temp_df,1):
                    conditional_cnt_13 += 1
        
        self.Write('\n')
        self.Write('Conditional probability using neighbour chamber as a tag')
        self.Write('Good_ch0 | Good_ch2: ' + str((conditional_cnt_02/ch_cnt[2])*100) + '%')
        self.Write('Good_ch2 | Good_ch0: ' + str((conditional_cnt_20/ch_cnt[0])*100) + '%')
        self.Write('Good_ch1 | Good_ch3: ' + str((conditional_cnt_13/ch_cnt[3])*100) + '%')
        self.Write('Good_ch3 | Good_ch1: ' + str((conditional_cnt_31/ch_cnt[1])*100) + '%')
        
                
        return
    
    
    def Hits_frequency_in_detectors(self):
        
        print('Creating 2D matrix of selected hits...')
        Raw_data_plots(self.events).Hit_Matrix_plot(self.stamp,self.save)
        plt.pause(0.1)
        
        
    def Hits_with_position_frequency(self):
        
        print('Creating 2D matrix of detected hits which position can be calculated...')
        Raw_data_plots(self.events[self.events['POSITION']!=0]).Hit_Matrix_plot(self.stamp,self.save)
        plt.pause(0.1)
        
    
    def Drift_times(self,nbins=60):
        
        print('Computing drift times...')
        Raw_data_plots(self.events).Drift_times_plot(nbins,self.stamp,self.save)

        return
    
    
    def Hits_per_event(self):
        
        print('Computing number of hits per event...')
        Raw_data_plots(self.events).Hits_per_event_plot(self.stamp,self.save)
        
        return
    
    
    def Single_track(self,num):
        
        Raw_data_plots(self.events).Hit_Matrix_single_plot(num,self.stamp,self.save)
        
        return
    
    
    def Hits_per_chamber(self):
        
        print('Computing number of hits per chamber...')
        Raw_data_plots(self.events).Hits_per_chamber(self.stamp,self.save)
        
        return
    
    
    def Hits_per_FPGA(self):
        
        print('Computing number of hits per FPGA...')
        Raw_data_plots(self.events).Hits_per_FPGA(self.stamp,self.save)
        
        return
    
    
    def Write(self,output):
        
        if self.stamp: print(output + '\n')
        
        if self.save:
            output_file = open('Trigger_efficiency.txt',"a")
            output_file.write(output + '\n')
            output_file.close()
            
        return
        
        
    
        
        
        
        
        
        
        
