#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 22:31:51 2019

@author: mattia
"""




import pandas as pd
import numpy as np
from itertools import combinations
from copy import copy

import constants as c
from pattern import meantimer_implementation

import time	
from tqdm import tqdm



class Trigger:
    
    def __init__(self,trigger):
        
        self.trigger = trigger
        self.events = pd.DataFrame()
                        
    
    def Complete_Dataframe(self,df):
                   
        # TIME
        # There is a problem with the precision of the measures, so we drop the orbit
        # Real time: data['TIME_NS'] = data["ORBIT_CNT"]*3564*25 + data["BX_COUNTER"]*25 + data["TDC_MEAS"]*25/30
        df['TIME_NS'] = df["ORBIT_CNT"]*3564*25 + df["BX_COUNTER"]*25 + df["TDC_MEAS"]*25/30
        
        # LAYER
        # To get the layer we must get the remainder of the TDC_CHANNEL with 4
        df['LAYER'] = df['TDC_CHANNEL']%4
        # Map 1 -> 4
        df.loc[df['LAYER'] == 1, 'LAYER'] = 4
        # Map 4 -> 1
        df.loc[df['LAYER'] == 0, 'LAYER'] = 1
        
        # CHAMBER
        df['CHAMBER'] = 0  # empty column
        # Detector 0
        df.loc[(df['FPGA']==0) & (df['TDC_CHANNEL']<=64),'CHAMBER'] = 0
        # Detector 1
        df.loc[(df['FPGA']==0) & (df['TDC_CHANNEL']>64) & (df['TDC_CHANNEL']<=128),'CHAMBER'] = 1
        # Detector 2
        df.loc[(df['FPGA']==1) & (df['TDC_CHANNEL']<=64),'CHAMBER'] = 2
        # Detector 3
        df.loc[(df['FPGA']==1) & (df['TDC_CHANNEL']>64) & (df['TDC_CHANNEL']<=128),'CHAMBER'] = 3
           
        # CELL
        df['CELL'] = ((df['TDC_CHANNEL']%64)/4).apply(np.ceil).astype(int)
        # TDC_CHANNEL%64=0 refers always to cell 16 of layer 1
        df.loc[df['CELL']==0,'CELL'] = 16
        
        # GLOBAL TDC CHANNEL
        df['TDC_CHANNEL_GLOB'] = (df['TDC_CHANNEL'] - 64 * (df['CHAMBER']%2)).astype(np.uint8)

        
        # Adding columns to be calculated
        nHits = df.shape[0]
        df['t0'] = np.zeros(nHits, dtype=np.float64)
        df['EVENT_NUMBER'] = np.ones(nHits, dtype=np.uint32) * -1
        
        return df
        
    
    def Create_Events(self,df):
                
        # Silence warning
        pd.options.mode.chained_assignment = None  # default = 'warn'
            
        # Grouping hits separated by large time gaps together
        events = df.sort_values('TIME_NS')
        grp = events['TIME_NS'].diff().fillna(0)
        events['TIME_GAP']=grp
        grp[grp <= 1.1*c.TMAX] = 0
        grp[grp > 0] = 1
        grp = grp.cumsum().astype(np.int32)
        events['EVENT_NUMBER'] = grp
                        
        if self.trigger != 0:
            
            # TDetecting events through trigger 139 (mean-time trigger) or 137-138 (Scintillator trigger)
            orbit = events.loc[df['TDC_CHANNEL']==self.trigger,'ORBIT_CNT']
            list_orbit = orbit.values.tolist()
            events = events.loc[df['ORBIT_CNT'].isin(list_orbit)]
            
        # Remove the hits related to the external trigger
        events = events[events['TDC_CHANNEL']<129]
        
        # Values will be sort according to their EVENT_NUMBER, CHAMBER and LAYER
        events.drop(['TIME_GAP'],axis=1,inplace=True)
        events = events.sort_values(by = ['EVENT_NUMBER','CHAMBER','LAYER'])   
        
        return events
    
    
    def Compute_time_pedestal(self,events):
                
        #Selecting only events with at least 3 hits
        gr_ev = events.groupby('EVENT_NUMBER')
        rejected_ev = gr_ev.size()
        rejected_ev = rejected_ev[rejected_ev < 6]
        ev_sel = events['EVENT_NUMBER'].isin(rejected_ev.index)
        events.loc[ev_sel, 'EVENT_NUMBER'] = -1
        events.drop(events.index[events['EVENT_NUMBER'] == -1], inplace=True) 
        events_accepted = []

        print('Computing time pedestal...')
        for event, df in tqdm(events.groupby('EVENT_NUMBER')):

            n_layers = df.groupby('CHAMBER')['LAYER'].agg('nunique')
            
            # Skipping events that have less than 2 chambers with at least 3 layers
            if n_layers[n_layers >= 3].shape[0] < 2: continue
        
            # Calculating numbers of hits in each chamber
            grp = df.groupby('CHAMBER')['TDC_CHANNEL_GLOB']
            nHits = grp.agg('nunique')

            tzeros_all = {}
            tzero_ch = 0
            # Starting from SLs with smallest N of hits
            chamber_indexes = nHits.loc[n_layers >= 3].sort_values().index
            for iCH, CH in enumerate(chamber_indexes):
                tzeros_all[CH],time_lst = meantimer_implementation(df[df['CHAMBER'] == CH])
                if len(tzeros_all[CH]) != 0:
                    tzero_ch = sum(tzeros_all[CH])/len(tzeros_all[CH])
                    
                    # Updating the TIME0 with meantimer result
                    pattern_time = (events['EVENT_NUMBER']==event) & (events['CHAMBER']==CH) & (events['TIME_NS'].isin(time_lst))
                    events.loc[pattern_time,'t0'] = tzero_ch  
            
            # Accepting the event
            events_accepted.append(event)

        events = events[events['EVENT_NUMBER'].isin(events_accepted)]
        
        return events
    
    
    def Compute_Position(self,events):
        
        # Compute the position
        events.loc[events['t0']!=0,'POSITION'] = (events['TIME_NS'] - events['t0'])*c.Vd
        events = events.fillna(0)
        # Now we have to filter false alignmens detected by the mask
        # pattern, which can be caused by noise
        # So we drop the rows with an unexpected result: such as Posi
        # tion bigger than 21 mm or smaller than 0 mm
        events.loc[(events['POSITION']<0) | (events['POSITION']>=21),['POSITION','t0']] = 0
        
        if self.trigger == 0:
            events = events[events['POSITION']!=0]
        
        
                
        return events

    
    
    def Get_Final_Dataframe(self,data):
        
        if not self.events.empty:
            return self.events
        
        dataframe = self.Complete_Dataframe(data)
        events = self.Create_Events(dataframe)
        ev_t0 = self.Compute_time_pedestal(events)
        ev_t0_pos = self.Compute_Position(ev_t0)

        # sorting events and removing useless columns       
        events_final = ev_t0_pos[['EVENT_NUMBER','ORBIT_CNT','FPGA','CHAMBER','LAYER','CELL','TIME_NS','t0','POSITION']]
        events_final.set_index(['EVENT_NUMBER','ORBIT_CNT','CHAMBER','LAYER'], inplace=True)
        events_final.sort_index(inplace=True)
        events_final.reset_index(inplace=True)
        
        self.events = copy(events_final)
        
        return events_final
        

    
    
    
    
    
