#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 20:42:15 2019

@author: mattia
"""



import pandas as pd
import numpy as np

import sys
import os
import shutil
import random
import copy

from Trigger import Trigger
from Eventstofile import Write_Events_tofile
from Raw_data_analysis import Raw_Data_Analysis
from dataframe_builder import DataFrameBuilder
from tracks_building import Tracks_Builder
from eventplot import Event_Plot
from Analysis import Analysis_class
import argparse

import time
from tqdm import tqdm
start = time.time()

####################################################################################################################
parser = argparse.ArgumentParser(description='Pre-processing and processing of data from LEMMA used for cosmic muons')
parser.add_argument('-i',   '--input',         default=False,                      help='The unpacked input file to analyze')
parser.add_argument('-d',   '--directory',     default=False,                      help='Directory where data file are stored')
parser.add_argument('-pp',  '--preprocess',    default=False, action='store_true', help='Pre-processing raw data')
parser.add_argument('-ra',  '--raw_analysis',  default=False,                      help='Raw analysis for reordered hits file')
parser.add_argument('-ef',  '--events_file',   default=False, 			           help='Write events file')
parser.add_argument('-t',   '--trigger',       default=0,     type=int,            help='Full analysis of single event by EVENT_NUMBER')
parser.add_argument('-p',   '--process',       default=False, action='store_true', help='Process events')
parser.add_argument('-c',   '--calibration',   default=False, action='store_true', help='Estimate calibration parameters')
parser.add_argument('-one', '--single',        default=False,                      help='Full analysis of single event by EVENT_NUMBER')
parser.add_argument('-vb',  '--verbose',       default=False, action='store_true', help='For a complete (but longer) analysis')
parser.add_argument('-v',   '--visualize',     default=False, action='store_true', help='To show analysis plots run-time')
parser.add_argument('-s',   '--save',          default=False, action='store_true', help='Save analysis results')
parser.add_argument('-n',   '--events_number', default=-1,    type=int,            help="Reduce the number of the events to be processed or hits to be preprocessed")
args = parser.parse_args()
####################################################################################################################

import mmap
def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

def data_memory_optimize(df):
    
    # Removing possible incomplete rows e.g. last line of last file
    df.dropna(inplace=True)
    # Converting to memory-optimised data types
    for name in ['HEAD', 'FPGA', 'TDC_CHANNEL', 'TDC_MEAS']:
        df[name] = df[name].astype(np.uint8)
    for name in ['BX_COUNTER']:
        df[name] = df[name].astype(np.uint16)
    for name in ['ORBIT_CNT']:
        df[name] = df[name].astype(np.uint32)
        
    # retain all words with HEAD=1
    df.drop(df.index[df['HEAD'] != 1], inplace=True)
    # Removing unused columns to save memory foot-print
    df.drop('HEAD', axis=1, inplace=True)  
    
    return df

print('\n')
print('#############################################################################')
print('\n')



def Read_input():
    
    if args.input:
        print('Reading file '+ args.input + '...')
        if args.preprocess or args.raw_analysis:
            data = data_memory_optimize(pd.read_csv(args.input)) 
            return data
        if args.process:
            datafile = open(args.input)
            return datafile
    
    if args.directory:
        print('Reading file stored in directory ' + args.directory + '...')
        data_lst = []
        file_lst = [file for file in os.listdir(args.directory) if os.path.isfile(os.path.join(args.directory,file))]
        for i in tqdm(range(0,len(file_lst))):
            df = data_memory_optimize(pd.read_csv(args.directory + file_lst[i]))
            data_lst.append(df)
        data = pd.concat(data_lst, ignore_index=True, copy=False)   
        return data
           


def Preprocess(data):
               
    if (args.events_number != -1): data = data.head(args.events_number)
       
    print('Data successfully acquired!')
    print('Total hits acquired: ' + str(data.shape[0]))

    # Triggers avaiable: 139 (mean-time), 0 (triggerless)
    events = Trigger(args.trigger).Get_Final_Dataframe(data)

    if args.events_file: 
        print('Creating csv output file...')
        events.to_csv('events_dataframe_' + args.events_file + '.csv')
        print('Creating events output file (txt)...')
        Write_Events_tofile(events,args.events_file + '.txt')



def Raw_data_analysis(data):
    
    print('Starting raw data analysis...')       

    events = pd.read_csv(args.raw_analysis)    
    
    if args.save:
        owd = os.getcwd()
        if args.input: dir_name = str('Raw data analysis plots' + args.input.replace('/','_'))
        elif args.directory: dir_name = str('Raw data analysis plots' + args.directory.replace('/','_'))
        if os.path.exists(dir_name): shutil.rmtree(dir_name) 
        os.mkdir(dir_name)
        os.chdir(dir_name)
        
    Analysis_Class = Raw_Data_Analysis(data,events,args.visualize,args.save)
#    Analysis_Class.Trigger_efficiency(args.verbose)
#    Analysis_Class.Hits_frequency_in_detectors()
    Analysis_Class.Drift_times()
#    Analysis_Class.Hits_per_event()
#    Analysis_Class.Hits_per_chamber()
#    Analysis_Class.Hits_per_FPGA()

    if args.save: os.chdir(owd)
 

def Random_event_search(datafile):
    
    datafile.seek(0)
    rnd = random.randint(1,get_num_lines(args.input))
    df = DataFrameBuilder(datafile).DataFrame_getn(rnd,stamp=False)
    Good = Tracks_Builder(df).GoodEvent
    print(rnd,Good)
    if not Good:
        return Random_event_search(datafile)
    if Good:
        print('Found ' + str(df.name))
        return df


def Process(datafile,Analysis):
        
    if args.single:
        
        print('Searching for event ' + str(args.single) + '...')
        if args.single == 'random': df = Random_event_search(datafile)
        else: df = DataFrameBuilder(datafile).DataFrame_getn(int(args.single),stamp=True)
        if df.empty:
            print('Event not found!')
            sys.exit()
        
        TB = Tracks_Builder(df)
        
        if args.visualize or args.save:
            Event_Plot(df,TB,args.visualize,args.save).Plot_2D(detectors_zoom=True,bkgplot=True)   
            Event_Plot(df,TB,args.visualize,args.save).Plot_3D()
        sys.exit() 
        
    counter = 0
    if args.events_number>0: tot_ev = args.events_number
    else: tot_ev = get_num_lines(args.input)
    
    for line in tqdm(datafile,total=tot_ev):
        if args.events_number>0 and counter>=args.events_number: break
        counter += 1
        df = DataFrameBuilder(datafile).DataFrame_getl(line,stamp=False)
        if df.empty: continue
        TB = Tracks_Builder(df)
        Analysis.update_analysis_object(df,TB,stats=True,angles=True,residuals=True)
    
    return


def Analysis_plot(Analysis):
    
    if args.save:
        owd = os.getcwd()
        dir_name = str('Analysis plots' + args.input.replace('/','_'))
        if os.path.exists(dir_name): shutil.rmtree(dir_name) 
        os.mkdir(dir_name)
        os.chdir(dir_name)
        
    Analysis.Residuals_plot_1D()
    Analysis.Residuals_plot_2D()
    Analysis.Good_events_informations()
    Analysis.Angle_resolution_plot()
    
    if args.save: os.chdir(owd)
    
    return          

        
        


if args.preprocess:
    data = Read_input()
    Preprocess(data)
    
if args.raw_analysis: 
    data = Read_input()
    Raw_data_analysis(data)
    
if args.process:
    
    print('Data successfully acquired!')
    print('Start processing events...')
    
    datafile = Read_input()
    Analysis = Analysis_class(args.visualize,args.save)
    Process(datafile,Analysis)
    Analysis_plot(Analysis)
    
    
    
            
        


    
    
    
    
    
    
    




end = time.time()
print("Time: " + str(end - start))
