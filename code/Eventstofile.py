#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 22:24:00 2019

@author: mattia
"""

import pandas as pd

import constants as c


def Write_Events_tofile(dataframe,filename):
    
    output_file = open(filename,"w+")
        
    for EVENT_NUMBER, temp_df in dataframe.groupby('EVENT_NUMBER'):
        
        event = []
        
        event.append(EVENT_NUMBER)
        event.append(temp_df.shape[0])
        
        for hit in temp_df.itertuples():
            
            if hit.POSITION == 0:
                continue
                        
            event.append(hit.CHAMBER)
            layer = hit.LAYER
            event.append(layer)
            
            # computing xleft and xright
            xleft = c.XCELL*(hit.CELL-1) + c.XCELL/2 - hit.POSITION
            xright = c.XCELL*(hit.CELL-1) + c.XCELL/2 + hit.POSITION
            
            # layers 1 and 3 starts from half cell from the left of the detector
            if layer == 1 or layer == 3:
                xleft += c.XCELL/2
                xright += c.XCELL/2
            
            event.append(xleft)
            event.append(xright)
            drift_time = hit.TIME_NS - hit.t0
            event.append(drift_time)
            
        for value in event:
            output_file.write(str(value) + ' ')
        output_file.write('\n')
        
        
    output_file.close()    
    
    return


