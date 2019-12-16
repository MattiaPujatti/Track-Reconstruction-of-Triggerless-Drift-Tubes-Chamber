#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 17:24:10 2019

@author: mattia
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 13:23:01 2019

@author: mattia
"""

import pandas as pd

import constants as c

# function to compute the change to global coordinates from detector's one
def To_global(value, coo, det):

    if coo == 'x':
        for chamber in range(4):
            if det == chamber: return (value + c.X_SHIFT[int(det)])
        
    elif coo == 'z':
        for chamber in range(4):
            if det == chamber: return (value + c.Z_SHIFT[int(det)])
            



class DataFrameBuilder:
    
    def __init__(self,datafile):
         self.datafile = datafile
       
        
    def GetData(self,line):

            # create a list of float with single event data
            event_lst = [float(i) for i in line.split()]

            # if there are not particles in the event, it is useless to continue the analisys    
            if ( event_lst[1] == 0 ):
                return []
            else:
                return event_lst
            
            
    def DF_builder(self,line):
            
        event_lst = self.GetData(line)

        # return an empty dataframe if there are no datas in the line
        if (event_lst == []):
            df = pd.DataFrame(data=None)
            return df
        
        nhits = int(event_lst[1])
        coordinates = ['CHAMBER', 'LAYER', 'XLEFT_L', 'X_RIGHT_L', 'DRIFT_TIME', 'XLEFT_G', 'XRIGHT_G', 'Z_G']
        hits = [i+1 for i in range(int(nhits))]

        # create a dataframe for every event with at least one particle's signal    
        df = pd.DataFrame(data=None, index=hits, columns=coordinates)
        df.name = 'Event ' + str(int(event_lst[0]))
        
        try:
            # analogue to do while in C++
            i,j = 1,2
            while True:
            
                chamber = event_lst[j]
                layer =   event_lst[j+1]
                time =    event_lst[j+4]
                xll = event_lst[j+2]
                xrl = event_lst[j+3]
                zl  = (layer-0.5)*c.ZCELL    #local z
            
                xlg = To_global(xll, 'x',chamber)
                xrg = To_global(xrl, 'x',chamber)
                zg  = To_global(zl, 'z',chamber)
                            
                df.loc[i] = [chamber,layer,xll,xrl,time,xlg,xrg,zg]
                i += 1
                j += 5
                if(j >= (2 + 5*nhits)):
                    break
        except: return pd.DataFrame(data=None)
        
        return df
    

        
        
# function that get a single event (by a line) and return the corresponding dataframe              
    def DataFrame_getl(self,line,stamp=False):
        
        df = self.DF_builder(line)
        if stamp:
            if not df.empty:
                print(df.name)
                print(df)
        return df
        
# function that get a single event (by its number) and return the corresponding dataframe                
    def DataFrame_getn(self,evn,stamp=False):
        
        found = False
        self.datafile.seek(0)
        for l in self.datafile:
            event_lst = self.GetData(l)
            if event_lst == []:
                continue                    
            else:
                if (event_lst[0] == evn):
                    df = self.DF_builder(l)
                    found = True
                else:
                    continue
        
        if found:
            if stamp:
                if not df.empty:
                    print(df.name)
                    print(df)
            return df
        else:
            return pd.DataFrame(data=None)
        
   


           
                
                




