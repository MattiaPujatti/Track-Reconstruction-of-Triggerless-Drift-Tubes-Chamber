#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 18:18:42 2019

@author: mattia
"""

import numpy as np
from scipy import odr
from itertools import combinations

import constants as c


def fitting_func(par,x):
    m = par[0]
    q = par[1]
    return m*x+q


class DataFrameOperate:
        
    def __init__(self,dataframe):
        
        self.df = dataframe
        
        check0 = 0 in self.df.chamber.values
        check1 = 1 in self.df.chamber.values
        check2 = 2 in self.df.chamber.values
        check3 = 3 in self.df.chamber.values
        
        self.checksx = check0 and check1
        self.checkdx = check2 and check3
        
        self.dic_layer_cnt = {0:[0,0,0,0],1:[0,0,0,0],2:[0,0,0,0],3:[0,0,0,0]}
        
        self.GoodEvent = self.Checkgoodevent()
        
        self.dict_parameters = {}
            
        
    
# function that returns the number of particles that passed through detectors        
    def TracksDetect(self):
        
        tracks = 0
        
        if self.checksx and self.checkdx:
            tracks = 2
        elif self.checksx or self.checkdx:
            tracks = 1
        else:
            tracks = 0
            
        return tracks
    
    
# check at least 3 hits in 3 different layer of every chamber (for the local linear fit)
# check if hits are in the known regions
    def CheckChamber(self,nch):
        
        bkg_cnt = 0       
        
        for event in self.df.itertuples():
                                    
            if event.chamber == nch:
                
                # don't want more than 2 bkg singals outside ranges
                if (event.xleft_l < c.REG[nch][0]) or (event.xright_l > c.REG[nch][1]):
                    bkg_cnt+=1
                if bkg_cnt>1:
                    return False
                
                for i in range(4):
                    if event.layer == (i+1):
                        self.dic_layer_cnt[nch][i]+=1
                    else:
                        continue
                            
        # don't want more than 3 hits (6 points) for layer
        for j in range(4):
            if self.dic_layer_cnt[nch][j] > 3:
                return False
            
        # we need at least 3 signal from different layers to perform a local fit
        if self.dic_layer_cnt[nch].count(0)>1:
            return False
        
        return True
    

# function that check if an event can be considered good: i.e. if we can provide
# a local linear fit (with 3 points) in 2 consecutives chambers
    def Checkgoodevent(self):
        
        ntracks = self.TracksDetect()
        
        if ntracks == 0:
            return False
        
        elif ntracks == 1:
            
            if self.checksx:
                if self.CheckChamber(0) and self.CheckChamber(1):
                    return True
                else:
                    return False
            
            elif self.checkdx:
                if self.CheckChamber(2) and self.CheckChamber(3):
                    return True
                else:
                    return False
                
        elif ntracks == 2:
            
            check01 = self.CheckChamber(0) and self.CheckChamber(1)
            check23 = self.CheckChamber(2) and self.CheckChamber(3)
            
            checkdecay = check01 and check23
            
            if checkdecay:
                return True
            else:
                return False
            
    
    def BestLocalFit(self,ch):
        
        nhits = 2*sum(self.dic_layer_cnt[ch])
                
        R = 100000
        best_parameters = []
        
        arr = np.empty(shape=(nhits,2))

        i=0
        for event in self.df.itertuples():
                    
            if event.chamber == ch:
                arr[i] = [event.xleft_g,event.z_g]
                arr[i+1] = [event.xright_g,event.z_g]
                i+=2    
            else:
                continue
        
        
        # to do fit for groups of 3 and 4 points (when is possible) 
        if 0 in self.dic_layer_cnt[ch]:
            k = 3
        else:
            k = 4
                                                
        p = np.empty(shape=(k,2))
        best_hits = np.empty(shape=(k,2))
        lst = list(range(nhits))    
        
        for case in combinations(lst,k):
            a=0
            for index in case:
                p[a] = arr[index]
                a+=1

                                   
            # to prevent fit on points in the same layer (same z)
            cnt_z = [0 for t in range(k)]
            for j in range(k):
                cnt_z[j] = list(p[:,1]).count(p[j,1])
            if cnt_z != [1 for l in range(k)]:
                continue
                
            
            linear_model = odr.Model(fitting_func)
            data = odr.RealData(p[:,0],p[:,1],sx=c.XERR)
            odr1 = odr.ODR(data=data,model=linear_model,beta0=[0,0])
            out = odr1.run()
            reduced_chisquare = out.res_var
            
            if reduced_chisquare < R:
                R = reduced_chisquare
                best_parameters = out.beta
                
                h=0
                for index in case:
                    best_hits[h] = arr[index]
                    h+=1                
                    
            else:
                continue                          
        
        results = (best_parameters[0],best_parameters[1],best_hits)
 
        return results
   
    
    
    
    def Global_fit(self,side):

        if side == 'sx':
            if not 0 in self.dict_parameters:
                self.Get_localfit_parameters(0)
            if not 1 in self.dict_parameters:
                self.Get_localfit_parameters(1)
              
            hits = np.vstack((self.dict_parameters[0][2],self.dict_parameters[1][2]))
                
        elif side == 'dx':
            if not 2 in self.dict_parameters:
                self.Get_localfit_parameters(2)
            if not 3 in self.dict_parameters:
                self.Get_localfit_parameters(3)

            hits = np.vstack((self.dict_parameters[2][2],self.dict_parameters[3][2]))
                
        linear_model = odr.Model(fitting_func)
        data = odr.RealData(hits[:,0],hits[:,1],sx=c.XERR)
        odr1 = odr.ODR(data=data,model=linear_model,beta0=[0,0])
        out = odr1.run()
        
        return out.beta

            
                
    def Get_localfit_parameters(self,chamber):
        
        if chamber in self.dict_parameters:
            return self.dict_parameters[chamber]
        
        else:
            fit_results = self.BestLocalFit(chamber)
            self.dict_parameters.update( {chamber:fit_results} )
            return self.dict_parameters[chamber]
    
    
    def Get_globalfit_parameters(self,side):
        
        if side in self.dict_parameters:
            return self.dict_parameters[side]
        
        else:
            fit_results = self.Global_fit(side)
            self.dict_parameters.update( {side:fit_results} )
            return self.dict_parameters[side]
        

    def Get_Results(self):
        
        if self.checksx:
            self.Get_globalfit_parameters('sx')
        
        if self.checkdx:
            self.Get_globalfit_parameters('dx')
            
        return
            
        
        

