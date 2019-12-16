#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 18:18:42 2019

@author: mattia
"""

import numpy as np
import pandas as pd

from scipy import odr, stats, interpolate
from itertools import combinations
from math import fabs

import constants as c


def fitting_func(par,x):
    m = par[0]
    q = par[1]
    return m*x+q      
        
    


class Tracks_Builder:
        
    def __init__(self,dataframe):
        
        self.dict_local_parameters = {}
        self.dict_global_parameters = {}
        self.vector3D = ()
        
        self.df = dataframe
        
        self.GoodChamber = {}        
        self.GoodEvent = self.Checkgoodevent() 
        
        self.Reco_tracks()
                    
        
    
# check at least 3 hits in 3 different layer of every chamber (for the local linear fit)
# check if hits are in the known regions
    def CheckChamber(self,nch):
        
        if nch in self.GoodChamber: return self.GoodChamber[nch]
        
        pd.options.mode.chained_assignment = None  # default='warn'
        
        ch_dataframe = self.df[self.df['CHAMBER']==nch]
        ch_dataframe.reset_index(inplace=True,drop=True)

#        bkg_cnt = 0
                
#        for hit in ch_dataframe.itertuples():
#            
#           # don't want more than 2 bkg singals outside ranges
#           if (hit.xleft_l < c.REG[nch][0]) or (hit.xright_l > c.REG[nch][1]):
#                bkg_cnt+=1
#           if bkg_cnt>1:
#                    return False
        
        # removing events with no enough points or with too much noise
        if (ch_dataframe.shape[0] < 3) or (ch_dataframe.shape[0] > 7):
            self.GoodChamber.update( {nch:False} )
            return False
         
        check_layers = [0,0,0,0]  
            
        # don't want more than 3 hits (6 points) for layer (1 signal + 2 bkg)
        for i in range(1,5):
            if i in ch_dataframe['LAYER'].values:
                check_layers[i-1] = 1
            if ch_dataframe[ch_dataframe['LAYER']==i].shape[0] > 3:
                self.GoodChamber.update( {nch:False} )
                return False
            
        # we need at least 3 signal from different layers to perform a local fit
        if sum(check_layers) < 3:
            self.GoodChamber.update( {nch:False} )
            return False

        self.GoodChamber.update( {nch:True} )
        
        return True
    
    

# function that check if an event can be considered good: i.e. if we can provide
# a local linear fit (with 3 points) in 2 consecutives chambers
    def Checkgoodevent(self):
        
        # Need at least two chambers with information in both x and y axis
        self.CheckChamber(0)
        self.CheckChamber(1)
        self.CheckChamber(2)
        self.CheckChamber(3)
        
#        check_x = self.GoodChamber[0] or self.GoodChamber[2]
#        check_y = self.GoodChamber[1] or self.GoodChamber[3]
#        
#        if check_x and check_y:
#            if self.Secondary_check_fit: return True
        
        if self.GoodChamber == {0:True,1:True,2:True,3:True}: return True
            
        return False           
        
        

# function to implement another method to check if a series of hits is a signal event
# checking at the level of the local fit, asking for interpolations with no
# eccessive slope            
    def Secondary_check_fit(self):
    
        for ch in self.dict_local_parameters:
            
            # Preventing eccessive slope values
            local_slope = abs(self.Get_localfit_parameters(ch).beta[0])
            if c.Slope_limits[0] <= local_slope <= c.Slope_limits[1]: return False
            # Checking chi square
            if self.Get_localfit_parameters(ch)[0].res_var > c.Chi_square_max: return False
        
        return True    
    
    
    
    def Reco_tracks(self):
        
        if not self.GoodEvent: return     

        for ch in range(4):
            if self.GoodChamber[ch]: self.Get_localfit_parameters(ch)
        
        self.Get_globalfit_parameters('xz')
        self.Get_globalfit_parameters('yz')
        #self.Global_track_3D()
        self.vector3D = self.Reco_track_projection()
        
        return 
        
        
        
    def BestLocalFit(self,ch):
        
        ch_dataframe = self.df[self.df['CHAMBER']==ch]
        ch_dataframe.reset_index(inplace=True,drop=True)
        
        nhits = 2*ch_dataframe.shape[0]
                
        R = 999999
        best_out = None
        
        l_arr = np.array(ch_dataframe[['XLEFT_G','Z_G']])
        r_arr = np.array(ch_dataframe[['XRIGHT_G','Z_G']])
        arr = np.vstack((l_arr,r_arr))
        
        
        # to do fit for groups of 3 and 4 points (when is possible)  
        if set([1,2,3,4]).issubset(ch_dataframe['LAYER'].tolist()):
            k = 4
        else:
            k = 3
                                                                    
        p = np.empty(shape=(k,2))
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
                best_out = out
                best_hits = np.array(p, copy=True)  
            
                    
            else:
                continue            

        results = (best_out,best_hits)
 
        return results
    
    
    
    def Reco_track_projection(self):
        
        # XZ section
        proj_xz = self.dict_global_parameters['xz'][0]
        proj_yz = self.dict_global_parameters['yz'][0]
        
        a = proj_xz[0]
        b = proj_xz[1]
        c = proj_yz[0]
        d = proj_yz[1]
        
        y = np.linspace(-1000,1000,100)
        x = (c/a)*y -(b/a) + (d/a)
        z = c*y + d
    
        
        return(x,y,z)
    
    
    
    
    
    
#    def Global_track_3D(self):
#        
#        hits = []
#        
#        for ch in [0,2]:
#            if self.GoodChamber[ch]:
#                hits_ch = pd.DataFrame(self.dict_local_parameters[ch][1],columns = ['x','z'])
#                hits_ch['y'] = c.LENGHT/2
#                hits.append(hits_ch)
#        
#        for ch in [1,3]:
#            if self.GoodChamber[ch]:
#                hits_ch = pd.DataFrame(self.dict_local_parameters[ch][1],columns = ['y','z'])
#                hits_ch['x'] = c.LENGHT/2
#                hits.append(hits_ch)
#            
#        points = pd.concat(hits,ignore_index=True, copy=False)
#    
#        return
   
    
    
    
    def Global_fit(self,side):

        if side == 'xz':   
            hits = np.vstack((self.dict_local_parameters[0][1],self.dict_local_parameters[2][1]))
        elif side == 'yz':
            hits = np.vstack((self.dict_local_parameters[1][1],self.dict_local_parameters[3][1]))
        
        slope,intercept,_,_,_ = stats.linregress(hits[:,0],hits[:,1])      
        self.dict_global_parameters.update( {side:((slope,intercept),hits)} )
        
        return ((slope,intercept),hits)

            
                
    def Get_localfit_parameters(self,chamber):
        
        if chamber in self.dict_local_parameters:
            return self.dict_local_parameters[chamber]
        
        else:
            self.dict_local_parameters.update( {chamber:self.BestLocalFit(chamber)} )
            return self.dict_local_parameters[chamber]
    
    
    def Get_globalfit_parameters(self,side):
        
        if side in self.dict_global_parameters:
            return self.dict_global_parameters[side]
        
        else:
            fit_results = self.Global_fit(side)
            self.dict_global_parameters.update( {side:fit_results} )
            return self.dict_global_parameters[side]
        

   
    def Get_stuff_for_residuals(self):
        
        parameters = self.dict_global_parameters
        points_xz = parameters['xz'][1]
        points_yz = parameters['yz'][1]
        
        self.df['WIRE'] = (self.df['XLEFT_G']+self.df['XRIGHT_G'])/2
        self.df['FITTED'] = 0
        
        self.df.loc[(self.df['CHAMBER'].isin([0,2])) & (self.df['XLEFT_G'].isin(points_xz[:,0])),'FITTED'] = self.df['XLEFT_G']
        self.df.loc[(self.df['CHAMBER'].isin([0,2])) & (self.df['XRIGHT_G'].isin(points_xz[:,0])),'FITTED'] = self.df['XRIGHT_G']
        self.df.loc[(self.df['CHAMBER'].isin([1,3])) & (self.df['XLEFT_G'].isin(points_yz[:,0])),'FITTED'] = self.df['XLEFT_G']
        self.df.loc[(self.df['CHAMBER'].isin([1,3])) & (self.df['XRIGHT_G'].isin(points_yz[:,0])),'FITTED'] = self.df['XRIGHT_G']
        self.df = self.df.loc[self.df['FITTED']!=0]
        
        df_xz = self.df[self.df['CHAMBER'].isin([0,2])]
        df_yz = self.df[self.df['CHAMBER'].isin([1,3])]
        
        arr_xz = np.array((df_xz['FITTED'],df_xz['Z_G'],df_xz['WIRE']))
        arr_yz = np.array((points_yz[:,0],points_yz[:,1],df_yz['WIRE']))
        arr_xz = np.transpose(arr_xz)
        arr_yz = np.transpose(arr_yz)
        
        final_dict = {'par_xz':parameters['xz'][0],'par_yz':parameters['yz'][0],'xz':arr_xz,'yz':arr_yz}
        
        return final_dict
    
    
    
    
            
        
        

