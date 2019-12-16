#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 18:44:10 2019

@author: mattia
"""



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
import constants as c
import seaborn as sns

# x = f(y)
def f(par,y):
    m = par[0]
    q = par[1]
    return (y-q)/m


class Analysis_class:
    
    def __init__(self,visualize,save):
        
        self.good_event_counter = 0
        self.bad_event_chamber = {'0':0,'1':0,'2':0,'3':0}
        self.total_events = 0
        
        self.residuals_stored = {0:[],1:[],2:[],3:[]}
        self.residuals_2D_stored = {'res':[],'dist_wire':[]}
        self.angles_stored = {'xz':[],'yz':[]}
        
        self.visualize = visualize
        self.save = save
        
        
        
    def Compute_Residuals(self,dataframe,track_par):        
        
        
        sides = ['xz','yz']        
        chambers = {'xz':[0,2],'yz':[1,3]}
        
        for side in sides:
            
            df = dataframe
            par = track_par[side][0]
            
            df['WIRE'] = (df['XLEFT_G']+df['XRIGHT_G'])/2
            df['HIT_X'] = 0
            df.loc[(df['CHAMBER'].isin(chambers[side])) & (df['XLEFT_G'].isin(track_par[side][1][:,0])),'HIT_X'] = df['XLEFT_G']
            df.loc[(df['CHAMBER'].isin(chambers[side])) & (df['XRIGHT_G'].isin(track_par[side][1][:,0])),'HIT_X'] = df['XRIGHT_G']
            df = df.loc[df['HIT_X']!=0]
            
            df['FIT_X'] = f(par,df['Z_G'])
            df['HIT-WIRE'] = abs(df['HIT_X']-df['WIRE'])
            df['FIT-WIRE'] = abs(df['FIT_X']-df['WIRE'])
            df['RESIDUAL'] = df['FIT-WIRE'] - df['HIT-WIRE']
            
            df = df.loc[abs(df['RESIDUAL']) < (c.XCELL/2)]
            
            for ch in chambers[side]:
                res_ch = df.loc[df['CHAMBER']==ch,'RESIDUAL']
                self.residuals_stored[ch].append(res_ch)
            
            self.residuals_2D_stored['res'].append(df['RESIDUAL'])
            self.residuals_2D_stored['dist_wire'].append(df['HIT-WIRE'])
            
        return
    
    
    def Compute_angles(self,track_par):
        
        sides = ['xz','yz']   
        
        for side in sides:
            slope = track_par[side][0][0]
            angle = math.pi/2 - math.atan2(1,1/slope)
            self.angles_stored[side].append(angle)
        
        return 
        
    
    def update_good_event_counter(self):
        self.good_event_counter += 1
        return
    
    def update_rejected_counter(self,TB):
        if TB.GoodChamber == {0:False,1:True,2:True,3:True}: self.bad_event_chamber['0']+=1
        if TB.GoodChamber == {0:True,1:False,2:True,3:True}: self.bad_event_chamber['1']+=1
        if TB.GoodChamber == {0:True,1:True,2:False,3:True}: self.bad_event_chamber['2']+=1
        if TB.GoodChamber == {0:True,1:True,2:True,3:False}: self.bad_event_chamber['3']+=1
        return
    
    
    def update_analysis_object(self,df,TB,stats=False,angles=False,residuals=False):
        
        self.total_events +=1        
        if not TB.GoodEvent:
            self.update_rejected_counter(TB)
            return
        par = TB.dict_global_parameters
        if stats: self.update_good_event_counter()
        if angles: self.Compute_angles(par)
        if residuals: self.Compute_Residuals(df,par)
        return
    
    
    
    def Residuals_plot_1D(self,bins=100):
        
        fig = plt.figure()
        subplot_index = [3,4,1,2]

        for ch in range(4):

            residuals = np.concatenate(self.residuals_stored[ch])
            plt.hist(residuals,bins,density=True)
            
            plt.subplot(2,2,subplot_index[ch])
            
            mean = np.mean(residuals)
            variance = np.var(residuals)
            sigma = np.sqrt(variance)
            x = np.linspace(min(residuals), max(residuals),bins)
            plt.plot(x, stats.norm.pdf(x, mean, sigma),label='Gauss fit')
            plt.title('Chamber ' + str(ch))
            plt.xlabel('Residuals (mm)')
            
            self.Write('Gaussian fit parameters of residuals 1D for chamber ' + str(ch) + ': ')
            self.Write('Mu: ' + str(mean))
            self.Write('Sigma: ' + str(sigma))
        
        fig.tight_layout()            
        
        if self.visualize:
                plt.show()
                plt.close()
        if self.save:
            plt.savefig('Residuals_1D.png')
        
        return
    
    
    def Residuals_plot_2D(self,bins=100):
        
        sns.set(style="white", color_codes=True)

        residuals = np.array(np.concatenate(self.residuals_2D_stored['res']),dtype=float)
        dist_hit_wire = np.array(np.concatenate(self.residuals_2D_stored['dist_wire']),dtype=float)
        scatter_points = pd.DataFrame({'res':residuals,'dist_wire':dist_hit_wire})        
        
        x,y,erry=np.array([]),np.array([]),np.array([])        
        for i in np.arange(0,20,1.5):
           x = np.append(x,i+1.5/2)
           y = np.append(y,scatter_points.loc[(scatter_points['dist_wire']>i) & (scatter_points['dist_wire']<=i+1.5),'res'].mean())
           erry = np.append(erry,scatter_points.loc[(scatter_points['dist_wire']>i) & (scatter_points['dist_wire']<=i+1.5),'res'].std())
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(dist_hit_wire,residuals)
        sns.jointplot(x=scatter_points['dist_wire'],y=scatter_points['res'],kind="reg", color="k", marker="+", line_kws={'label':"y={0:.3f}x+{1:.3f}".format(slope,intercept)})
        plt.errorbar(x, y, yerr=erry, fmt='o',color="r")   
        
        plt.xlabel('Distance from wire (mm)')
        plt.ylabel('Residuals (mm)')
        
        self.Write('Parameters of the 2D residuals linear interpolation:')
        self.Write('Slope: ' + str(slope))
        self.Write('Intercept: ' + str(intercept))
        
        if self.visualize:
                plt.show()
                plt.close()
        if self.save: plt.savefig('Residuals_2D.png')
        
        return
    
    
    def Angle_resolution_plot(self,bins=50):
        
        plt.figure()  
        titles = ['Chambers 0,2' , 'Chambers 1,3'] 
        sides = ['xz','yz']

        for j in range(2):   
            plt.subplot(1,2,j+1)
            angles = np.array(self.angles_stored[sides[j]])
            plt.hist(angles,bins,density=True)
            plt.title(titles[j])
            plt.xlabel('Angle respect to z axis (rad)')
            
            mean = np.mean(angles)
            variance = np.var(angles)
            sigma = np.sqrt(variance)
            x = np.linspace(min(angles), max(angles),bins)
            plt.plot(x, stats.norm.pdf(x, mean, sigma),label='Gauss fit')
            
            self.Write('Gaussian fit parameters of angle resolution for side ' + sides[j] + ': ')
            self.Write('Mu: ' + str(mean))
            self.Write('Sigma: ' + str(sigma))
            
        if self.visualize:
                plt.show()
                plt.close()
        if self.save: plt.savefig('Angles.png')
        
        return
    
    
    def Good_events_informations(self):
        
        self.Write('Total events analized: ' + str(self.total_events))
        self.Write('Good events selected: ' + str(self.good_event_counter))
        self.Write('Percentage of events judged good: ' + str((self.good_event_counter/self.total_events)*100) + '%')
        for ch in range(4):
            self.Write('Events rejected (only) by bad chamber ' + str(ch) + ': ' + str(self.bad_event_chamber[str(ch)]) + '(' + str((self.bad_event_chamber[str(ch)]/self.total_events)*100) + '%)')
        return
    
    
    def Write(self,output):
        if self.visualize: print(output + '\n')
        if self.save:
            output_file = open('Analysis_parameters.txt',"a")
            output_file.write(output + '\n')
            output_file.close()
        return
    
    
        
        
        