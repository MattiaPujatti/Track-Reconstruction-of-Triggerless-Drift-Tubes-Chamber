#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 18:34:10 2019

@author: mattia
"""



import mpl_toolkits.mplot3d.art3d as art3d
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


import constants as c

# chambers global coordinates
chambers_height = c.Z_SHIFT
x0 = 0
z0 = c.Z_SHIFT[0]
x1 = 0
z1 = c.Z_SHIFT[1]
x2 = 0
z2 = c.Z_SHIFT[2]
x3 = 0
z3 = c.Z_SHIFT[3]


# function that take single hit's coordinates and define a patch for the 
# corresponding cell in which signal has been found
def BuildCell(hit):
    wire_x = (hit.XLEFT_G + hit.XRIGHT_G)/2
    xc = wire_x - c.XCELL/2
    zc = hit.Z_G - c.ZCELL/2
    cell = patches.Rectangle((xc,zc), c.XCELL, c.ZCELL, fill = False, color = "lightgrey")
    return cell


# get hit's cell and print it
def StampCell(hit,ax):
    cell = BuildCell(hit)
    ax.add_patch(cell)
    return
    

def fitting_func(x,par):
    return par[1]+x*par[0]




class Event_Plot:
    
    def __init__(self,dataframe,TB,stamp=False,save=False):
        self.df = dataframe
        self.tb = TB
        self.visualize = stamp
        self.save = save
        

    def Plot_2D(self,detectors_zoom=False,bkgplot=False):
      
        # if the event is not good, it won't be plotted if bkgplot is not set True
        # bkgplot is False by default to prevent useless plots in cycles
        Good = True 
        if not self.tb.GoodEvent:
            if bkgplot:
                Good = False
            else:
                return
        
        fig = plt.figure()  
        
        det0 = patches.Rectangle((x0,z0), c.LENGHT, c.WIDTH, fill=False, color = "lightgrey")
        det1 = patches.Rectangle((x1,z1), c.LENGHT, c.WIDTH, fill=False, color = "lightgrey")
        det2 = patches.Rectangle((x2,z2), c.LENGHT, c.WIDTH, fill=False, color = "lightgrey")
        det3 = patches.Rectangle((x3,z3), c.LENGHT, c.WIDTH, fill=False, color = "lightgrey")
        
        chs = [[0,2],[1,3]]
        detectors = [det0,det1,det2,det3]
        xlabel = ['x [mm]', 'y[mm]']
        side = ['xz','yz']
        titles = ['Chambers 0,2' , 'Chambers 1,3'] 
        
        for j in range(2):
            
            ax = plt.subplot(1,2,j+1)
        
            # plotting chamber positions
            ax.add_patch(detectors[chs[j][0]])
            ax.add_patch(detectors[chs[j][1]])
                
            # drawing points on chambers 0,2
            df = self.df[(self.df['CHAMBER']==chs[j][0]) | (self.df['CHAMBER']==chs[j][1])]
            ax.scatter(df['XLEFT_G'],df['Z_G'],c='blue', s=7 ,label='left_hits')
            ax.scatter(df['XRIGHT_G'],df['Z_G'],c='red', s=7, label='right_hits')
        
       
            # drawing cells
            for hit in df.itertuples():
                StampCell(hit,ax)  
        
        
            ax.set_xlim(-100,1000)
            ax.set_ylim(-100,1000)
            ax.set_xlabel(xlabel[j])
            ax.set_ylabel('z [mm]')
            ax.set_title(str(self.df.name) + ' - ' + titles[j])
            ax.set_xticks(list(i for i in range(-100,1000,25)),minor=True)
            ax.set_yticks(list(i for i in range(-100,1000,25)),minor=True)
            ax.tick_params(which='both',direction='in')
        
            fig.tight_layout()
        
            if Good:
                self.Local_fit_plot(chs[j])
                self.Global_fit_plot([side[j]])
                plt.legend()
            
        if self.visualize:
            plt.show()
            plt.close()
        if self.save: plt.savefig('track.png')
        
        # zoom to detectors as subplots
        if detectors_zoom:
            fig1 = plt.figure(1, figsize=(8,5))
            for n in range(4):
                self.Plot_detector(n)
                plt.scatter(self.df['XLEFT_G'],self.df['Z_G'],c='blue')
                plt.scatter(self.df['XRIGHT_G'],self.df['Z_G'],c='red')
                if Good:
                    self.Local_fit_plot()
                    self.Global_fit_plot()
            fig1.tight_layout()
            
            if self.visualize:
                plt.show()
                plt.close()
            if self.save: plt.savefig('detectors.png')
        
           
        return
    
    
    def Plot_3D(self,projection=False):
        

        fig = plt.figure()
        ax=fig.gca(projection='3d')

        # vertices of a chamber
        v0 = np.array([[0,0,0],[0,700,0],[700,0,0],[700,700,0],[0,0,52],[0,700,52],[700,0,52],[700,700,52]])
        v1 = np.array([[0,0,62],[0,700,62],[700,0,62],[700,700,62],[0,0,114],[0,700,114],[700,0,114],[700,700,114]])
        v2 = np.array([[0,0,782],[0,700,782],[700,0,782],[700,700,782],[0,0,834],[0,700,834],[700,0,834],[700,700,834]])
        v3 = np.array([[0,0,844],[0,700,844],[700,0,844],[700,700,844],[0,0,896],[0,700,896],[700,0,896],[700,700,896]])        

        chamber_verts = [v0,v1,v2,v3]
        colors = ['black','grey','black','grey']
        
        i=0
        for v in chamber_verts:
            ch = [ [v[0],v[2]],[v[4],v[6]],[v[0],v[4]],[v[2],v[6]],[v[0],v[1]],[v[4],v[5]],[v[2],v[3]],[v[6],v[7]],[v[1],v[5]],[v[3],v[7]],[v[1],v[3]],[v[5],v[7]]]
            ax.add_collection3d(art3d.Poly3DCollection(ch,facecolors=colors[i], linewidths=1, edgecolors=colors[i], alpha=.25))
            i+=1
        

        # Getting vector of 3D recostruction
        xline,yline,zline = self.tb.vector3D
        ax.plot3D(xline, yline, zline, 'blue',label='Muon Track')
        
        # Printing 2D projection on axes3D
        if projection:
            
            # XZ
            par = self.tb.Get_globalfit_parameters('xz')[0]
            x_xz = np.linspace(-100,900,100)
            z_xz = fitting_func(x_xz,par)
            y_xz = np.full((x_xz.shape),-100)
            ax.plot3D(x_xz, y_xz, z_xz,'darkgrey',label='xz_proj')
            # YZ
            par = self.tb.Get_globalfit_parameters('yz')[0]
            y_yz = np.linspace(-100,900,100)
            z_yz = fitting_func(x_xz,par)
            x_yz = np.full((x_xz.shape),-100)
            ax.plot3D(x_yz, y_yz, z_yz,'black',label='yz_proj')
        
        ax.set_title(str(self.df.name))         
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_xlim(-100,1000)
        ax.set_ylim(-100,1000)
        ax.set_zlim(-100,1000)
        plt.legend()
        
        if self.visualize:
                plt.show()
                plt.close()
        if self.save: plt.savefig('3D_track.png')
        
        return
    
    
    
    
    
    
    
    def Plot_detector(self,n):
        
        # for events in the middle of the detectors        
        
        df_ch = self.df[self.df['CHAMBER']==n]
        limits = [0,0]
        
        if df_ch.empty: limits = [-100,1000]
        else:
            limits[0] = min(df_ch['XLEFT_G']) - c.XCELL
            limits[1] = max(df_ch['XRIGHT_G']) + c.XCELL
        
         
        subplot_index = [3,4,1,2]
        
        det=patches.Rectangle((c.X_SHIFT[n],c.Z_SHIFT[n]), c.LENGHT, c.WIDTH, fill=False, color = "lightgrey")
    
        # building subplot of the chamber n
        ax = plt.subplot(2,2,subplot_index[n])
        
        for hit in self.df.itertuples():
            StampCell(hit,ax)
            
        ax.add_patch(det)
        plt.xlim(limits[0],limits[1])
        plt.ylim([c.Z_SHIFT[n]-2*c.ZCELL,c.Z_SHIFT[n]+c.WIDTH+2*c.ZCELL])
        if n in [0,2]: plt.xlabel("x [mm]")
        if n in [1,3]: plt.xlabel("y [mm]")
        plt.ylabel("z [mm]")
        plt.title('Chamber ' + str(n))
        
        return 
    
        
    def Local_fit_plot(self,which=[]):
        
        if which:
            chambers = which
        else: chambers = [0,1,2,3]
                
        for ch in chambers:
            
            if self.tb.GoodChamber[ch]:
                fit_results = self.tb.Get_localfit_parameters(ch)[0].beta
                   
            x=np.linspace(-1000,1000,100)
            plt.plot(x,fitting_func(x,fit_results),color = "lightblue",linestyle='--',label = 'local_fit')
                
        return 
    
    
    
    
    def Global_fit_plot(self,which=[]):
             
        if which:
            tracks = which
        else: tracks = ['xz','yz']
        
        for side in tracks:
            par = self.tb.Get_globalfit_parameters(side)[0]
            
            x=np.linspace(-1000,1000,100)
            plt.plot(x,fitting_func(x,par),color = "black",label = 'projection')
                
        return
 
    
    
    
    
