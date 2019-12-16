#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 22:30:22 2019

@author: mattia
"""

from __future__ import division

import matplotlib.pyplot as plt

from Trigger import Trigger
from Process_plots import Raw_data_plots


def Print_Useful_Stuff(data,events,trigger_class=0):
   
    plot_class = Raw_data_plots(events)
    
    # if events are not selected, will print only information about datas
    if trigger_class == 0:
        
        print 'Total hits detected (signal + background): ' + str(data[data['TDC_CHANNEL']<129].shape[0])
        
        trigger_used = set()
        for x in data[data['TDC_CHANNEL']>129]['TDC_CHANNEL']:
            trigger_used.add(x)
        print 'Triggers avaiable: ' + str(list(trigger_used))
        
        print '\n'
        print '2D matrix of detected hits: '
        plot_class.Hit_Matrix_plot()
        
        return
        
    
    ev_det = trigger_class.Get_detected_events()

    # Visualize all the hit from orbits marked by trigger 
    print 'Total hits detected: '+ str(data[data['TDC_CHANNEL']<129].shape[0])
    print 'Number of hits belonging to orbit marked by selected trigger: ' + str(events.shape[0]/data[data['TDC_CHANNEL']<129].shape[0]*100) + '%'
    print 'Number of events detected using the trigger: ' + str(ev_det)
    print '\n'
    print '2D matrix of detected hits: '
    plot_class.Hit_Matrix_plot()
    plt.pause(0.1)
    
    # Visualize only hit with a legit calculated position
    print 'Number of hit which show a legit calculated position: ' + str(events[events['POSITION']!=0].shape[0]/data[data['TDC_CHANNEL']<129].shape[0]*100) + '%'
    new_plot_class = Processing_Hits_Plots(events[events['POSITION']!=0])
    print '\n'
    print '2D matrix of detected hits which position can be calculated: '
    new_plot_class.Hit_Matrix_plot()
    plt.pause(0.1)
    
    print '\n'
    print 'Distribution of drift times of events detected'
    plot_class.Drift_times_plot(nbins=60)
    plt.pause(0.1)
    
    print '\n'
    print 'Histogram and distribution of the number of hits per event'
    plot_class.Hits_per_event_plot()
    plt.pause(0.1)
    
    
    return





