#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 11:14:22 2019

@author: mattia
"""

import numpy as np
from itertools import combinations
import constants as c
    

def Spot_Alignment(dataframe):
    
    n_points = dataframe.shape[0]
    t0 = 0
        
    # to consider all possible combinations of hits
    for case in combinations(list(range(n_points)),3):
        
        i=0
        possible_alignment = np.empty((3,3))
        for index in case:
            possible_alignment[i] = np.array([dataframe.iloc[index].LAYER,dataframe.iloc[index].CELL,dataframe.iloc[index].TIME_NS])
            i+=1
                        
        # Necessary conditions to apply Talete's theorem
            
        # Condition(1): layer patterns are 123 or 234
        check_123 = ([1,2,3] == possible_alignment[:,0].tolist())
        check_234 = ([2,3,4] == possible_alignment[:,0].tolist())
        align = check_123 or check_234
        if not align:
            continue
                
        # Condition(2): CELL_1 = CELL_3
        if possible_alignment[0][1] != possible_alignment[2][1]:
            continue
        
        # Condition(3): cells must differ at most by 1
        if abs(possible_alignment[0][1]-possible_alignment[1][1]) > 1:
            continue
        if abs(possible_alignment[1][1]-possible_alignment[2][1]) > 1:
            continue
                
        # Compute time pedestal with Talete's theorem
        t0 = (possible_alignment[0][2] + possible_alignment[2][2] + 2*possible_alignment[1][2] - 2*c.TMAX)/4
                
        return t0        
                
                
                
                