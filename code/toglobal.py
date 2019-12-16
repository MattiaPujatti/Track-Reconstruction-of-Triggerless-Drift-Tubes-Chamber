#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 23:43:47 2019

@author: mattia
"""

import constants as c

xs0 = c.X_SHIFT[0]
xs1 = c.X_SHIFT[1]
xs2 = c.X_SHIFT[2]
xs3 = c.X_SHIFT[3]

zs0 = c.Z_SHIFT[0]
zs1 = c.Z_SHIFT[1]
zs2 = c.Z_SHIFT[2]
zs3 = c.Z_SHIFT[3]


# function to compute the change to global coordinates from detector's one
def To_global(value, coo, det):
    
    if coo == 'x':
        if det == 0:
            return value + xs0
        elif det == 1:
            return value + xs1
        elif det == 2:
            return value + xs2
        elif det == 3:
            return value + xs3
    elif coo == 'z':
        if det == 0:
            return value + zs0
        elif det == 1:
            return value + zs1
        elif det == 2:
            return value + zs2
        elif det == 3:
            return value + zs3         
    
        
    