#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 14:07:36 2019

@author: mattia
"""

"""CELL PATTERNS FOR MEANTIMER PATTERN MATCHING
patterns  =  dictionary of patterns (list of triplets of hits time in bx relative to orbit: BX + TDC/30), arranged by key being a string identifier used to select the proper mean-timing eq to be solved 
"""

import itertools
import numpy as np
import pandas as pd
import constants as c


############################################## MEANTIMER EQUATIONS

def meantimer_implementation(df):
    
    # Getting a TIME column as a Series with TDC_CHANNEL_GLOB as index
    df_time = df.loc[:, ['TDC_CHANNEL_GLOB', 'TIME_NS', 'LAYER']]
    df_time.sort_values('TIME_NS', inplace=True)
    # Split hits in groups where time difference is larger than maximum event duration
    grp = df_time['TIME_NS'].diff().fillna(0)
    grp[grp <= 1.1*c.TMAX] = 0
    grp[grp > 0] = 1
    grp = grp.cumsum().astype(np.uint16)
    df_time['grp'] = grp
    #Removing groups with less than 3 unique hits
    df_time = df_time[df_time.groupby('grp')['TDC_CHANNEL_GLOB'].transform('nunique') >= 3]

    # Determining the TIME0 using triplets [no external trigger]
    tzeros = []
    # Processing each group of hits
    patterns = PATTERN_NAMES.keys()
    time_lst = []
    for grp, df_grp in df_time.groupby('grp'):
        df_grp.set_index('TDC_CHANNEL_GLOB', inplace=True)
        # Selecting only triplets present among physically meaningful hit patterns
        channels = set(df_grp.index.astype(np.int16))
        triplets = set(itertools.permutations(channels, 3))
        triplets = triplets.intersection(patterns)

        # Grouping hits by the channel for quick triplet retrieval
        times = df_grp.groupby(df_grp.index)['TIME_NS']

        # Analysing each triplet
        for triplet in triplets:
            triplet_times = [times.get_group(ch).values for ch in triplet]
            for t1 in triplet_times[0]:
                for t2 in triplet_times[1]:
                    for t3 in triplet_times[2]:
                        timetriplet = (t1, t2, t3)
                        if max(timetriplet) - min(timetriplet) > 1.1*c.TMAX:
                            continue
                        pattern = PATTERN_NAMES[triplet]
                        mean_time = meantimereq(pattern, timetriplet)
                        tzeros.append(mean_time)
                        for j in range(3):
                            time_lst.append(timetriplet[j])
    return tzeros,time_lst



def meantimereq(pattern, timelist):
    """Function returning the expected t0 out of hits triples. None by default 
    tkey is the identifier of the univoque equation to be used given the pattern of hits in the triplet
    timelist is a len=3 list of hits time
    """
    patt = pattern[:-1]
    if patt in ('ABC','BCD'): 
        tzero = 0.25 * ( timelist[0] + 2*timelist[1] + timelist[2] - 2*c.TMAX)
    elif patt == 'ABD':
        tzero = 0.25 * ( 3*timelist[1] + 2*timelist[0] - timelist[2] - 2*c.TMAX)
    elif patt == 'ACD':
        tzero = 0.25 * ( 3*timelist[1] + 2*timelist[2] - timelist[0] - 2*c.TMAX)
    else:
        return None

    return tzero





NCHANNELS = 64
############################################# POSSIBLE HIT PATTERNS
PATTERNS = {}
### 3 ABC RIGHT
PATTERNS['ABCr']  = [ (1+x, 3+x,  2+x) for x in range(0, NCHANNELS, 4) ]
#A |1   x    |5   o    |9   o    |
#B     |3   x    |7   o    |
#C |2   x    |6   o    |10  o    |
#D     |4   o    |8   o    |
### 3 ABC LEFT
PATTERNS['ABCl'] = [ (5+x, 3+x,  6+x) for x in range(0, NCHANNELS, 4)[:-1] ]
#A |1   o    |5   x    |9   o    |
#B     |3   x    |7   o    |
#C |2   o    |6   x    |10  o    |
#D     |4   o    |8   o    |

### 3 BCD RIGHT
PATTERNS['BCDr']  = [ (3+x, 6+x,  4+x) for x in range(0, NCHANNELS, 4)[:-1] ]
#A |1   o    |5   o    |9   o    |
#B     |3   x    |7   o    |
#C |2   o    |6   x    |10  o    |
#D     |4   x    |8   o    |
### 3 BCD LEFT
PATTERNS['BCDl'] = [ (3+x, 2+x,  4+x) for x in range(0, NCHANNELS, 4) ]
#A |1   o    |5   o    |9   o    |
#B     |3   x    |7   o    |
#C |2   x    |6   o    |10  o    |
#D     |4   x    |8   o    |

### 3 ACD RIGHT
PATTERNS['ACDr']  = [ (1+x, 2+x,  4+x) for x in range(0, NCHANNELS, 4) ]
#A |1   x    |5   o    |9   o    |
#B     |3   o    |7   o    |
#C |2   x    |6   o    |10  o    |
#D     |4   x    |8   o    |
### 3 ACD LEFT
PATTERNS['ACDl'] = [ (5+x, 6+x,  4+x) for x in range(0, NCHANNELS, 4)[:-1] ]
#A |1   o    |5   x    |9   o    |
#B     |3   o    |7   o    |
#C |2   o    |6   x    |10  o    |
#D     |4   x    |8   o    |

### 3 ABD RIGHT
PATTERNS['ABDr']  = [ (1+x, 3+x,  4+x) for x in range(0, NCHANNELS, 4) ]
#A |1   x    |5   o    |9   o    |
#B     |3   x    |7   o    |
#C |2   o    |6   o    |10  o    |
#D     |4   x    |8   o    |
### 3 ABD LEFT
PATTERNS['ABDl'] = [ (5+x, 3+x,  4+x) for x in range(0, NCHANNELS, 4)[:-1] ]
#A |1   o    |5   x    |9   o    |
#B     |3   x    |7   o    |
#C |2   o    |6   o    |10  o    |
#D     |4   x    |8   o    |

# Transposed dictionary to quickly find pattern name
PATTERN_NAMES = {}
for name, patterns in PATTERNS.items():
    for pattern in patterns:
        PATTERN_NAMES[pattern] = name

# Lists of channels where hits from good events are expected AUGUST
ACCEPTANCE_CHANNELS = {
    0: range(1, NCHANNELS+1),
    1: range(1, NCHANNELS+1),
    2: range(1, NCHANNELS+1),
    3: range(1, NCHANNELS+1),
}