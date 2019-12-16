#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 13:58:00 2019

@author: mattia
"""

import itertools
import numpy as np
import pandas as pd
import constants as c
from pattern import PATTERN_NAMES, meantimereq


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


def mean_tzero(tzeros):
    """Calculates the most probable t0 from multiple candidates of different meantimer solutions"""
    df = tzero_clusters(tzeros)
    if df is None:
        return -1, [], 0
    # Selecting the largest cluster
    gr = df.groupby('cluster')
    nSols = gr.agg('size')
    nSLs = gr['ch'].agg('nunique')
    cluster_id = nSols.idxmax()
    cluster = df.loc[df['cluster'] == cluster_id, 't0'].values
    return cluster.mean(), cluster, nSLs[cluster_id]


# Difference between t0 candidates that should be clustered together for the mean [ns]
MEAN_TZERO_DIFF = 4

def tzero_clusters(tzeros):
    """Groups different tzero candidates into clusters"""
    clusters = []
    # Creating a dataframe with the data
    vals = []
    for ch, tz in tzeros.items():
        for tzero in tz:
            vals.append({'t0': tzero, 'ch': ch})
    if not vals:
        return None
    df = pd.DataFrame.from_dict(vals)
    df.sort_values('t0', inplace=True)
    # Adding a cluster column
    gr = df.sort_values('t0')['t0'].diff().fillna(0)
    gr[gr <= MEAN_TZERO_DIFF] = 0
    gr[gr > MEAN_TZERO_DIFF] = 1
    gr = gr.cumsum().astype(np.int8)
    df['cluster'] = gr
    # Counting number of layers and solutions in each cluster
    clusters = df.groupby('cluster')
    nSols = clusters.agg('size')
    nLayers = clusters['ch'].agg('nunique')
    clusters_valid = nSols[(nLayers >= 2) & (nSols >= 2)].index
    if len(clusters_valid) < 1:
        return None
    return df[df['cluster'].isin(clusters_valid)]














