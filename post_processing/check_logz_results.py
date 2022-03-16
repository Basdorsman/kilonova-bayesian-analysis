#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 11:00:39 2022

@author: bas
"""

from heatmaps import getLogZ

models = ['kilonova','shock']
data = 'shock'
delays = [0, 12, 24]
dist = 40

log10zs = []
for model in models:
    for delay in delays:
        log10z, dlogz, model, data = getLogZ(model,data,dist,delay=delay)
        print(log10z)
        log10zs.append(log10z)
        
log10bs = []
for index in [0,1,2]:
    log10b = log10zs[index]-log10zs[index+3]
    log10bs.append(log10b)