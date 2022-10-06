#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 11:00:39 2022

@author: bas
"""

from heatmaps import getLogZ

models = ['shock','kilonova']
datas = ['shock','kilonova']
delay = 0
dist = 40
optical_band='r'
uv_band='no_uv'

log10zs = []
for model in models:
    for data in datas:        
        log10z, dlogz, model, data = getLogZ(model,data,dist,optical_band=optical_band,uv_band=uv_band,delay=delay)
        print("model:",model, ", data:",data)
        print(log10z)
        log10zs.append(log10z)

# shock model is strongly ruled out for the kilonova data in r

# log10bs = []
# for index in [0,1,2]:
#     log10b = log10zs[index]-log10zs[index+3]
#     log10bs.append(log10b)