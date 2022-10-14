#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 11:59:37 2021

@author: basdorsman
"""

import pickle
import numpy as np
from heatmaps import getLogZ

# uv vs no uv given r band kilonova default and shock data at 40 mpc
# models = ['kilonova','kilonova','shock', 'shock']
# read_data = 'kilonova_optical'
# dist = 40
# include_optical = ['r']
# include_uvs = ['no_uv','uv', 'no_uv','uv'] #uv, no_uv
# timestamps = ['21-08-18 0854','21-08-17 1956','21-08-19 1142', '21-08-19 1019'] #g+r
# timestamps = ['21-08-27 1012','21-08-30 1032','21-08-27 0949', '21-08-27 0956'] #g+r
# bs_string = ''.join(include_optical)
# resultsstrings = [f"{model}model_{read_data}data_{dist}Mpc_{bs_string}band_{include_uv}/{timestamp} results" for model, include_uv, timestamp in zip(models, include_uvs, timestamps)]


# confusion matrices at 160 mpc
# models = ['shock','kilonova','kilonova_uvboost','shock','kilonova','shock','kilonova_uvboost']
# read_datas = ['shock','shock','shock','kilonova_optical','kilonova_optical','kilonova_uvboost_optical','kilonova_uvboost_optical']
# include_optical = 'no_optical'
# dist = 160
# include_uv = 'uv'
# timestamps = ['21-10-04 1135','21-10-06 0828','21-10-04 1529','21-10-06 1417','21-10-04 1224','21-10-06 1437','21-10-04 1443']
# bs_string = ''.join(include_optical)
# resultsstrings = [f"{model}model_{read_data}data_{dist}Mpc_{bs_string}band_{include_uv}/{timestamp} results" for model, read_data,  timestamp in zip(models, read_datas, timestamps)]

# leftover kn uvboost vs kn models
# models = ['kilonova','kilonova_uvboost','kilonova','kilonova_uvboost','kilonova','kilonova_uvboost']
# read_datas = ['kilonova_uvboost_optical','kilonova_optical','kilonova_uvboost_optical','kilonova_optical','kilonova_uvboost_optical','kilonova_optical']
# include_optical = 'no_optical'
# dists = [40, 40, 100, 100, 160, 160]
# include_uv = 'uv'
# timestamps = ['21-10-07 1018','21-10-08 1002','21-10-07 1856','21-10-20 0942','21-10-07 2018','21-10-08 0915']
# bs_string = ''.join(include_optical)
# resultsstrings = [f"{model}model_{read_data}data_{dist}Mpc_{bs_string}band_{include_uv}/{timestamp} results" for model, read_data, dist, timestamp in zip(models, read_datas, dists, timestamps)]

# for resultsstring in resultsstrings:
#     with open(resultsstring,'rb') as analysis_results: 
#         results = pickle.load(analysis_results)
#     logz = results.logz[-1]
#     log10z = logz*np.log10(np.e)
#     print()
#     print(resultsstring)
#     results.summary()
#     print('log10z',log10z)

# redden vs not redden
model = 'kilonova'
data = 'kilonova'
dist = 40
reddens = [True, False]

evidences = [getLogZ(model,data,dist,optical_band='no_optical', uv_band='NUV_D', delay=0, redden=redden) for redden in reddens]
print(evidences)
