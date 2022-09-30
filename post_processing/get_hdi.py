#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 16:34:03 2021

@author: basdorsman
"""

import pickle
import matplotlib.pyplot as plt
import arviz as az

#HDI
with open("kilonovamodel_kilonova_opticaldata_40Mpc_rband_uv/21-08-30 1032 results",'rb') as analysis_results: 
    samples_40Mpc_uv = pickle.load(analysis_results).samples
    
with open("kilonovamodel_kilonova_opticaldata_40Mpc_rband_no_uv/21-08-27 1012 results",'rb') as analysis_results: 
#with open("shockmodel_shock_opticaldata_40Mpc_rband_uv/21-09-06 1201 results",'rb') as analysis_results: 
    samples_40Mpc_no_uv = pickle.load(analysis_results).samples
    

samples = [samples_40Mpc_uv, samples_40Mpc_no_uv]
parameters = ['M_ejecta','v_min','v_k','v_max','k_high','k_low','n']
hdis = [az.hdi(sample, hdi_prob=.95) for sample in samples]

hdis_dict = {}
for index in range(len(parameters)):
    hdis_dict[parameters[index]] = [hdi[index] for hdi in hdis]
    
for key in hdis_dict: 
    print(key)
    for index in range(len(hdis_dict[key])):
        print(f'interval {index}: {hdis_dict[key][index]}')
        print(f'interval length {hdis_dict[key][index][1]-hdis_dict[key][index][0]}')
    
    
