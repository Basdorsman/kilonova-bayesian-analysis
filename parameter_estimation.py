#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 09:28:48 2021

@author: basdorsman
"""

from dynesty import NestedSampler
import numpy as np
from produce_lightcurve import Lightcurve
import synphot as sp
import dorado.sensitivity
import astropy.units as u
from astropy import constants as c
import pickle 
import time
import os
from parameters import getParameters

# get parameters
parameters = getParameters(osargs_list=['read_data','model','delay','dist','include_optical','include_uv','print_progress','method','max_time','sample'])

model = parameters['model'] #shock, kilonova, kilonova_uvboost
delay = parameters['delay'] #hours
dist = parameters['dist'] #mpc
include_optical = parameters['include_optical'].split(',') # 'r', ['u', 'g','r', 'I', 'z'], ['False']
print_progress=parameters['print_progress']=='True'
method = parameters['method'] #'test', 'timeout', 'pool'
max_time = int(parameters['max_time']) # seconds, parameter for 'timeout' method
include_uv = parameters['include_uv'].split(',')
read_data = parameters['read_data']
sample_method = parameters['sample']

######## MORE PARAMETERS, DONT TOUCH ##########
distance = dist * u.Mpc
heating = 'beta'

if model == 'shock':
    radiation = 'shock' 
elif model == 'kilonova' or model == 'kilonova_uvboost':
    radiation = 'kilonova'


bs = {}
# Include uv bands
if include_uv == ['False']:
    uv_string = 'no_uv'
else:
    uv_string = ''.join(include_uv)
    if type(include_uv) == list:
        for key in include_uv:
            bs[key] = getattr(dorado.sensitivity.bandpasses, key)
    elif type(include_uv) == str:
        bs[include_uv] = getattr(dorado.sensitivity.bandpasses, include_uv)

# Include optical bands
if include_optical == ['False']:
    optical_string = 'no_optical'
else:
    optical_string = ''.join(include_optical)  
    if type(include_optical) == list:
        for key in include_optical:
            bs[key] = sp.SpectralElement.from_file(f'input_files/bands/SLOAN_SDSS.{key}.dat')
    elif type(include_optical) == str:
        bs[include_optical] = sp.SpectralElement.from_file(f'input_files/bands/SLOAN_SDSS.{include_optical}.dat')
     


#### READ DATA #########
with open(f'input_files/data/SNR_fiducial_{read_data}_{dist}Mpc_opticalbands_ugri_uvbands_NUV_DD1D2_{delay}h_delay.pkl','rb') as tf:
    data_list = pickle.load(tf)
ts_data, abmags_data, snrs, abmags_error = data_list



########## LOG PROBABILITIES ##########

if model == 'kilonova' or model == 'kilonova_uvboost':
    ndim = 7 
    limits = np.ndarray((ndim,2),dtype=object)
    limits[0] = [0.01, 0.1]  # mass (solar masses)
    if model == 'kilonova_uvboost':
        limits[1:4] = np.array(([0.05, 0.2], ['vmin', 'vmax'], [0.21, 0.8]),dtype=object)  # velocities (cm/s)
        limits[4:6] = np.array(([1,10],[0.01,0.1]))  # opacities (cm^2/g)
    elif model == 'kilonova':
        limits[1:4] = np.array(([0.05, 0.2], ['vmin', 'vmax'], [0.3, 0.8]),dtype=object)  # velocities (cm/s)
        limits[4:6] = np.array(([1,10],[0.1,1]))  # opacities (cm^2/g)
    limits[6] = [4,5]
elif model == 'shock':
    #limits=np.array(([0.01,10],[0.01,5],[0.01,3],[0.1,10])) #old broad
    limits=np.array(([1,10],[0.5,5],[1,3],[1,10])) #lim[[lower,upper],..
    # k in 0.1 cm^2/g, M in 0.01 solar masses, v in 0.1c, R_0 in 10^10 cm
    ndim = limits.shape[0]
    

########## DEFINE MODEL ##########

lightcurve_object = Lightcurve(distance, heating_function=heating)

if model == 'kilonova' or model == 'kilonova_uvboost':

    def priortransform(uniform): 
        mass = (limits[0,1]-limits[0,0])*uniform[0]+limits[0,0]
        #velocities = (limits[1:4,1]-limits[1:4,0])*uniform[1:4]+limits[1:4,0]
        v_min = (limits[1,1]-limits[1,0])*uniform[1]+limits[1,0]
        v_max = (limits[3,1]-limits[3,0])*uniform[2]+limits[3,0]
        v_k = (v_max-v_min)*uniform[3]+v_min
        #velocities = np.array([v_min, v_k, v_max])
        opacities = (limits[4:6,1]-limits[4:6,0])*uniform[4:6]+limits[4:6,0]
        n = (limits[6,1]-limits[6,0])*uniform[6]+limits[6,0]
        theta = np.array([mass, v_min, v_k, v_max, opacities[0], opacities[1], n]) #for postprocessing
        return theta    

    def lightcurve_model(t,theta_reshaped,bandpasses):
        abmags = lightcurve_object.calc_abmags_combined(t,theta_reshaped,bandpasses,radiation = radiation)
        return abmags

    def loglikelihood(theta):
        theta_reshaped = np.array((theta[0] * u.Msun, np.array((theta[1], theta[2], theta[3])) * c.c, np.array((theta[4], theta[5])) * u.cm**2 / u.g, theta[6]), dtype= 'object')
        sigmas2 = {}
        loglikelihood = {}
        abmags_model = lightcurve_model(ts_data, theta_reshaped, bs)
        for key in bs:
            sigmas2[key] = abmags_error[key]**2
            loglikelihood[key] = -0.5 * np.sum((abmags_data[key] - abmags_model[key]) ** 2 / sigmas2[key] + np.log(2*np.pi*sigmas2[key]))
        return sum(loglikelihood.values())
    
elif model == 'shock':

    def priortransform(u):
        theta = (limits[:,1]-limits[:,0])*u+limits[:,0]
        return theta
    
    def lightcurve_model(t, theta, bandpasses):
        abmags = lightcurve_object.calc_abmags_combined(t,theta,bandpasses,radiation=radiation)
        return abmags
    
    def loglikelihood(theta):
        sigmas2 = {}
        loglikelihood = {}
        abmags_model = lightcurve_model(ts_data, theta, bs)
        for key in bs:
            sigmas2[key] = abmags_error[key]**2
            loglikelihood[key] = -0.5 * np.sum((abmags_data[key] - abmags_model[key]) ** 2 / sigmas2[key] + np.log(2*np.pi*sigmas2[key]))
        return sum(loglikelihood.values())

####### TESTS / PARAMETER ESTIMATION #############
if method == 'test':
    print(f'testing...{model} model, {read_data} data')
    import timeit
    if model == 'kilonova' or model == 'kilonova_uvboost':
        #uniform_random = np.random.rand(ndim)
        #theta = priortransform(uniform_random)
        #### parameters
        if model == 'kilonova':
            theta = np.array((0.05,0.1,0.2,0.4,3.0,0.5,4.5))
        elif model == 'kilonova_uvboost':
            theta = np.array((0.05,0.1,0.2,0.23,3.0,0.04,4.5))
        logl = loglikelihood(theta)
        print('theta: ',theta)
        print('log(likelihood): ', logl)
        timing_loglikelihood = int(np.round(1e3 * np.median(timeit.repeat(
            'loglikelihood(theta)',
            globals=globals(), number=1, repeat=10))))
        print(f'One loglikelihood calculation = {timing_loglikelihood} ms')
        theta_reshaped = np.array((theta[0] * u.Msun, np.array((theta[1], theta[2], theta[3])) * c.c, np.array((theta[4], theta[5])) * u.cm**2 / u.g, theta[6]), dtype= 'object')
        abmags_model = lightcurve_model(ts_data, theta_reshaped, bs)
    elif model == 'shock':
        k = 10 # 0.1 cm^2/g
        M = 0.5 #0.01 solar masses
        v = 2 #0.1c
        R = 10 #10^10 cm #Initial radius for shock
        theta = k, M, v, R
        # uniform_random = np.random.rand(ndim)
        # theta = priortransform(uniform_random)
        logl = loglikelihood(theta)
        print('theta: ',theta)
        print('log(likelihood): ', logl)
        timing_loglikelihood = int(np.round(1e3 * np.median(timeit.repeat(
            'loglikelihood(theta)',
            globals=globals(), number=1, repeat=10))))
        print(f'One loglikelihood calculation = {timing_loglikelihood} ms')
        abmags_model = lightcurve_model(ts_data, theta, bs)
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots()
    for key in bs:
        ax.plot(ts_data[key].to_value('day'),abmags_data[key],'x')
        ax.plot(ts_data[key].to_value('day'),abmags_model[key])
    print_string = f'output_files/plots/test_{model}model_{read_data}data_{delay}h_delay_{dist}Mpc_{optical_string}band_{uv_string}band.png'    
    fig.savefig(print_string)
    print(f'saved in {print_string}')
elif method == 'timeout' or method == 'pool':
    ########## NESTED SAMPLER #########
    start_time = time.time()
    folderstring = f'output_files/results/{model}model_{read_data}data_{delay}h_delay'
    try:
        os.mkdir(folderstring)
        print(f'Created directory: {folderstring}')
    except:
        print(f'Directory already exists: {folderstring}')

    filestring = f'{dist}Mpc_{optical_string}band_{uv_string}band'
    if not os.path.exists(folderstring+f'/{filestring}_results'):
        with open(folderstring+f'/{filestring}_priorlims','wb') as kilonova_limits :
            pickle.dump(limits, kilonova_limits)
    
        if method == 'pool':
            from schwimmbad import MultiPool
            print('initiating sampler...')
            with MultiPool() as pool:
                sampler = NestedSampler(loglikelihood, priortransform, ndim, pool=pool)
                print('poolsize = ',pool.size)
                sampler.run_nested(print_progress=print_progress)
        
        elif method == 'timeout':
            from Timeout import timeout
            print('initiating sampler...')
            sampler_start = time.time()
            sampler = NestedSampler(loglikelihood, priortransform, ndim,sample=sample_method)
            sampler_time = int(np.ceil(time.time()-sampler_start))
            print(f'sampler initiated in {sampler_time} seconds')
            try:
                with timeout(seconds=max_time-sampler_time):
                    sampler.run_nested(print_progress=print_progress)
            except TimeoutError:
                pass
        else:
            print('error in method')
    
        print("--- %s seconds ---" % (time.time() - start_time))
        with open(folderstring + f'/{filestring}_results', 'wb') as kilonova_results :
            pickle.dump(sampler.results, kilonova_results)
    else:
        print(f'{filestring}_results already exists, skipping...')
