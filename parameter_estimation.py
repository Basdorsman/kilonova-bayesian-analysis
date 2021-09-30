#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 09:28:48 2021

@author: basdorsman
"""

from dynesty import NestedSampler
import numpy as np
from astropy.table import QTable
from produce_lightcurve import Lightcurve
import synphot as sp
import dorado.sensitivity
import astropy.units as u
from astropy import constants as c
import datetime
import pickle 
import time
import os

########### PARAMETERS, CHANGE HERE ##########
model = 'shock' #shock, kilonova, kilonova_uvboost
read_data = 'shock_optical' #kilonova_optical, shock_optical
dist = 40 #mpc
include_uv = 'uv' #uv, no_uv
include_optical = ['u', 'g','r', 'i'] # ['u', 'g','r', 'i'], False
print_progress=True
method = 'test' #test, timeout, pool
max_time = 100*60*60 # seconds, parameter for 'timeout' method

######## MORE PARAMETERS, DONT TOUCH ##########
if include_optical == False:
    bs_string = 'no_optical'
else:
    bs_string = ''.join(include_optical)
distance = dist * u.Mpc
heating = 'beta'
time_data = 8 # parameter for '(...)_time' data

if model == 'shock':
    radiation = 'shock' 
elif model == 'kilonova' or model == 'kilonova_uvboost':
    radiation = 'kilonova'

if include_uv == 'uv':
    b_dorado = dorado.sensitivity.bandpasses.NUV_D
    bs = {'dorado' : b_dorado}
elif include_uv == 'no_uv':
    bs = {}
if not include_optical == False:
    for key in include_optical: 
        bs[key] = sp.SpectralElement.from_file(f'input_files/bands/SLOAN_SDSS.{key}.dat')


#### READ DATA #########

with open(f'input_files/data/SNR_fiducial_{read_data}_{dist}Mpc_ugriband.pkl','rb') as tf:
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
    limits=np.array(([0.01,10],[0.01,5],[0.01,3],[0.1,10])) #lim[[lower,upper],..
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
        theta = np.array((0.05,0.1,0.2,0.4,3.0,0.5,4.5))
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
    fig.savefig('output_files/plots/test.png')

elif method == 'timeout' or method == 'pool':
    ########## NESTED SAMPLER #########
    start_time = time.time()
    folderstring = f'output_files/{model}model_{read_data}data_{dist}Mpc_{bs_string}band_{include_uv}'
    try:
        os.mkdir(folderstring)
        print(f'Created directory: {folderstring}')
    except:
        print(f'Directory already exists: {folderstring}')

    daytimestring = datetime.datetime.now().strftime("%y-%m-%d %H%M")
    with open(folderstring + f'/{daytimestring} priorlims','wb') as kilonova_limits :
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
        sampler = NestedSampler(loglikelihood, priortransform, ndim)
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
    with open(folderstring + f'/{daytimestring} results', 'wb') as kilonova_results :
        pickle.dump(sampler.results, kilonova_results)