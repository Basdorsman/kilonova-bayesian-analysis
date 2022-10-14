import numpy as np
from produce_lightcurve import Lightcurve
import synphot as sp
import dorado.sensitivity
import astropy.units as u
from astropy import constants as c
import dill as pickle 
import os
from dynesty_sampler import getSampler, wrappedSampler, find
from parameters import getParameters

# get parameters
parameters = getParameters(osargs_list=['read_data','model','delay','dist','include_optical','include_uv','print_progress','method','resume_previous','sample','save_after_seconds','parallel','dlogz_threshold','redden'])
# parameters = {
#     'model' : 'kilonova_uvboost',
#     'delay' : 0,
#     'dist' : 40,
#     'include_optical' : 'False',
#     'print_progress' : 'True',
#     'method' : 'sample',
#     'include_uv' : 'NUV_D',
#     'read_data' : 'shock',
#     'sample' : 'auto',
#     'resume_previous' : 'False',
#     'save_after_seconds' : 600, #float/int or will resort to False
#     'parallel' : 'True',
#     'dlogz_threshold' : 10000
#     }

model = parameters['model'] #shock, kilonova, kilonova_uvboost
delay = parameters['delay'] #hours
dist = int(parameters['dist']) #mpc
include_optical = parameters['include_optical'].split(',') # 'r', ['u', 'g','r', 'I', 'z'], ['False']
print_progress=parameters['print_progress']=='True'
method = parameters['method'] #'plot', 'sample'
include_uv = parameters['include_uv'].split(',')
read_data = parameters['read_data']
sample = parameters['sample']
resume_previous = parameters['resume_previous']=='True'
redden = parameters['redden']=='True'
try:
    save_after_seconds = int(parameters['save_after_seconds'])
except:
    save_after_seconds = False
try:
    parallel = int(parameters['parallel'])
except:
    parallel = False
dlogz_threshold=float(parameters['dlogz_threshold'])

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
with open(f'input_files/data/SNR_fiducial_{read_data}_{dist}Mpc_opticalbands_ugri_uvbands_NUV_D_{delay}h_delay_redden_{redden}.pkl','rb') as tf:
    data_list = pickle.load(tf)
ts_data, abmags_data, snrs, abmags_error, extinction_curves = data_list



########## LOG PROBABILITIES ##########

if model == 'kilonova' or model == 'kilonova_uvboost':
    ndim = 7 
    limits = np.ndarray((ndim,2),dtype=object)
    limits[0] = [0.01, 0.1]  # mass (solar masses)
    if model == 'kilonova_uvboost':
        limits[1:4] = np.array(([0.05, 0.2], ['vmin', 'vmax'], [0.21, 0.8]), dtype=object)  # velocities (cm/s)
        limits[4:6] = np.array(([1,10],[0.01,0.1]))  # opacities (cm^2/g)
    elif model == 'kilonova':
        limits[1:4] = np.array(([0.05, 0.2], ['vmin', 'vmax'], [0.3, 0.8]), dtype=object)  # velocities (cm/s)
        limits[4:6] = np.array(([1,10],[0.1,1]))  # opacities (cm^2/g)
    limits[6] = [4,5]
elif model == 'shock':
    #limits=np.array(([0.01,10],[0.01,5],[0.01,3],[0.1,10])) #old broad
    limits=np.array(([1,10],[0.5,5],[1,3],[1,10])) #lim[[lower,upper],.. #arxiv limits
    #limits=np.array(([1,10],[0.5,5],[1,3],[4.9,5.1])) #lim[[lower,upper],.. #restrained limits
    # k in 0.1 cm^2/g, M in 0.01 solar masses, v in 0.1c, R_0 in 10^10 cm
    ndim = limits.shape[0]


########## DEFINE MODEL ##########
lightcurve_object = Lightcurve(distance, heating_function=heating)

def lightcurve_model(t,theta_reshaped,bandpasses,extinction = False):
    abmags = lightcurve_object.calc_abmags_combined(t,theta_reshaped,bandpasses,radiation = radiation, extinction = extinction)
    return abmags

folderstring = f'output_files/results/{model}model_{read_data}data_{delay}h_delay'
filestring = f'{dist}Mpc_{optical_string}band_{uv_string}band_redden_{redden}'

if not (resume_previous == True and isinstance(find(filestring+'_sampler_dlogz=*', folderstring), str)):
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
    
        def loglikelihood(theta):
            theta_reshaped = np.array((theta[0] * u.Msun, np.array((theta[1], theta[2], theta[3])) * c.c, np.array((theta[4], theta[5])) * u.cm**2 / u.g, theta[6]), dtype= 'object')
            sigmas2 = {}
            loglikelihood = {}
            abmags_model = lightcurve_model(ts_data, theta_reshaped, bs, extinction=extinction_curves)
            for key in bs:
                sigmas2[key] = abmags_error[key]**2
                loglikelihood[key] = -0.5 * np.sum((abmags_data[key] - abmags_model[key]) ** 2 / sigmas2[key] + np.log(2*np.pi*sigmas2[key]))
            return sum(loglikelihood.values())
        
    elif model == 'shock':
    
        def priortransform(u):
            theta = (limits[:,1]-limits[:,0])*u+limits[:,0]
            return theta
        
        def loglikelihood(theta):
            sigmas2 = {}
            loglikelihood = {}
            abmags_model = lightcurve_model(ts_data, theta, bs, extinction=extinction_curves)
            for key in bs:
                sigmas2[key] = abmags_error[key]**2
                loglikelihood[key] = -0.5 * np.sum((abmags_data[key] - abmags_model[key]) ** 2 / sigmas2[key] + np.log(2*np.pi*sigmas2[key]))
            return sum(loglikelihood.values())

####### TESTS / PARAMETER ESTIMATION #############
if method == 'plot':
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
        abmags_model = lightcurve_model(ts_data, theta_reshaped, bs, extinction=extinction_curves)
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
        abmags_model = lightcurve_model(ts_data, theta, bs, extinction=extinction_curves)
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots()
    for key in bs:
        ax.plot(ts_data[key].to_value('day'),abmags_data[key],'x')
        ax.plot(ts_data[key].to_value('day'),abmags_model[key])
    try:    
        os.mkdir('./output_files/plots')
        print('created folder for plots')
    except:
        print('folder for plots exists')
    print_string = f'./output_files/plots/test_{model}model_{read_data}data_{delay}h_delay_{dist}Mpc_{optical_string}band_{uv_string}band_redden_{redden}.png'    
    fig.savefig(print_string)
    print(f'saved in {print_string}')

    ########## NESTED SAMPLER #########
elif method == 'sample':
    try:
        os.mkdir(folderstring)
        print(f'Created directory: {folderstring}')
    except:
        print(f'Opened directory: {folderstring}')
    #if not os.path.exists(folderstring+f'/{filestring}_results'):
    if not isinstance(find(filestring+'_results_dlogz=*', folderstring), str):
        with open(folderstring+f'/{filestring}_priorlims','wb') as prior_limits :
            pickle.dump(limits, prior_limits)
        with open(folderstring+f'/{filestring}_parameters','wb') as input_parameters :
            pickle.dump(parameters,input_parameters)
        try:
            priortransform
        except NameError:
            sampler, previous_dlogz = getSampler(ndim, folderstring, filestring, parallel=parallel, sample=sample, resume_previous=resume_previous)
            priortransform=sampler.prior_transform.func
            loglikelihood=sampler.loglikelihood.func
        else:
            sampler, previous_dlogz = getSampler(ndim, folderstring, filestring, loglikelihood=loglikelihood, priortransform=priortransform, parallel=parallel, sample=sample, resume_previous=resume_previous)
        wrappedSampler(sampler, folderstring, filestring, previous_dlogz=previous_dlogz, sample=sample, save_after_seconds=save_after_seconds, print_progress=print_progress, parallel=parallel, dlogz_threshold=dlogz_threshold)
    else:
        print(f'{filestring}_results already exists, skipping...')
        
        
