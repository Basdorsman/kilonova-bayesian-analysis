from dynesty import NestedSampler
import numpy as np
from astropy.table import QTable
from produce_lightcurve import Lightcurve
import astropy.units as u
from astropy import constants as c
import datetime
import pickle 
import time
import os

########### PARAMETERS ##########
dist = 40
distance = dist * u.Mpc
heating = 'beta'
model = 'kilonova_uvboost'
read_data = 'shock'
time_data = 12
method = 'pool'
max_time = 100*60*60 # seconds, 'timeout' method



############ READ DATA ###################
with open(f'input_files/SNR_fiducial_shock_40Mpc_opticalbands_ugri_uvbands_NUV_DD1D2_0h_delay.pkl','rb') as tf:
    data_list = pickle.load(tf)
ts_data, abmags_data, snrs, abmags_error = data_list

t_data = ts_data['NUV_D']
AB_mag_data=abmags_data['NUV_D']
AB_error = abmags_error['NUV_D']



########## DEFINE MODEL ##########
lightcurve_object = Lightcurve(distance, heating_function=heating)

def kilonova_model(t,theta_reshaped):
    abmags = lightcurve_object.calc_abmags_kilonova(t,theta_reshaped)
    return abmags

########## LOG PROBABILITIES ##########
ndim = 7 
limits = np.ndarray((ndim,2),dtype=object)
limits[0] = [0.01, 0.1]  # mass (solar masses)

if read_data == 'kilonova_uvboost':
    limits[1:4] = np.array(([0.05, 0.2], ['vmin', 'vmax'], [0.21, 0.8]),dtype=object)  # velocities (cm/s)
    limits[4:6] = np.array(([1,10],[0.01,0.1]))  # opacities (cm^2/g)
else:
    limits[1:4] = np.array(([0.05, 0.2], ['vmin', 'vmax'], [0.3, 0.8]),dtype=object)  # velocities (cm/s)
    limits[4:6] = np.array(([1,10],[0.1,1]))  # opacities (cm^2/g)
limits[6] = [4,5]


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

# log_likelihood
def loglikelihood(theta):
    theta_reshaped = np.array((theta[0] * u.Msun, np.array((theta[1], theta[2], theta[3])) * c.c, np.array((theta[4], theta[5])) * u.cm**2 / u.g, theta[6]), dtype= 'object')
    AB_mag_model = kilonova_model(t_data,theta_reshaped)
    sigma2 = AB_error**2
    return -0.5 * np.sum((AB_mag_data - AB_mag_model) ** 2 / sigma2 + np.log(2*np.pi*sigma2))


if method == 'test':
    ########## QUICKTEST #############
    if model == 'kilonova':
        theta = np.array((0.05,0.1,0.2,0.4,3.0,0.5,4.5))
    elif model == 'kilonova_uvboost':
        theta = np.array((0.05,0.1,0.2,0.23,3.0,0.04,4.5))
    logl = loglikelihood(theta)
    print('theta: ',theta)
    print('log(likelihood): ', logl)
    import timeit
    timing_loglikelihood = int(np.round(1e3 * np.median(timeit.repeat(
        'loglikelihood(theta)',
        globals=globals(), number=1, repeat=10))))
    print(f'One loglikelihood calculation = {timing_loglikelihood} ms')

elif method == 'timeout' or method == 'pool':
    ########## NESTED SAMPLER #########
    start_time = time.time()
    try:
        os.mkdir(f'output_files/kilonovamodel_{read_data}data_{dist}Mpc')
        print('Created directory')
    except:
        print('Directory already exists')

    daytimestring = datetime.datetime.now().strftime("%y-%m-%d %H%M")
    with open(f'output_files/kilonovamodel_{read_data}data_{dist}Mpc/{daytimestring} priorlims {time_data}','wb') as kilonova_limits :
        pickle.dump(limits, kilonova_limits)

    if method == 'pool':
        from schwimmbad import MultiPool
        print('initiating sampler...')
        with MultiPool() as pool:
            sampler = NestedSampler(loglikelihood, priortransform, ndim, pool=pool, dlogz=10)
            print('poolsize = ',pool.size)
            sampler.run_nested(print_progress=True)
    
    elif method == 'timeout':
        from Timeout import timeout
        print('initiating sampler...')
        sampler_start = time.time()
        sampler = NestedSampler(loglikelihood, priortransform, ndim)
        sampler_time = int(np.ceil(time.time()-sampler_start))
        print(f'sampler initiated in {sampler_time} seconds')
        try:
            with timeout(seconds=max_time-sampler_time):
                sampler.run_nested(print_progress=True,dlogz=10)
        except TimeoutError:
            pass
    else:
        print('error in method')

    print("--- %s seconds ---" % (time.time() - start_time))
    with open(f'output_files/kilonovamodel_{read_data}data_{dist}Mpc/{daytimestring} results {time_data}', 'wb') as kilonova_results :
        pickle.dump(sampler.results, kilonova_results)
