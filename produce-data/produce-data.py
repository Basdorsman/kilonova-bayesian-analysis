import numpy as np
import synphot as sp
import pandas as pd
from dorado.scheduling import Orbit, FOV
from astropy.time import Time
from astropy.coordinates import ICRS, SkyCoord
from importlib import resources
from astropy.table import QTable
from astropy_healpix import HEALPix
import dorado.sensitivity
import astropy.constants as c
import astropy.units as u

###### .py files in parent directory
import os
import sys
sys.path.insert(0, './')
from produce_lightcurve import Lightcurve
from bol_to_band import mag_AB_error
from parameters import getParameters

# get parameters
parameters = getParameters(osargs_list=['delay','dist','read_data','redden','optical_delay'])

dist = int(parameters['dist'])
read_data = parameters['read_data']
optical_delay = int(parameters['optical_delay']) # hours
delay = int(parameters['delay']) # hours
redden = parameters['redden']=='True'
# print('redden:',redden)

#### parameters
if read_data == 'kilonova':
    mass = 0.05 * u.Msun
    velocities = np.asarray([0.1, 0.2, 0.4]) * c.c
    opacities = np.asarray([3.0, 0.5]) * u.cm**2 / u.g
    n = 4.5
    theta = np.array((mass, velocities, opacities, n),dtype=object)
    radiation = 'kilonova'
elif read_data == 'kilonova_uvboost':
    mass = 0.05 * u.Msun
    velocities = np.asarray([0.1, 0.2, 0.23]) * c.c
    opacities = np.asarray([3.0, 0.04]) * u.cm**2 / u.g
    n = 4.5
    theta = np.array((mass, velocities, opacities, n),dtype=object)
    radiation = 'kilonova'
elif read_data == 'shock':
    # k in 0.1 cm^2/g, M in 0.01 solar masses, v in 0.1c, R_0 in 10^10 cm
    theta = [5, 1, 2, 5]
    radiation = 'shock'

distance = dist * u.Mpc

b_NUVD = dorado.sensitivity.bandpasses.NUV_D
b_D1 = dorado.sensitivity.bandpasses.D1
b_D2 = dorado.sensitivity.bandpasses.D2
bs_uv = [b_NUVD, b_D1, b_D2]
bs_uv_name = ['NUV_D','D1','D2']

#D1 and D2 break dust map
bs_uv = [b_NUVD]
bs_uv_name = ['NUV_D']


b_u = sp.SpectralElement.from_file('./input_files/bands/SLOAN_SDSS.u.dat')
b_g = sp.SpectralElement.from_file('./input_files/bands/SLOAN_SDSS.g.dat')
b_r = sp.SpectralElement.from_file('./input_files/bands/SLOAN_SDSS.r.dat')
b_i = sp.SpectralElement.from_file('./input_files/bands/SLOAN_SDSS.i.dat')
bs_optical = [b_u, b_g, b_r, b_i]
bs_optical_name = ['u', 'g', 'r', 'i']


bs_optical_string  = ''
for band in bs_optical_name:
    bs_optical_string += band

bs_uv_string  = ''
for band in bs_uv_name:
    bs_uv_string += band


# define optical observation time
t_start = optical_delay
t_end = int(36)+optical_delay
t_optical = np.linspace(t_start,t_end,int((t_end-t_start)/12+1))*u.hour



###### import event data, event number 10 is an example of successful detection
colnames = ['alpha','alpha1','alpha2','alpha3','alpha4','alpha5','alpha6','amp_order','bandpass','beta','coa_phase','distance','eff_dist_g','eff_dist_h','eff_dist_l','eff_dist_t','eff_dist_v','end_time_gmst','eta','f_final','f_lower','g_end_time','g_end_time_ns','geocent_end_time','geocent_end_time_ns','h_end_time','h_end_time_ns','inclination','l_end_time','l_end_time_ns','latitude','longitude','mass1','mass2','mchirp','numrel_dat','numrel_mode_max','numrel_mode_min','phi0','polarization','process_id','psi0','psi3','simulation_id','source','spin1x','spin1y','spin1z','spin2x','spin2y','spin2z','t_end_time','t_end_time_ns','taper','theta0','v_end_time','v_end_time_ns','waveform','filler']
events_data = pd.read_csv('./input_files/sim_inspiral table.txt',delimiter=",",names=colnames)
end_times = events_data['geocent_end_time']
latitudes = events_data['latitude']*u.rad
longitudes = events_data['longitude']*u.rad
distances = events_data['distance']
masses1 = events_data['mass1']
masses2 = events_data['mass2']

idx = 10

GW_detection_time = Time(end_times[idx],format='gps').utc
healpix = HEALPix(nside=32, frame=ICRS())
target = SkyCoord(longitudes[idx],latitudes[idx])
target_pixel = healpix.skycoord_to_healpix(target)
schedule = QTable.read('./input_files/{}.ecsv'.format(idx))
fov = FOV.from_rectangle(7.1 * u.deg)
footprints = [fov.footprint_healpix(healpix, row['center'], row['roll']) for row in schedule]
schedule['found'] = [target_pixel in footprint for footprint in footprints]

# Fetch schedule UV
with resources.path('dorado.scheduling.data','dorado-625km-sunsync.tle') as path:
    orbit = Orbit(path)
T_orbit = orbit.period

# calc exposure schedule UV
t_half_exposure = 5*u.min
t_schedule = schedule['time'][schedule['found']] 
t_schedule_concatenated = t_schedule+t_half_exposure
t_exposure = t_schedule-GW_detection_time+t_half_exposure
t_exposure_concatenated = t_exposure.value

coord = schedule['center']
coord_concatenated = coord
coord_found = schedule['center'][schedule['found']]
coord_found = coord_found[0]


orbits = 14
start_time = 0.02 #days

for i in range(orbits-1):
    delta_t = (i+1)*T_orbit
    t_exposure_to_append = t_exposure.value+delta_t.to_value(u.day)
    t_schedule_to_append = t_schedule+delta_t
    t_exposure_concatenated = np.concatenate([t_exposure_concatenated,t_exposure_to_append])
    t_schedule_concatenated = np.concatenate([t_schedule_concatenated,t_schedule_to_append])
    coord_concatenated = np.concatenate([coord_concatenated, coord])

# Removing very early times, I'm aware small time data sometimes introduces issues
t_UV_data = t_exposure_concatenated[t_exposure_concatenated>start_time]*u.day
len_to_remove = len(t_schedule_concatenated)-len(t_UV_data)
t_UV_object = t_schedule_concatenated[len_to_remove:]
coord_concatenated = coord_concatenated[len_to_remove:]
# Implement delay, on t_UV_data only, not the t_UV_object because I want to keep the satellite-object relations the same.
t_UV_data = np.asarray([element.value + delay/24 for element in t_UV_data])*u.day


# Producing AB mags and corresponding SNRs
lightcurve_object = Lightcurve(distance,heating_function='beta')

# Initialize empty dictionaries
t_data = {}
abmags = {}
snrs = {}  
AB_error = {}
extinction_curves = {}
# print('tests')
# print(len(t_UV_data))
# print(len(coord_concatenated))

# print("tests")

# print(coord_found[0])
# print(coord_concatenated[0])


# generate extinction curves
for b, b_name in zip(bs_uv+bs_optical, bs_uv_name+bs_optical_name):
    if redden:
        extinction_curves[b_name] = lightcurve_object.get_extinctioncurve(coord_found, b) #for a single detection per orbit
    else:
        extinction_curves[b_name] = False

# generate data
for b, b_name in zip(bs_uv, bs_uv_name):

    t_data[b_name] = t_UV_data
    abmags[b_name] = lightcurve_object.calc_abmags(t_UV_data,theta,b,b_name,radiation=radiation,extinction=extinction_curves[b_name])
    snrs[b_name] = lightcurve_object.calc_snrs_dorado(t_UV_data,theta,t_UV_object,coord_found,bandpass=b,radiation=radiation, extinction=extinction_curves[b_name]) #only one coord
    AB_error[b_name] = mag_AB_error(snrs[b_name])


for b, b_name in zip(bs_optical, bs_optical_name):
    t_data[b_name] = t_optical
    abmags[b_name] = lightcurve_object.calc_abmags(t_optical,theta,b, b_name,radiation=radiation, extinction=extinction_curves[b_name])
    snrs[b_name] = lightcurve_object.calc_snrs_optical(t_optical, theta, b, b_name, radiation=radiation, extinction=extinction_curves[b_name])
    AB_error[b_name] = mag_AB_error(snrs[b_name])

# print(abmags)

# SAVE DATA
print(f"saving data to ./input_files/data/SNR_fiducial_{read_data}_{dist}Mpc_opticalbands_{bs_optical_string}_uvbands_{bs_uv_string}_{delay}h_delay_redden_{redden}_optical_delay_{optical_delay}.pkl'")
import pickle
data = [t_data, abmags, snrs, AB_error, extinction_curves]
with open(f'./input_files/data/SNR_fiducial_{read_data}_{dist}Mpc_opticalbands_{bs_optical_string}_uvbands_{bs_uv_string}_{delay}h_delay_redden_{redden}_optical_delay_{optical_delay}.pkl', 'wb') as tf:
    pickle.dump(data,tf)

# PLOT DATA
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

for b_name in bs_uv_name:
    ax.errorbar(t_UV_data.to_value('day'), abmags[b_name],yerr=AB_error[b_name],label=b_name)
for b_name in bs_optical_name:
    ax.errorbar(t_optical.to_value('day'), abmags[b_name],yerr=AB_error[b_name],label=b_name)

ax.invert_yaxis()
ax.set_ylabel('Magnitude (AB)')
ax.set_xlabel('t (days)')
ax.legend()
fig.show()

# SAVE PLOTS
print_string = f'produce-data/plots/data_{read_data}_{dist}_{delay}h_delay_optical_delay_{optical_delay}.png'  
try:    
    os.mkdir('produce-data/plots')
    print('Created folder for plots')
except:
    print('No new folder')
fig.savefig(print_string)
print(f'saved plot in {print_string}')




