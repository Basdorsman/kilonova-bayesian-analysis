#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 14:26:13 2022

@author: bas
"""
import dill as pickle
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
import astropy.units as u
import numpy as np
import astropy.constants as c
import synphot as sp
import dorado.sensitivity

import sys
sys.path.insert(0, '../')
from produce_lightcurve import Lightcurve

read_datas=['shock','kilonova']
colors=['tab:blue','tab:orange']
linestyles=['none','none']
bands=('r','NUV_D')
dist=40
bs_optical_string='ugri'
bs_uv_string='NUV_DD1D2'
delay=0
fontsize=11

ts_data = {}
abmags_data = {}
snrs = {}
abmags_error = {}

fig, ax = plt.subplots(figsize=(4.5,3))


for read_data in read_datas:
    with open(f'../input_files/data/SNR_fiducial_{read_data}_{dist}Mpc_opticalbands_{bs_optical_string}_uvbands_{bs_uv_string}_{delay}h_delay.pkl','rb') as tf:
        ts_data[read_data], abmags_data[read_data], snrs[read_data], abmags_error[read_data]=pickle.load(tf)
    # UV and optical
    # for band,linestyle in zip(bands,linestyles):
    #     ax.errorbar(ts_data[read_data][band].to_value('hour'), abmags_data[read_data][band],yerr=abmags_error[read_data][band], color='k', linestyle=linestyle)
    # UV
    # ax.errorbar(ts_data[read_data]['NUV_D'].to_value('hour'), abmags_data[read_data]['NUV_D'],yerr=abmags_error[read_data]['NUV_D'], color='k', linestyle='none')
    # optical
    ax.errorbar(ts_data[read_data]['r'].to_value('hour'), abmags_data[read_data]['r'],yerr=abmags_error[read_data]['r'], color='k', linestyle='none')

# t_uv = np.logspace(0, 1.4)*u.hour
t_r = np.logspace(1, 1.7)*u.hour
# ts = [t_uv, t_r]
distance = dist * u.Mpc
# b_NUVD = dorado.sensitivity.bandpasses.NUV_D
b_r = sp.SpectralElement.from_file('../input_files/bands/SLOAN_SDSS.r.dat')
# bs = [b_NUVD, b_r]
# b_names = ['NUV_D','r']
# lines = ['-','--']
legend_names = ['shock data','kilonova data']



lightcurve_object = Lightcurve(distance,heating_function='beta')

abmags={}
for radiation, color, legend_name in zip(read_datas,colors, legend_names):
    abmags[radiation]={}
    
    # #### parameters
    # if radiation == 'kilonova':
    #     mass = 0.05 * u.Msun
    #     velocities = np.asarray([0.1, 0.2, 0.4]) * c.c
    #     opacities = np.asarray([3.0, 0.5]) * u.cm**2 / u.g
    #     n = 4.5
    #     theta = np.array((mass, velocities, opacities, n),dtype=object)
    # elif radiation == 'shock':
    #     mass = 0.05 * u.Msun
    #     velocities = np.asarray([0.12, 0.31, 0.61]) * c.c
    #     opacities = np.asarray([1.8, 0.39]) * u.cm**2 / u.g
    #     n = 4.47
    #     theta = np.array((mass, velocities, opacities, n),dtype=object)
        
    # abmags[radiation]['r'] = lightcurve_object.calc_abmags(t_r,theta,b_r,'r',radiation='kilonova')
    # ax.plot(t_r.to_value('hour'),abmags[radiation]['r'], color=color,linewidth=1,zorder=1, label=legend_name, linestyle='-')
    
    #### parameters
    if radiation == 'kilonova':
        # k in 0.1 cm^2/g, M in 0.01 solar masses, v in 0.1c, R_0 in 10^10 cm
        theta = [3, 0.5, 2, 5]
    elif radiation == 'shock':
        # k in 0.1 cm^2/g, M in 0.01 solar masses, v in 0.1c, R_0 in 10^10 cm
        theta = [5, 1, 2, 5]
        
    abmags[radiation]['r'] = lightcurve_object.calc_abmags(t_r,theta,b_r,'r',radiation='shock')
    ax.plot(t_r.to_value('hour'),abmags[radiation]['r'], color=color,linewidth=1,zorder=1, label=legend_name, linestyle='-')
    


ax.vlines(15, ymin=17, ymax=18.5, color = 'grey', linestyle='-.', label='AT2017gfo')
ax.invert_yaxis()
ax.set_xlabel('t (hours)',fontsize=fontsize)
ax.set_ylabel('Apparent Magnitude (AB)',fontsize=fontsize)
ax.legend(fontsize=fontsize-2)
fig.tight_layout()