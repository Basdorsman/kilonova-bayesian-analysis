#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 13:34:57 2022

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

read_datas=['shock','kilonova','kilonova_uvboost']
colors=['tab:blue','tab:orange','tab:green']
linestyles=['none','none']
bands=('r','NUV_D')
dist=40
bs_optical_string='ugri'
bs_uv_string='NUV_D'#'NUV_DD1D2'
delay=0
redden='True'
optical_delay=12
fontsize=11

ts_data = {}
abmags_data = {}
snrs = {}
abmags_error = {}
extinction_curves = {}

fig, ax = plt.subplots(figsize=(4.5,3))


for read_data in read_datas:
    with open(f'../input_files/data/SNR_fiducial_{read_data}_{dist}Mpc_opticalbands_{bs_optical_string}_uvbands_{bs_uv_string}_{delay}h_delay_redden_{redden}_optical_delay_{optical_delay}.pkl','rb') as tf:
        ts_data[read_data], abmags_data[read_data], snrs[read_data], abmags_error[read_data], extinction_curves[read_data]=pickle.load(tf)
    # for band,linestyle in zip(bands,linestyles):
    #     ax.errorbar(ts_data[read_data][band].to_value('hour'), abmags_data[read_data][band],yerr=abmags_error[read_data][band], color='k', linestyle=linestyle)
    ax.errorbar(ts_data[read_data]['NUV_D'].to_value('hour'), abmags_data[read_data]['NUV_D'],yerr=abmags_error[read_data]['NUV_D'], color='k', linestyle='none')


t_uv = np.logspace(0, 1.4)*u.hour
t_r = np.logspace(1, 1.7)*u.hour
ts = [t_uv, t_r]
distance = dist * u.Mpc
b_NUVD = dorado.sensitivity.bandpasses.NUV_D
b_r = sp.SpectralElement.from_file('../input_files/bands/SLOAN_SDSS.r.dat')
bs = [b_NUVD, b_r]
b_names = ['NUV_D','r']
lines = ['-','--']
legend_names = ['Shock Interaction','Default','Lower Early Opacity']

lightcurve_object = Lightcurve(distance,heating_function='beta')

abmags={}
for radiation, color, legend_name in zip(read_datas, colors, legend_names):
    abmags[radiation]={}
    
    #### parameters
    if radiation == 'kilonova':
        mass = 0.05 * u.Msun
        velocities = np.asarray([0.1, 0.2, 0.4]) * c.c
        opacities = np.asarray([3.0, 0.5]) * u.cm**2 / u.g
        n = 4.5
        theta = np.array((mass, velocities, opacities, n),dtype=object)
    elif radiation == 'kilonova_uvboost':
        mass = 0.05 * u.Msun
        velocities = np.asarray([0.1, 0.2, 0.23]) * c.c
        opacities = np.asarray([3.0, 0.04]) * u.cm**2 / u.g
        n = 4.5
        theta = np.array((mass, velocities, opacities, n),dtype=object)
    elif radiation == 'shock':
        # k in 0.1 cm^2/g, M in 0.01 solar masses, v in 0.1c, R_0 in 10^10 cm
        theta = [5, 1, 2, 5]

    # for b, b_name, t, line in zip(bs, b_names, ts, lines):
    #     abmags[radiation][b_name] = lightcurve_object.calc_abmags(t,theta,b,b_name,radiation=radiation)
    #     ax.plot(t.to_value('hour'),abmags[radiation][b_name], color=color,linewidth=1,zorder=1, label=f'{radiation} {b_name}', linestyle=line)
    abmags[radiation]['NUV_D'] = lightcurve_object.calc_abmags(t_uv,theta,b_NUVD,'NUV_D',radiation=radiation, extinction=extinction_curves[radiation]['NUV_D'])
    ax.plot(t_uv.to_value('hour'),abmags[radiation]['NUV_D'], color=color,linewidth=1,zorder=1, label=legend_name, linestyle='-')
    
ax.vlines(15, ymin=17, ymax=23, color = 'grey', linestyle='-.', label='AT2017gfo')
ax.invert_yaxis()
ax.set_xlabel('t (hours)',fontsize=fontsize)
ax.set_ylabel('Apparent Magnitude (AB)',fontsize=fontsize)
ax.legend(fontsize=fontsize-2)
fig.tight_layout()
fig.savefig('plots/comparison_threeway.png',dpi=300)