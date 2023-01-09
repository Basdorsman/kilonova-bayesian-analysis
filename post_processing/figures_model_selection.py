#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 14:10:40 2022

@author: bas
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from heatmaps import heatmap, annotate_heatmap, getLogZ, bayesPlot

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
cmaps= ("RdYlBu","Greens")
yticklabelvisibles = (True, False)

analysisModels = ['shock','kilonova_uvboost','kilonova']
analysisLabels=['Shock','Kilonova Lower\nEarly Opacity','Kilonova Default']

dataModels = ['shock','kilonova_uvboost','kilonova']
dataLabels=['Data: Shock','Data: Kilonova\nLower Early Opacity','Data: Kilonova']

cbarlabels=["Evidence: Log$_{10}$($\mathcal{Z}$)","Bayes' factor: Log$_{10}$($\mathcal{B}$)"]
redden=True

fontsize_fields = 14
fontsize_bayesplot=13

threshold_textcolor = (0.75, 0.75)


#%% Figure 2

optical_bands=['no_optical','r','r']
uv_bands=['NUV_D','no_uv','NUV_D']
bands = ['UV band', 'r band', 'UV and r band']
markers = ['o', 'P', '*']
bayesColors=['tab:orange','tab:blue','tab:green']
dists=[40, 70, 100, 130, 160]
kilonovamodels = ['kilonova','kilonova_uvboost']
bayesTitles = ['Default Kilonova Model', 'Lower Early Opacity Kilonova Model']
legends = [False, True]
ylabels = [True, False]

fig, axes = plt.subplots(1,2, sharey=True, figsize=(11.5,3.5))
for ax, kilonovamodel, title, legend, ylabel in zip(axes, kilonovamodels, bayesTitles, legends, ylabels):
    for uv_band, optical_band, band, marker, color in zip(uv_bands, optical_bands, bands, markers, bayesColors):
        fig, ax = bayesPlot(fig, ax, models=['shock', kilonovamodel], datas=['shock',kilonovamodel], dists=dists, optical_band=optical_band, uv_band=uv_band, redden=redden, legend_labels=[f'Shock data, {band}',f'Kilonova data, {band}'], marker=marker, color=color)
    ax.set_title(f'{title}',fontsize=fontsize_bayesplot-1)
    ax.set_yscale('symlog', linthresh=2)
    ax.hlines(2,40,160,color='k',label='$\log_{10}(\mathcal{B})=2$')
    ax.tick_params(axis='x', labelsize=fontsize_bayesplot-1)
    ax.set_xlabel('$d_L$ (Mpc)',fontsize=fontsize_bayesplot-1)
    if legend:
        ax.legend(bbox_to_anchor=(1,1),loc='upper left',fontsize=fontsize_bayesplot-2)
    if ylabel:
        ax.set_ylabel('$\log_{10}(\mathcal{B})$',fontsize=fontsize_bayesplot-1)
        ax.tick_params(axis='y', labelsize=fontsize_bayesplot-1)
fig.tight_layout()
plt.show()
fig.savefig('./plots/bayes_plot_bands.png',dpi=300)

#%% Figure 6

optical_band='no_optical'
uv_band='NUV_D'
delays = [0, 2, 4, 12]
first_datas=['1.2','3.2','5.2','13.2']
markers = ['o', 'P', '*', 'D', 'v']
bayesColors=['tab:orange','tab:blue','tab:green', 'tab:red', 'tab:purple']
dists=[40, 70, 100, 130, 160]
kilonovamodels = ['kilonova','kilonova_uvboost']
bayesTitles = ['Default Kilonova Model', 'Lower Early Opacity Kilonova Model']
legends = [False, True]
ylabels = [True, False]

fig, axes = plt.subplots(1,2, sharey=True, figsize=(11.5,3.5))
for ax, kilonovamodel, title, legend, ylabel in zip(axes, kilonovamodels, bayesTitles, legends, ylabels):
    for delay, first_data, marker, color in zip(delays, first_datas, markers, bayesColors):
        fig, ax = bayesPlot(fig, ax, models=['shock', kilonovamodel], datas=['shock',kilonovamodel], delay=delay, dists=dists, optical_band=optical_band, uv_band=uv_band, redden=redden, legend_labels=[f'Shock data, {first_data}h',f'Kilonova data, {first_data}h'], marker=marker, color=color)
    ax.set_title(f'{title}',fontsize=fontsize_bayesplot-1)
    ax.set_yscale('symlog', linthresh=2)
    ax.hlines(2,40,160,color='k',label='$\log_{10}(\mathcal{B})=2$')
    ax.tick_params(axis='x', labelsize=fontsize_bayesplot-1)
    ax.set_xlabel('$d_L$ (Mpc)',fontsize=fontsize_bayesplot-1)
    if legend:
        ax.legend(bbox_to_anchor=(1,1),loc='upper left',fontsize=fontsize_bayesplot-2)
    if ylabel:
        ax.set_ylabel('$\log_{10}(\mathcal{B})$',fontsize=fontsize_bayesplot-1)
        ax.tick_params(axis='y', labelsize=fontsize_bayesplot-1)
fig.tight_layout()
plt.show()
fig.savefig('./plots/bayes_plot_delays.png',dpi=300)