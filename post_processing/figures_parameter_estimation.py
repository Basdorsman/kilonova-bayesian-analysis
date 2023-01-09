#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 14:21:31 2022

@author: bas
"""

from cornerplots import cornerPlot, getSamples
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

redden = True
optical_delays=12

kilonova_fontsize=22
shock_fontsize = 20
labelpad=0.08

#%% Figure 3

model = 'kilonova' # kilonova, kilonova_uvboost, shock
datas = 'kilonova'
delay = '0h'
distance = ['160','100','40']
band = 'no_opticalband_NUV_Dband'

legend_texts= ['160 Mpc','100 Mpc','40 Mpc']
colors=['blue','orange','green']
linestyles = ['solid','dashdot','dashed']
linewidths = [3,3,3]

samples,legend_texts = getSamples(model, datas, delay, distance, band, redden, optical_delays, legend_texts=legend_texts)
figure = cornerPlot(model, samples, legend_texts, colors, linestyles, linewidths, fill_contours=True)

plt.show()
figure.savefig(f'plots/cornerplot_{band}_{model}.png',dpi=300,pad_inches=0.3,bbox_inches='tight')

#%% Figure 4

model = 'shock' # kilonova, kilonova_uvboost, shock
datas = 'shock'
delay = '0h'
distance = ['160','100','40']
band = 'no_opticalband_NUV_Dband'
redden=True
legend_texts= ['160 Mpc','100 Mpc','40 Mpc']
colors=['blue','orange','green']
linestyles = ['solid','dashdot','dashed']
linewidths = [3,3,3]

samples,legend_texts = getSamples(model, datas, delay, distance, band, redden, optical_delays, legend_texts=legend_texts)
figure = cornerPlot(model, samples, legend_texts, colors, linestyles, linewidths)

plt.show()
figure.savefig(f'plots/cornerplot_{band}_{model}.png',dpi=300, pad_inches=0.3, bbox_inches='tight')

#%% Figure 5

model = 'kilonova' # kilonova, kilonova_uvboost, shock
datas = 'kilonova'
delay = '0h'
distance = '40'
band = ['rband_no_uvband','no_opticalband_NUV_Dband','rband_NUV_Dband']
redden=True
legend_texts=['r-band', 'UV-band', 'r-band and UV-band']
colors=['blue','orange','green']
linestyles = ['solid','dashdot','dashed']
linewidths = [3,3,3]

samples,legend_texts = getSamples(model, datas, delay, distance, band, redden, optical_delays, legend_texts=legend_texts)
figure = cornerPlot(model, samples, legend_texts, colors, linestyles, linewidths)

plt.show()
figure.savefig(f'plots/cornerplot_{band}_{model}.png',dpi=300,pad_inches=0.3,bbox_inches='tight')

#%% Figure 7

model = 'kilonova' # kilonova, kilonova_uvboost, shock
datas = 'kilonova'
delay = ['24h','12h','0h']
distance = '40' #['160','100','40']
#band = ['no_opticalband_NUV_Dband','no_opticalband_NUV_DD2band','rband_NUV_Dband']
band = 'no_opticalband_NUV_Dband'
redden=True
legend_texts=['Data starting at 25.2h','Data starting at 13.2h','Data starting at 1.2h']
colors=['blue','orange','green']
linestyles = ['solid','dashdot','dashed']
linewidths = [3,3,3]

samples,legend_texts = getSamples(model, datas, delay, distance, band, redden, optical_delays, legend_texts=legend_texts)
figure = cornerPlot(model, samples, legend_texts, colors, linestyles, linewidths)
plt.show()
figure.savefig(f'plots/cornerplot_{delay}_{band}_{model}.png',dpi=300,pad_inches=0.3,bbox_inches='tight')