# -*- coding: utf-8 -*-
"""
Created on Fri May 14 16:04:57 2021

@author: super
"""

from getdist import plots, MCSamples
import pickle

with open("kilonovamodel_kilonovadata 21-05-11 0433", 'rb') as analysis_results: 
    results1 = pickle.load(analysis_results)

with open("kilonovamodel_kilonovadata 21-05-12 0050", 'rb') as analysis_results:
    results2 = pickle.load(analysis_results)


ndim = results1.samples.shape[1]
names = ["x%s"%i for i in range(ndim)]
labels =  ["x_%s"%i for i in range(ndim)]
samples1 = MCSamples(samples=results1.samples,names = names, labels = labels,label='kilonova 40 Mpc')
samples2 = MCSamples(samples=results2.samples,names = names, labels = labels,label='kilonova 100 Mpc')

# Triangle plot
g = plots.get_subplot_plotter()
g.triangle_plot([samples1,samples2], filled=False)