# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 18:13:31 2021

@author: super
"""

import numpy as np
import dorado.sensitivity
import synphot
from kilonova_heating_rate import lightcurve
import astropy.units as u
from shock_model import Teff_f,R_ph_f
from bol_to_band import get_abmag
from ETC import ETC

class Lightcurve():
    def __init__(self,distance, heating_function = 'beta'):
        self.distance = distance # u.Mpc
        self.heating_function = heating_function

    def get_blackbody(self, t_rel, theta, radiation = 'kilonova'):
        if radiation == 'kilonova' or radiation == 'kilonova_uvboost':
            L, T, r = lightcurve(t_rel, theta[0], theta[1], theta[2], theta[3], heating_function = self.heating_function)
        elif radiation == 'shock':
            k,M,V,R_0=theta
            T = Teff_f(t_rel.to_value(u.day),k,M,V,R_0)*u.K
            r = R_ph_f(t_rel.to_value(u.day),k,M,V,R_0)*u.cm
        return T, r

    def get_sed(self,t_rel,theta, bandpass, radiation = 'kilonova'):
        T, r = self.get_blackbody(t_rel, theta, radiation = radiation)
        seds = [synphot.SourceSpectrum(synphot.BlackBody1D,temperature=TT).taper(bandpass.waveset) * np.pi *(rr / self.distance).to(u.dimensionless_unscaled)**2 for TT, rr in zip(T, r)]
        return seds

    def calc_abmags(self,t_rel,theta, bandpasses, bandpasses_name, radiation = 'kilonova'):
        T, r = self.get_blackbody(t_rel, theta, radiation = radiation)
        if isinstance(bandpasses, dict):
            abmags = {}
            for key in zip(bandpasses):
                abmags[key] = get_abmag(T, r, self.distance, bandpasses[key])  #list of bandpasses
        else:
            abmags = get_abmag(T, r, self.distance, bandpasses) #single bandpass
        return abmags
     
    def insort(self, t_dict, kind='mergesort'):
        '''Create a single data series from a dictionary of various
        data series. Strips duplicate data. Mergesort seemed a tiny bit
        faster for my sorted large array, so that's my default.
        
        Parameters
        ----------
        t_dict : dictionary of data series
        
        Returns
        -------
        combination : stripped array of data points.
        
        '''
        
        arrays = [t_dict[key].to_value('day') for key in t_dict]
        combination = np.array([])
        for array in arrays:
            combination = np.concatenate([combination, array])
        combination.sort(kind=kind)
        flag = np.ones(len(combination), dtype=bool)
        np.not_equal(combination[1:], combination[:-1], out=flag[1:])
        return combination[flag]
    
    def is_present(self,my_list,compare_list):
        '''For each item in my_list, check if it is present in compare_list.
        
        Parameters
        ----------
        my_list : array
        compare_list : array
        
        Returns
        -------
        output : list, same shape as my_list
        '''

        output = []
        for my_item in my_list:
                is_present = np.any(compare_list==my_item*np.ones(len(compare_list)))
                output.append(is_present)
        return output

    def calc_abmags_combined(self,t_dict,theta,bandpasses_dict, radiation = 'kilonova'):
        '''Efficiently calculate abmags for multiple sets of time series.
        
        By leveraging the insort function, the expensive lightcurve function
        is called only a mininimum amount of times. This is efficient if t_dict
        contains duplicate time data. 
        
        Parameters
        ----------
        t_dict : dictionary of time series
        theta : light curve parameters
        bandpasses_dict : dictionary of bandpasses
        
        Returns
        -------
        abmags : array of abmags
        '''

        t_combined = self.insort(t_dict)
        T, r = self.get_blackbody(t_combined*u.day, theta, radiation = radiation)
        abmags = {}
        for key in bandpasses_dict:
            flag = self.is_present(t_combined,t_dict[key].to_value('day'))
            abmags[key] = get_abmag(T[flag], r[flag], self.distance, bandpasses_dict[key])
        return abmags

    def calc_snrs_dorado(self,t_rel,theta,ts_time_object,coords,bandpass=dorado.sensitivity.bandpasses.NUV_D,radiation = 'kilonova'): #dorado only
        seds = self.get_sed(t_rel,theta,bandpass,radiation = radiation)
        SNRs=[dorado.sensitivity.get_snr(sed, exptime=10*u.min,
                                          coord=coord, night=True, time =
                                          t_time_object) for sed,
               coord, t_time_object in zip(seds,coords,
                                           ts_time_object)]
        return np.array(SNRs)

    def calc_snrs_optical(self, t_rel, theta, bandpass, bandpass_name, radiation = 'kilonova'):
        abmags = self.calc_abmags(t_rel, theta, bandpass, bandpass_name, radiation = radiation)
        SNRs = [ETC(mag = abmag,etime=180, filter_string=bandpass_name)['S_Nt'] for abmag in abmags]
        return np.array(SNRs)
