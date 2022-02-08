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

class Lightcurve():
    def __init__(self,distance, bandpass = dorado.sensitivity.bandpasses.NUV_D,
                 heating_function = 'korobkin'):
        self.distance = distance # u.Mpc
        self.heating_function = heating_function
        self.bandpass = bandpass

    def get_sed_kilonova(self,t_rel,theta):
        L, T, r = lightcurve(t_rel, theta[0], theta[1], theta[2], theta[3], heating_function = self.heating_function)
        seds = [synphot.SourceSpectrum(synphot.BlackBody1D,temperature=TT).taper(self.bandpass.waveset) * np.pi *(rr / self.distance).to(u.dimensionless_unscaled)**2 for TT, rr in zip(T, r)]
        return seds

    def calc_abmags_kilonova(self,t_rel,theta):
         L, T, r = lightcurve(t_rel, theta[0], theta[1], theta[2], theta[3], heating_function = self.heating_function)
         abmags = get_abmag(T, r, self.distance, self.bandpass)
         return abmags
     
    def calc_snrs_kilonova(self,t_rel,theta,ts_time_object,coords,exptimes = 10*u.min, nighttimes=True):
        seds = self.get_sed_kilonova(t_rel,theta)
        SNRs = [dorado.sensitivity.get_snr(sed, exptime=exptimes,
                                          coord=coord, night=nighttimes, time =
                                          t_time_object) for sed,
               coord, t_time_object in zip(seds,coords,
                                           ts_time_object)] 
        return np.array(SNRs)

    
    def get_sed_shock(self,t_rel,theta):
        k,M,V,R_0=theta
        T_sh = Teff_f(t_rel.to_value(u.day),k,M,V,R_0)*u.K
        R_sh = R_ph_f(t_rel.to_value(u.day),k,M,V,R_0)*u.cm
        seds = [synphot.SourceSpectrum(synphot.BlackBody1D,temperature=TT).taper(self.bandpass.waveset) * np.pi *(rr / self.distance).to(u.dimensionless_unscaled)**2 for TT, rr in zip(T_sh, R_sh)]
        return seds
    
    def calc_abmags_shock(self,t_rel,theta):
        '''Calculate AB mags for shock interaction model.

        The input parameters must be these two: t, theta, because this function
        will be used for the parameter estimation.
        '''
        k,M,V,R_0=theta
        T_sh = Teff_f(t_rel.to_value(u.day),k,M,V,R_0)*u.K
        R_sh = R_ph_f(t_rel.to_value(u.day),k,M,V,R_0)*u.cm
        abmags = get_abmag(T_sh, R_sh, self.distance, self.bandpass)
        return abmags
    
    def calc_snrs_shock(self,t_rel,theta,ts_time_object,coords,exptimes = 
                        10*u.min, nighttimes=True):    
        seds = self.get_sed_shock(t_rel,theta)
        SNRs = [dorado.sensitivity.get_snr(sed, exptime=exptimes,
                                          coord=coord, night=nighttimes, time =
                                          t_time_object) for sed, coord, 
                t_time_object in zip(seds,coords, ts_time_object)] 
        return np.array(SNRs)
