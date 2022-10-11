#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 10:33:41 2022

@author: bas
"""

import numpy as np
np.random.seed(0)
from astropy import units as u
from astropy.coordinates import GeocentricTrueEcliptic, get_sun, SkyCoord
from astropy.time import Time
import synphot
import dorado.sensitivity
from synphot.units import PHOTLAM
import matplotlib.pyplot as plt
from dorado.sensitivity import _get_reddening_law
from dorado.sensitivity import _get_dust_query
reddening_law = _get_reddening_law()
dust_query = _get_dust_query()


time = Time('2020-10-31 12:33:12')
sun = get_sun(time).transform_to(GeocentricTrueEcliptic(equinox=time))

source_spectrum = synphot.SourceSpectrum(synphot.ConstFlux1D, amplitude=21 * u.ABmag)
bandpass = dorado.sensitivity.bandpasses.NUV_D
coord = SkyCoord(sun.lon + 180*u.deg, 0*u.deg, frame=GeocentricTrueEcliptic(equinox=time))

def get_reddened_mag(coord, time, source_spectrum, bandpass):    
    ebv = dust_query(coord)
    extinction_curve = reddening_law.extinction_curve(ebv, bandpass.waveset)
    reddened_curve = source_spectrum * extinction_curve
    magnitude = synphot.Observation(reddened_curve, bandpass).effstim(u.ABmag)
    return magnitude

magnitude = get_reddened_mag(coord, time, source_spectrum, bandpass)
print(magnitude)

#%%

lons = 180+180*np.random.rand(1000)
lats = -90+180*np.random.rand(1000)

magnitude = []
for lon, lat in zip(lons, lats):
    coord = SkyCoord(sun.lon + lon*u.deg, lat*u.deg, frame=GeocentricTrueEcliptic(equinox=time))
    magnitude.append(get_reddened_mag(coord, time, source_spectrum, bandpass).value)


fig, ax = plt.subplots(1,1)

ax.plot(lons,magnitude)
ax.hist(magnitude,bins=100)
ax.set_xlim([21,24])
ax.set_ylabel('frequency')
ax.set_xlabel('reddened AB mag (input = 21 AB mag) ')
#ax.set_xticks([i+21 for i in range(20)])

#%%
import astropy.constants as const

def get_extinctioncurve(coord, time, bandpass):
    from dorado.sensitivity import _get_reddening_law
    from dorado.sensitivity import _get_dust_query
    reddening_law = _get_reddening_law()
    dust_query = _get_dust_query()
    ebv = dust_query(coord)
    extinction_curve = reddening_law.extinction_curve(ebv, bandpass.waveset)
    return extinction_curve

extinction_curve = get_extinctioncurve(coord, time, bandpass)

def get_flux(T, r, distance, bandpass, z=0, extinction=False):
    h = const.h.cgs
    c = const.c.cgs
    k_B = const.k_B.cgs
    
    wav = bandpass.waveset/(1+z)
    # wav_nounits = wav.value
    
    # Calculate spectral radiance : np.array(len(T),len(wav))
    B_l = 2*h*c**2/np.multiply(wav**5,np.exp(h*c/(k_B*np.outer(T,wav)))-1)
    
    # Scale spectral density = flux
    f_l = B_l.T*np.pi*(r / distance)**2/(1+z)
    
    # redden
    if not isinstance(extinction, bool):
        f_l = f_l.T * extinction._model.lookup_table
        f_l = f_l.T
    
    # Convert flux to photlam and multiply by bandpass
    f_PHOTLAM = (f_l.T/((h*c/bandpass.waveset))).to(1/u.cm**2/u.s/u.Angstrom)*bandpass(bandpass.waveset)
    return f_PHOTLAM.T

T = 10000 * u.Kelvin
distance = 40 * u.Mpc
r = 10 * u.pc

flux_source = get_flux(T, r, distance, bandpass, 0)
flux_extinction = get_flux(T, r, distance, bandpass, 0, extinction=extinction_curve)
z = 0.03
flux_reddened = get_flux( T, r, distance, bandpass, z)
flux_both = get_flux( T, r, distance, bandpass, z, extinction=extinction_curve)

fig, ax = plt.subplots(1,1)
ax.plot(bandpass.waveset,flux_source,label='neither')
ax.plot(bandpass.waveset,flux_extinction,label='MW extinction')
ax.plot(bandpass.waveset,flux_reddened,label='40 mpc redshifted')
ax.plot(bandpass.waveset,flux_both,label='both')

ax.set_xlim([1500,3000])
ax.legend()

#%%

def get_abmag(T, r, distance, bandpass, z=0, extinction=False):
    '''calculate AB magnitude faster.
    
    This code is largely copied from Synphot, using: synphot.BlackBody1D,
    synphot.SourceSpectrum and synphot.Observation. However, it's vectorized. 
    For multiple bandpasses, I use loops.

    Parameters
    ----------
    T : 1darray * quantity temperature
        Blackbody Temperature
    r : 1darray * quantity length
        Blackbody radius
    distance : quantity length
        Luminosity distance
    bandpass : object synphot.spectrum.SpectralElement
    z : quantity redshift
    extinction : boolean (False) or object synphot.reddening.ExtinctionCurve
    
    Returns
    -------
    abmag : list of 1darrays 
        list of AB magnitudes. One entry per bandpass
    '''

    h = const.h.cgs
    c = const.c.cgs
    k_B = const.k_B.cgs
    
    wav = bandpass.waveset/(1+z)
    #wav_nounits = wav.value
    
    # Calculate spectral radiance : np.array(len(T),len(wav))
    B_l = 2*h*c**2/np.multiply(wav**5,np.exp(h*c/(k_B*np.outer(T,wav)))-1)
    
    # Scale spectral density = flux
    f_l = B_l.T*np.pi*(r / distance)**2/(1+z)
    
    # redden
    if not isinstance(extinction, bool):
        f_l = f_l.T * extinction._model.lookup_table
        f_l = f_l.T
    
    # Convert flux to photlam and multiply by bandpass
    f_PHOTLAM = (f_l.T/((h*c/bandpass.waveset))).to(1/u.cm**2/u.s/u.Angstrom)*bandpass(bandpass.waveset)
    num = abs(np.trapz(f_PHOTLAM.value * bandpass.waveset.value, x=bandpass.waveset.value)) 
    den = np.ones(len(T))*abs(np.trapz(bandpass(bandpass.waveset) * bandpass.waveset.value,
                                       x=bandpass.waveset.value))

    # Convolve with bandpass for total photon flux
    val = (num / den)*synphot.units.PHOTLAM

    # Convert flux from PHOTLAM to AB magnitude.
    abmag = synphot.units.convert_flux(bandpass.pivot(),val,u.ABmag).value
    return abmag


T = [6000 * u.Kelvin, 4000 * u.Kelvin]
distance = 40 * u.Mpc
r = 0.003 * u.pc

ABmag = get_abmag(T, r, distance, bandpass, z=0, extinction = False)