#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 15:00:52 2021

@author: basdorsman
"""

import numpy as np


def ETC(S_N = None, mag = None, etime=None, itel=1, ires=0, filter_string ='r', moonp=1, air=3):

#    # SEVEN INPUT PARAMETERS (ires depends on itel)
#    S_N = parseFloat(input.S_N)
#    mag = parseFloat(input.mag)
#    etime = input.etime !== "" ? parseInt(input.etime) : 0
#    itel = parseInt(input.itel)
#    ires = itel == 8 ? 1 : 0  # if itel is 8 (2m0 FLOYDS), set ires to 1
#    filter = input.filter
#    moonp = parseInt(input.moonp) # new = 0 half = 1 full = 2
#    air= parseFloat(input.air) # airmass

    adiam  = 3.0   # Aperture diameter in arcseck
    ifilt = 2      # Index for filter array default (2) is 'V' filter
#    Nlambda        # Nphotons/A/cm^2/sec
#    Nobj           # Nphotons/sec from object
#    g_e
#    MAGinst        # Calculated instrumental mag (eg. in Sextractor)
#    DNinst         # Calculated instrumental Flux in DN
#    Npix           # Number of pixels in aperture
#    Nas            # Number of square arcsec in aperture
#    Nbkgd          # Number of photons/sec in background
#    NbDN           # Number of DN in background in etime
#    throughput     # total QE
    endloop = 0
    S_Nt = S_N     # S/N achieved in exp time
    S_Nlim = S_N
    magt = mag     # Mag achieved in exp time
    t = etime      # exposure time in sec
#    NeObj
#    NeBkgd
#    NeDark         # Nelect due to dark current
#    NeRon          # Nelect due to readout
#    PkDN           # estimate of peak DN in aperture (1/3 of stated aperture)
#    result         # Options are: s => calculate S/N, m => calc mag, e => calc exposure time
    resoln = 5
    dark = 0.0     # dark current (e- per pixel per second)
    ron = 0.0      # readout noise (e- per pix)
    gain = 0.0     # gain (e- per DN)

    spfrac = 0.4   # ADDITIONAL transmission through spectrometers
    if itel == 1: 
        saturated = 100000.0
    else:
        saturated = 40000.0        # Saturation counts. Only Sinistro (itel=1) has higher limit.

    # Names of filters
    Filt= ['U', 'B', 'V', 'R', 'I', 'u', 'g', 'r', 'i', 'Z', 'Y']

    # Collecting areas per ITEL: 0m4, 1m0, 2m0, 2m0, 0m8, Context, 1m0, 0m8
    Carea = [1200.0, 6260.0, 27000.0, 27000.0, 4380.0, 160.0, 6260.0, 4380.0]

    # Pixel scale per ITEL: 0m4 SBIG, 1m0 Sinistro, 2m0 Merope, 2m0 Spectral, 2m0 MuSCAT3, Context, 1m0 SBIG, 0m8 SBIG
    # Apixel = [0.57, 0.39, 0.278, 0.309, 0.43, 9.5, 0.47, 0.57] # original
    Apixel = [0.57, 0.389, 0.278, 0.300, 0.27, 9.5, 0.464, 0.57] # from website

    # Filter central wavelengths in microns
    Fcent = [0.350, 0.437, 0.549, 0.653, 0.789, 0.354, 0.476, 0.623, 0.760, 0.853, 0.975] # 20200921
    #Fcent = [0.359, 0.437, 0.549, 0.653, 0.789, 0.355, 0.476, 0.623, 0.760, 0.853, 0.975] # original
    #Fcent = [0.350, 0.436, 0.545, 0.641, 0.798, 0.354, 0.477, 0.622, 0.755, 0.870, 1.004] # from website
    #Fcent = [0.350, 0.434, 0.535, 0.629, 0.883, 0.349, 0.474, 0.622, 0.765, 0.869, 1.005] # from ConfigDB

    # Filter nominal (effective) bandwidths in microns
    Fband = [0.050, 0.107, 0.083, 0.137, 0.128, 0.057, 0.140, 0.135, 0.148, 0.113, 0.118] # 20200921
    #Fband = [0.068, 0.107, 0.083, 0.137, 0.128, 0.062, 0.140, 0.135, 0.148, 0.113, 0.118] # original
    #Fband = [0.050, 0.089, 0.084, 0.158, 0.154, 0.057, 0.150, 0.139, 0.129, 0.104, 0.112] # from website
    #Fband = [0.050, 0.095, 0.085, 0.126, 0.312, 0.067, 0.151, 0.134, 0.139, 0.103, 0.108] # from ConfigDB

    # Standard flux (in Jy) for 0 mag (Vega = UBVRI or AB = ugriz) for each filter
    Fjy = [1755, 4050, 3690, 3060, 2540, 3680, 3631, 3631, 3631, 3631, 3631] # original
    # Fjy = [1810, 4260, 3640, 3080, 2550, 3631, 3631, 3631, 3631, 3631, 3631] # traditional(?) standard fluxes


    # Hayes & Latham extinction coeffs for 2200m (mag per airmass)
    Ext= [0.54, 0.23, 0.12, 0.09, 0.04, 0.59, 0.14, 0.08, 0.06, 0.04, 0.03] # 20200921
    #Ext= [0.49, 0.23, 0.12, 0.09, 0.04, 0.54, 0.14, 0.08, 0.06, 0.04, 0.03] # updated by S. Valenti Aug 2015
    #Ext=[0.45, 0.21, 0.12, 0.07, 0.03, 0.48, 0.16, 0.09, 0.04, 0.03, 0.03] # H&L ext coeffs for FTN 3000m
    #Ext=[0.63, 0.32, 0.18, 0.13, 0.07, 0.70, 0.26, 0.15, 0.08, 0.06, 0.03] # H&L ext coeffs for FTS 1100m
    #Ext=[0.49, 0.24, 0.13, 0.09, 0.04, 0.54, 0.19, 0.11, 0.05, 0.04, 0.03] # original (for 2200m)

    # Theoretical 100% M1 for 2m, 27000 cm^2. Reference for these?
    TM1 = [25.33, 26.52, 25.89, 26.07, 25.57, 26.05, 26.60, 26.27, 26.16, 25.60, 25.45]

    # M1 zero point magnitudes for each telescope/filter  - ADJUST ZERO-POINTS here
    M1 = [
        [18.0, 20.3,  20.7,  21.2,  20.3,  18.7,  21.0,  20.7,  20.1,  18.9,  17.8], #0m4K (1mK-1.8) (SBIG)
        [19.7, 22.77, 23.10, 23.3,  22.5,  20.4,  23.42, 23.32, 22.95, 22.1,  20.0], #1mF (Sinistro)
        [22.4, 24.1,  24.6,  24.3,  23.5,  23.0,  24.5,  24.6,  24.5,  23.7,  21.6], #2mE (Merope)
        [21.3, 24.4,  24.6,  24.9,  24.1,  22.0,  24.5,  24.6,  24.4,  23.7,  21.6], #2mF (Spectral)
        [22.4, 24.1,  24.6,  24.3,  23.5,  23.0,  24.3,  24.2,  23.65,  22.5,  21.6], #2mE (MuSCAT3)
        [16.5, 18.1,  18.5,  19.0,  18.1,  17.2,  18.7,  18.5,  17.9,  16.7,  15.6], #Context (1mK-3.98) (SBIG)
        [20.5, 22.15, 22.50, 23.0,  22.1,  21.2,  22.75, 22.50, 21.90, 20.7,  19.6], #1mK (SBIG)
        [20.1, 21.7,  22.10, 22.6,  21.7,  20.8,  22.3,  22.1,  21.5,  20.3,  19.2]  #0m8K (1mK-0.39) (SBIG)
   ]

    # Estimated throughput of 0m4/0m8/1m telescope/SBIG
    TPUT = [
        [0.05, 0.08, 0.19, 0.20, 0.17, 0.05, 0.14, 0.14, 0.09, 0.05, 0.02], #0m4K (SBIG)
        [0.05, 0.14, 0.35, 0.35, 0.25, 0.05, 0.24, 0.28, 0.22, 0.17, 0.03], #1mF (Sinistro)
        [0.07, 0.11, 0.30, 0.20, 0.15, 0.06, 0.15, 0.22, 0.22, 0.17, 0.03], #2mE (Merope)
        [0.05, 0.14, 0.35, 0.35, 0.25, 0.05, 0.24, 0.28, 0.22, 0.17, 0.03]  #2mF (Spectral)
    ]

    # Sky brightness in mags/sq-as (Vega or AB) for each filter at new, half, full moon
    Msky = [
        [23.0, 22.5, 21.6, 20.6, 19.8, 23.5, 22.0, 21.1, 20.6, 20.2, 19.4], # new moon
        [20.0, 20.5, 20.3, 20.0, 18.8, 21.0, 20.3, 20.2, 19.7, 19.2, 18.0], # half
        [17.0, 17.8, 17.5, 17.4, 17.0, 18.0, 17.6, 17.5, 17.5, 16.8, 16.5]  # full
    ]
    

    if (ires == 1): # if FLOYDS is selected
        resoln = 500.0 # R for FLOYDS is between 400 and 700
        itel = 2  # assume FLOYDS is like 2m0 Merope
        spfrac = 0.8 # Estimated fraction, selected by Pickles, to match exposure times derived by Sand.
    elif (ires == 2): #added to accommodate NRES, but this can't yet be selected.
        resoln = 50000.0
        spfrac = 0.4
        if itel != 4:
            itel= 1
    

    if ((S_Nt is not None) and (mag is not None)): # if values for S/N and Mag are entered, then find Exp Time.
        t = 1.0  # set initial value for exp time
        result = 'e'  
    elif ((S_Nt is not None) and (etime is not None)): # if values for S/N and Exp Time are entered, then find Mag.
        magt = 30.0 # set initial value for magnitude
        result = 'm'
    else: #Last option: Find S/N
        result = 's'
    

    # Here are the values for itel: 0 -> 0m4 SBIG, 1 -> 1m0 Sinistro, 2 -> 2m0 Merope, 3 -> 2m0 Spectral,
    #  4 -> 0m8 Merope, 5 -> Context, 6 -> 1m0 SBIG, 7 -> 0m8 SBIG
    carea = Carea[itel] # retrieve collecting area for selected telescope
    APM = 2.5 * np.log10(carea / 27000.0)  # mag correction (scale factor) for selected telescope
    apixel = Apixel[itel] # retrieve pixel scale for selected telescope/instrument

    ifilt = Filt.index(filter_string) # assign the Filt-array index of the selected filter to "ifilt"
    if (ifilt == -1): # How does this ever happen?
        ifilt = 2
    if (filter == 'w'): # But the 'w' filter can't be selected at the interface!
        ifilt= 7
    
    fband = Fband[ifilt] # Retrieve the band width for the selected filter.
    if (filter == 'w'):
        fband=0.423   #gp+rp+ip band. When does this get selected?

    extco = Ext[ifilt] # assign value of extinction coeff appropriate for selected filter to "extco"
    acorr = 1000.0 * (air - 1.0) * extco # acorr is [(Mag at air)-(Mag at zenith)] x 1000 recall: air x extco = mag
    maga = ((1000.0 * magt) + acorr) / 1000.0 # selected mag modified for airmass

    Msky = Msky[moonp][ifilt] # retrieve appropriate sky brightness for selected moon fraction and filter

    if (itel==0) :       # 0m4 SBIG
        throughput = TPUT[0][ifilt] 
        dark = 0.02
        ron  = 14.0
        gain = 2.2 # where does this come from?
    elif(itel==1):  # 1m0 Sinistro
        throughput = TPUT[1][ifilt]
        dark = 0.002
        ron  = 12.0 # original value: 14.0
        gain = 1.0
        #saturated = 140000.0
    elif(itel==2):  # 2m0 Merope
        throughput = TPUT[2][ifilt]
        dark = 0.04
        ron  = 4.0
        gain = 2.0
        if (ires==1):
            dark=0.00007 # change dark current for FLOYDS
        
    elif(itel==3):  # 2m0 Spectral
        throughput = TPUT[3][ifilt]
        dark = 0.002 # estimate?
        ron  = 11.3 # original value: 12.0 (fs02 = 11.3, fs03 = 8.2)
        gain = 8.0 # (fs02 = 8.0, fs03 = 5.2)
    elif (itel==4):  # 2m0 MuSCAT3
        throughput = TPUT[3][ifilt]
        dark = 0.005
        ron  = 12.0
        gain = 1.85
    elif (itel==5):  # Context Camera SBIG
        throughput = TPUT[0][ifilt]
        dark = 0.02 # estimate?
        ron  = 14.0 # estimate?
        gain = 2.2 # estimate?
    elif (itel==6):  # 1m0 SBIG
        throughput = TPUT[0][ifilt]
        dark = 0.02
        ron  = 13.5 # original value: 14.0
        gain = 1.4 # original value: 2.2
    else:        # 0m8 SBIG
        throughput = TPUT[0][ifilt]
        dark = 0.02
        ron  = 14.0
        gain = 2.2
    

    throughput = np.power(10.0, 0.4*(M1[itel][ifilt]-TM1[ifilt]-APM) )  # Throughput is ratio of fluxes.

    if(resoln>10):
        throughput *= spfrac # This only gets invoked if FLOYDS (or NRES) is selected.

    while (endloop < 1):
        Nlambda = Fjy[ifilt]/Fcent[ifilt]/6.6256 # Nphotons/cm^2/sec/A. Fjy in Jy Fcent in Âµm. 6.6256 conversion factor.
        if (resoln<10):  # all cameras except spectrographs
            Nobj    = Nlambda*carea*fband*10000/np.power(10.0, 0.4*maga) # Nphotons/sec from object through atmosphere. The 10000 is A -> Âµm.
        else: 
            Nobj    = Nlambda*carea*5000/resoln/np.power(10.0, 0.4*maga)
        if (ires==2 and maga>15):
            maga=15.0  # Sensible limits for NRES spectrograph. Currently (2 Sep 2015) unavailable.
        #   if(S_Nlim>50) S_Nlim=50.0  # sensible limits for spectrograph
      
        # g_e, MAGinst, DNinst are never used?
        g_e = 2.5*np.log10(gain/t)     # g_e in mags. 
        MAGinst= maga-TM1[ifilt] + APM + g_e     # MAGinst = instrumental mag. Not used again? 
        DNinst = np.power(10.0, (-0.4*MAGinst) ) # DNinst = a flux ratio. Not used again?

        Nas  = 0.7854*adiam*adiam # Number of sq.arcsec in aperture, simply PI*[D/2]^2.
        Npix = Nas/apixel/apixel  # Number of pix in aperture.

        Nbkgd = Nobj*np.power(10.0, -0.4*(Msky-maga)) # Nphotons/sec per sq.arcsec from sky bkgd.
        NbDN  = Nbkgd*apixel*apixel/gain # (Nphotons/sec per sq.arcsec)*([arcsec/pix]^2)*(DN/elect) = DN/pix per sec from sky bkgd
        Nbkgd = Nbkgd*Nas # Nphotons/sec over all sq.arcsec in aperture.

        Nobj = Nobj*throughput # Nphotons/sec, weighted by throughput
        Nbkgd= Nbkgd*throughput # Nphotons/sec, weighted by throughput
        NbDN = NbDN*throughput # DN/pix per sec, weighted by throughput

        NeDark= Npix*dark # Nelect/sec due to dark curren
        NeRon = Npix*ron*ron # Nelect due to readout ([e-*pix]^2)

        NeObj  = Nobj*t # Nelect from source over exp time
        NeBkgd = Nbkgd*t # Nelect from bkgd over exp time
        NbDN   = NbDN*t # DN/pix from sky skgd
        NeDark = NeDark*t # Nelect from dark current
        PkDN   = 9*NeObj/Npix/gain # peak DN in aperture (DN/pix)

        S_Nt = NeObj / np.sqrt( NeObj + NeBkgd + NeDark + NeRon ) # Here's the actual S/N equation for t seconds.

        if( result=='e'):  # If exp time requested, and
            if(S_Nt < S_Nlim):  # If calculated S/N hasn't exceeded input S/N,
                if (resoln<10):
                    t += 1 # Then add 1 second, for all cameras except spectrographs,
                else:
                    t += 20 # For spectrographs, add 20 seconds
            else:
                endloop=1 # If input S/N reached, then quit.

        elif( result=='m'):  # If mag limit requested, and
            if (S_Nt < S_Nlim):
                maga -= 0.1 # If calculated S/N hasn't exceeded input S/N, then brighten by 0.1 mag
            else:
                endloop=1 # If input S/N reached, then quit.

        else:
            endloop=1 # If S/N requested, then quit.
    

    # Rounding
    S_Nt = np.round(10.0*S_Nt) 
    S_Nt = S_Nt/10.0
    mag = ((1000.0*maga)-acorr)/1000.0
    mag = np.round(10.0*mag) 
    mag = mag/10.0
    NeObj = np.round(10.0*NeObj) 
    NeObj = NeObj/10.0
    NeBkgd = np.round(10.0*NeBkgd) 
    NeBkgd = NeBkgd/10.0
    NeDark = np.round(10.0*NeDark) 
    NeDark = NeDark/10.0
    NeRon = np.round(10.0*NeRon) 
    NeRon = NeRon/10.0
    NbDN = np.round(10.0*NbDN) 
    NbDN = NbDN/10.0
    PkDN = np.round(10.0*PkDN) 
    PkDN = PkDN/10.0
    Npix = np.round(10.0*Npix) 
    Npix = Npix/10.0
    Nas = np.round(10.0*Nas) 
    Nas = Nas/10.0
    throughput = np.round(100.0*throughput) 
    throughput = throughput/100.0

    Saturated=0
    if (PkDN>40000):
        Saturated=1000000

    r = {# These are the values sent to the screen as output?
        'S_Nt': S_Nt,
        'mag': mag,
        't': t,
        'dqe': throughput,
        'apixel': apixel,
        'Nas': Nas,
        'Msky': Msky,
        'NbDN': NbDN,
        'PkDN': PkDN,
        'extco': extco,
        'NeObj': NeObj,
        'NeBkgd': NeBkgd,
        'NeDark': NeDark,
        'NeRon': NeRon,
        'Npix': Npix,
        'dark': dark,
        'ron': ron,
        'gain': gain,
        'adiam': adiam,
        'carea':  Carea[itel],
        'fcent': Fcent[ifilt],
        'saturated': saturated}
    
    return r
