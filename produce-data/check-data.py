#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 14:05:52 2022

@author: bas
"""

import dill as pickle
import astropy.units as u

read_data='shock'
dist='40'
bs_optical_string='ugri'
bs_uv_string='NUV_DD1D2'
delay=0


with open(f'../input_files/data/SNR_fiducial_{read_data}_{dist}Mpc_opticalbands_{bs_optical_string}_uvbands_{bs_uv_string}_{delay}h_delay.pkl','rb') as tf:
    ts_data, abmags_data, snrs, abmags_error=pickle.load(tf)
    
print('time data:',ts_data['r'].to(u.hr))
print('abmags data',abmags_data['r'])