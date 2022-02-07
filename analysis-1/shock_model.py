# k in 0.1 cm^2/g
# M in 0.01 solar masses
# v in 0.1c
#R_0 in 10^10 cm
import numpy as np


sigma_SB = 5.670373e-5 # g/s^3/K^4 
h= 6.626e-27 # erg s 
kb = 1.38e-16 # erg / K
c = 2.99792458e10 # cm / s
Mpc = 3.08568e24 # cm

def Lum_f(t,k,M,v,R_0):
    L = 9.5e40*k**-0.6*M**0.4*v**(8/5)*R_0*t**-0.8 #erg/s
    return(L)

def Teff_f(t,k,M,v,R_0):
    T = 6.2e3*k**(-7/30)*M**(1/60)*v**(1/15)*R_0**0.25*t**(-8/15) #K
    return(T)

def R_ph_f(t,k,M,v,R_0):
    R_ph = 3e14*k**(1/6)*M**(1/6)*v**(2/3)*t**(2/3) #cm
    return(R_ph)

def B_f(T,f):
    B = 2*h*f**3/c**2/(np.exp(h*f/(kb*T))-1)
    return(B)

def get_light_curve(t,k,M,v,R_0,D,zr):
    D = D*Mpc
    L = 9.5e40*k**-0.6*M**0.4*v**(8/5)*R_0*t**-0.8 #erg/s
    T = 6.2e3*k**(-7/30)*M**(1/60)*v**(1/15)*R_0**0.25*t**(-8/15) #K
    R_ph = 3e14*k**(1/6)*M**(1/6)*v**(2/3)*t**(2/3) #cm
    f_dorado = (1+zr)*c/(0.2036e-4) #pivot wavelength 2036 amstrong
    corr_dorado = np.pi*B_f(T,f_dorado)/(sigma_SB*T**4) #correction factor
    f_dorado=corr_dorado*L/(4*np.pi*D**2)/(1+zr) #erg/s/cm^2/Hz
    mag_dorado = -2.5*np.log10(f_dorado)-48.6 #AB mag
    return(L,T,R_ph,mag_dorado)
    
    