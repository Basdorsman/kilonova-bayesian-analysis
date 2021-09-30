import numpy as np
import astropy.constants as const
import synphot

h = const.h.cgs
c = const.c.cgs
k_B = const.k_B.cgs

def B_of_lambda(l,T):
    B=2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))
    return(B)

def BB_curves(T,n_l,l_min,l_max):
    l = np.linspace(l_min,l_max,n_l)
    B = B_of_lambda(l,T)
    return(l,B)

def Lum_Band(T,R,min_lambda,max_lambda):
    def B_of_lambda_int(l):
        B = B_of_lambda(l,T)
        return(B)
    value,Error = quad(B_of_lambda_int,min_lambda,max_lambda)
    return((4*np.pi*R)**2*value)
    
def Lum(R,B):
    L = (4*np.pi*R)**2*B
    return(L)
 
def R_ph(L,T):
    R = np.sqrt(L/(4*np.pi*sigma*T**4))
    return(R)

def Flux_Energy(Distance,L):
    F = L/(Distance**2*4*np.pi)
    return(F)
    
def Flux_Photon(F_E,wavelength):
    F_p = F_E*wavelength/(c*h)
    return(F_p)

def mag_AB_error(SNR):
    e = np.where(SNR > 0.0001, 2.5*np.log10(1+1/SNR), 10) #if SNR is very low, error is large
    return(e)

from astropy import units as u

def B_of_nu(T,nu): #input with astropy.units please
    h = 6.626e-27*u.erg*u.s
    k = 1.38e-16*u.erg/u.K
    c = 2.99792458e10*u.cm/u.s
    return(2.0*h*nu**3.0)/(c*c)/(np.exp(h*nu/(k*T))-1)

def get_abmag_pivot(D,lambda_pivot,LBol,T): #all in astropy units
    sigma_SB = 5.670373e-5*u.erg/u.cm**2/u.s/u.K**4 #erg⋅cm−2⋅s−1⋅K−4.
    c = 2.99792458e10*u.cm/u.s
    zr = get_redshift(D)
    nu = (1.+zr)*c/lambda_pivot # observed wavelength to emmitted frequency
    corr = np.pi*B_of_nu(T,nu)/(sigma_SB*T**4) #spectral density B / emmited power per surface
    #print('L',LBol)
    #print('T',T)
    #print('corr',corr)
    f = corr*LBol/(4*np.pi*D*D)/(1.+zr) #observed flux f_nu (f_nu since B_nu was used)
    f = f.to(u.erg/u.cm**2) #removes a factor of some power of cm/Angstrom 
    #print('flux',f)
    ABmag = -2.5*np.log10(f.value)-48.6
    return(ABmag)

def get_flux_nu(D,lambda_pivot,LBol,T): #Expecting D, LBol, T in astropy.units
    sigma_SB = 5.670373e-5*u.erg/u.cm**2/u.s/u.K**4 #erg⋅cm−2⋅s−1⋅K−4.
    c = 2.99792458e10*u.cm/u.s
    zr = get_redshift(D)
    nu = (1.+zr)*c/lambda_pivot # observed wavelength to emmitted frequency
    corr = np.pi*B_of_nu(T,nu)/(sigma_SB*T**4) #spectral density B / emmited power per surface
    f = corr*LBol/(4*np.pi*D*D)/(1.+zr) #observed flux f_nu (f_nu since B_nu was used)
    f = f.to(u.erg/u.cm**2)
    return(f)

from astropy.cosmology import Planck13,z_at_value

def get_redshift(distance): #distance in u.Mpc
    return(z_at_value(Planck13.luminosity_distance,distance))

def correct_mass_redshift(mass,redshift):
    m_corrected = mass/(1+redshift)
    return(m_corrected)

def get_abmag(T,r,distance,bandpass):
    '''calculate AB magnitude faster.
    
    This code is largely copied from Synphot, using: synphot.BlackBody1D,
    synphot.SourceSpectrum and synphot.Observation. However, it's vectorized. 
    For multiple bandpasses, I use loops. Redshift has not yet been 
    implemented here.

    Parameters
    ----------
    T : 1darray * quantity temperature
        Blackbody Temperature
    r : 1darray * quantity length
        Blackbody radius
    distance : quantity length
        Luminosity distance
    bandpass : object synphot.spectrum.SpectralElement
    
    Returns
    -------
    abmag : list of 1darrays 
        list of AB magnitudes. One entry per bandpass
    '''
   
    wav = bandpass.waveset
    wav_nounits = wav.value

    # Calculate spectral radiance : np.array(len(T),len(wav))
    B_l = 2*h*c**2/np.multiply(wav**5,np.exp(h*c/(k_B*np.outer(T,wav)))-1)

    # Scale spectral density = flux
    f_l = B_l.T*np.pi*(r / distance)**2 

    # Convert flux to photlam and multiply by bandpass
    f_PHOTLAM = (f_l.T/((h*c/wav))).to(1/u.cm**2/u.s/u.Angstrom)*bandpass(wav) 
    num = abs(np.trapz(f_PHOTLAM.value * wav_nounits, x=wav_nounits)) 
    den = np.ones(len(T))*abs(np.trapz(bandpass(wav) * wav_nounits,
                                       x= wav_nounits))

    # Convolve with bandpass for total photon flux
    val = (num / den)*synphot.units.PHOTLAM

    # Convert flux from PHOTLAM to AB magnitude.
    abmag = synphot.units.convert_flux(bandpass.pivot(),val,u.ABmag).value
    return abmag

def get_abmag_synphot(T, r, distance, bandpass):
    seds = [synphot.SourceSpectrum(synphot.BlackBody1D, temperature=TT) * 
            np.pi * (rr / distance).to(u.dimensionless_unscaled)**2 for TT,rr
            in zip(T, r)]
    abmag_synphot = [synphot.Observation(sed, bandpass).effstim(u.ABmag).value
                     for sed in seds]
    return abmag_synphot
