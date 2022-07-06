import matplotlib.pyplot as plt
import plotting as myP
from scipy import interpolate
from astropy.io import fits
import numpy as np

# list wavelengths of the photometry filters (in Angstroms) -- Used lambda_ref from: http://svo2.cab.inta-csic.es/theory/fps3/
lam_A = np.array([
          1535.08,      # 'GALEX_FUV_FLUX'
          2300.79,      # 'GALEX_NUV_FLUX'
          4798.21,      # 'HSC_g_FLUX_APER2'
          6218.44,      # 'HSC_r_FLUX_APER2'
          7727.02,      # 'HSC_i_FLUX_APER2'
          8908.50,      # 'HSC_z_FLUX_APER2'
          9775.07,      # 'HSC_y_FLUX_APER2'
         12524.83,      # 'UVISTA_J_FLUX_APER2'
         16432.40,      # 'UVISTA_H_FLUX_APER2'
         21521.52,      # 'UVISTA_Ks_FLUX_APER2'
         35378.41,      # 'SPLASH_CH1_FLUX' 
         44780.49,      # 'SPLASH_CH2_FLUX' 
         56961.77,      # 'SPLASH_CH3_FLUX' 
         77978.39,      # 'SPLASH_CH4_FLUX' 
        235934.61,      # 'FIR_24_FLUX'      
             1.01e+6,   # 'FIR_100_FLUX'
             1.62e+6,   # 'FIR_160_FLUX'
             2.48e+6,   # 'FIR_250_FLUX'
             3.48e+6,   # 'FIR_350_FLUX'
             5.00e+6    # 'FIR_500_FLUX'
    ])

# returns array of observed wavelengths in Angstroms
def GetObservedWavelengths_A() : 
    PrintShape(lam_A)
    return lam_A

# prints shape of array
def PrintShape(arr) : 
    print('Array shape:\t', np.shape(arr))

# returns array of all observed photometry and their IDs. Bad values and measurements with fracErr are set to NaN.
def GetPhotometry(data, fracErr=0.5) :

    # get photometry table values
    flux_nu_uJy = np.array([   
        # The Galaxy Evolution Explorer (GALEX)
        np.array(data['GALEX_FUV_FLUX']), 
        np.array(data['GALEX_NUV_FLUX']),
        # Subaru / Hyper Suprime-Cam
        np.array(data['HSC_g_FLUX_APER2']),
        np.array(data['HSC_r_FLUX_APER2']),
        np.array(data['HSC_i_FLUX_APER2']),
        np.array(data['HSC_z_FLUX_APER2']),
        np.array(data['HSC_y_FLUX_APER2']),
        # Paranal VISTA / VIRCAM
        np.array(data['UVISTA_J_FLUX_APER2']),
        np.array(data['UVISTA_H_FLUX_APER2']),
        np.array(data['UVISTA_Ks_FLUX_APER2']),
        # Spitzer Large Area Survey with Hyper-Suprime-Cam (SPLASH) 
        np.array(data['SPLASH_CH1_FLUX']),
        np.array(data['SPLASH_CH2_FLUX']),
        np.array(data['SPLASH_CH3_FLUX']),
        np.array(data['SPLASH_CH4_FLUX']),
        # FIR
        np.array(data['FIR_24_FLUX']),    # Spitzer/MIPS.24mu
        np.array(data['FIR_100_FLUX']),   # Herschel/Pacs.green
        np.array(data['FIR_160_FLUX']),   # Herschel/Pacs.red
        np.array(data['FIR_250_FLUX']),   # Herschel/SPIRE.PSW
        np.array(data['FIR_350_FLUX']),   # Herschel/SPIRE.PMW
        np.array(data['FIR_500_FLUX'])    # Herschel/SPIRE.PLW
    ])  

    # get photometry table values
    flux_err_nu_uJy = np.array([   
        # The Galaxy Evolution Explorer (GALEX)
        np.array(data['GALEX_FUV_FLUXERR']), 
        np.array(data['GALEX_NUV_FLUXERR']),
        # Subaru / Hyper Suprime-Cam
        np.array(data['HSC_g_FLUXERR_APER2']),
        np.array(data['HSC_r_FLUXERR_APER2']),
        np.array(data['HSC_i_FLUXERR_APER2']),
        np.array(data['HSC_z_FLUXERR_APER2']),
        np.array(data['HSC_y_FLUXERR_APER2']),
        # Paranal VISTA / VIRCAM
        np.array(data['UVISTA_J_FLUXERR_APER2']),
        np.array(data['UVISTA_H_FLUXERR_APER2']),
        np.array(data['UVISTA_Ks_FLUXERR_APER2']),
        # Spitzer Large Area Survey with Hyper-Suprime-Cam (SPLASH) 
        np.array(data['SPLASH_CH1_FLUXERR']),
        np.array(data['SPLASH_CH2_FLUXERR']),
        np.array(data['SPLASH_CH3_FLUXERR']),
        np.array(data['SPLASH_CH4_FLUXERR']),
        # FIR
        np.array(data['FIR_24_FLUXERR']),    # Spitzer/MIPS.24mu
        np.array(data['FIR_100_FLUXERR']),   # Herschel/Pacs.green
        np.array(data['FIR_160_FLUXERR']),   # Herschel/Pacs.red
        np.array(data['FIR_250_FLUXERR']),   # Herschel/SPIRE.PSW
        np.array(data['FIR_350_FLUXERR']),   # Herschel/SPIRE.PMW
        np.array(data['FIR_500_FLUXERR'])    # Herschel/SPIRE.PLW
    ]) 

    # transpose so [i] is a source, not a column
    flux_nu_uJy = flux_nu_uJy.T
    flux_err_nu_uJy = flux_err_nu_uJy.T

    # change bad values to NaN
    flux_nu_uJy = np.where(flux_nu_uJy <= float(0), float('nan'), flux_nu_uJy)
    flux_nu_uJy = np.where(flux_err_nu_uJy/flux_nu_uJy >= fracErr, float('nan'), flux_nu_uJy) 

    # print info
    PrintShape(flux_nu_uJy)

    # return table of observed photometry with bad values set to nan
    return flux_nu_uJy

# returns array of ID numbers
def GetID(data) :
    # get id numbers
    ids = np.array(data['ID_COSMOS2015'])
    PrintShape(ids)
    return ids

# Returns observed wavelength array
def ConvertToRestWavelength(redshifts, lamObs=lam_A) : 
    # get table shape
    row = np.shape(redshifts) [0]
    col = np.shape(lamObs) [0]

    # initialize table of rest frame wavelengths
    lamRest = np.zeros( (row,col) )

    # get rest wavelength 
    for i,z in enumerate(redshifts) : 
        # lambda_rest = lambda_observed / (1+z)
        lamRest[i] = lamObs / (1+z)

    # print info
    PrintShape(lamRest)

    # return array of rest wavelengths
    return lamRest

# returns lambda*F_lambda after converting F_nu
def ConvertToEnergyDensity(lamRest_A, Fnu_uJy) : 

    # comvert units from angrstroms to centimeters
    lamRest_cm = lamRest_A * 1E-8

    # c[cm/s]/lambda^2
    c_lamR2 = 2.99792458 * 1E10  / lamRest_cm**2

    # convert Jy to cgs units: 1 Jy = 10^-23 erg/s/cm2/Hz, uJy = 10^-6 Jy
    Fnu_ergscm2Hz = Fnu_uJy * 1E-23 * 1E-6

    # multiply F_nu by c/lam^2 to get F_lam
    Flam_ergscm3 = np.multiply(Fnu_ergscm2Hz, c_lamR2)

    # lam * F_lam [erg/s/cm]
    lamFlam_ergscm2 = np.multiply(Flam_ergscm3, lamRest_cm)

    # verify same shape 
    PrintShape(lamFlam_ergscm2)

    # return lamFlam energy density
    return lamFlam_ergscm2

# returns normalized SED at one micron
def NormalizeSED_1um(lamRest_A, lamFlam_ergscm2) :

    # convert angrstrom to microns
    lamRest_um = lamRest_A * 1E-4

    # initialize list
    lamFlam_ergscm2_NORM = []

    for row in range(np.shape(lamRest_um)[0]) : 
        # get x and y list
        x = lamRest_um[row]
        y = lamFlam_ergscm2[row]
        # get log scale and exclude NaN
        logx = np.log10(x[~np.isnan(y)])
        logy = np.log10(y[~np.isnan(y)])
        # interpolate flux curve 
        f = interpolate.interp1d(logx, logy, kind='linear', fill_value='extrapolate')
        # normalize
        at1um = 10**f(np.log10(1))
        lamFlam_ergscm2_NORM.append( y / at1um)

    # convert list to array 
    lamFlam_ergscm2_NORM = np.array(lamFlam_ergscm2_NORM)

    # print shape 
    PrintShape(lamFlam_ergscm2_NORM)

    # return normalized lam*F_lam
    return lamFlam_ergscm2_NORM
    
# plot SED curves from x and y arrays 
def PlotSED(
        x,                  # x-axis data   lam
        y,                  # y-axis data:  lamFlam
        num=0,              # number of SED curves to plot 
        offset=0,           # index offset when plotting SED curves
        xunit = '$\mu m$',  # unit for x-axis 
        title='',           # plot title
        save=''             # filename to save
    ) : 

    # get number of sources to plot 
    n = 0
    if(num) :   n = num
    else :      n = np.shape(x)[0] # row

    # plot SED curves
    for i in range(n) : 
        plt.plot(x[i+offset],y[i+offset])

    # plot settings
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(10**-2.5, 10**3.5)
    plt.ylim(10**-3.5, 10**3.5)
    plt.xlabel('$\lambda_{rest} \; [$'+xunit+'$]$') 
    plt.ylabel('$Normalized \; \lambda F_{\lambda} \; [erg \; s^{-1} \; cm^{-2}]$')
    plt.grid()
    myP.addtext_n(n)

    # square axis
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_adjustable('box')

    # set title
    if(title!='') : 
        plt.title(title)

    # save plot as image 
    if(save) : 
        myP.save(save)
    
    # display 
    plt.show()   
