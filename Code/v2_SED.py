# all imports 
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import v2_AGN_DataAndPlotting as adp

# ====================== COLORS ======================

# custom colormaps 
red_cmap = mpl.colors.ListedColormap(['#6A040F','#9D0208','#D00000','#DC2F02','#E85D04','#F48C06','#FAA307','#FFBA08'])
yel_cmap = mpl.colors.ListedColormap(['#3C8802','#6FB505','#A1D907','#D1F520','#FFF705','#F7E200','#EECA00','#DDB100'])
blu_cmap = mpl.colors.ListedColormap(['#003E6E','#0461A7','#0699C6','#0FCBD1','#24E3C0','#47F08B','#94F76F','#BEF768'])

# returns the colormap corresponding to the key 
def GetCmap(key='turbo', n_ticks=9) :
    # red
    if( key=='red' or 
        key=='r'
      ) : 
        return red_cmap
    # yellow
    if( key=='yel' or 
        key=='yellow' or 
        key=='y'
      ) : 
        return yel_cmap
    # blue 
    if( key=='blu' or 
        key=='blue' or 
        key=='b'
      ) : 
        return blu_cmap
    # rainbow (default)
    return plt.cm.get_cmap('turbo', n_ticks-1)

# plot colorbar 
def PlotColorbar(cmap, min, max, n_ticks, label):
    # get tick marks
    interval = (max - min) / (n_ticks - 1)
    mult = np.arange(0,n_ticks)
    ticks = min+(interval*mult)
    # setup colorbar 
    norm = mpl.colors.Normalize(vmin=min, vmax=max)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    # plot 
    plt.colorbar(sm, ticks=ticks, label=label, format='%.1f')

# normalizes data for cmap() in log scale 
def NormalizeForLogCmap(z) : 
    # normalize
    znorm = z / np.nanmax(abs(z))   
    # get log of normalized 24um
    znorm_log = np.log10(znorm) 
    # normalize the log-normalized-z and +1 to make positive range 0-1        
    znorm_log_norm = (znorm_log / np.nanmax(abs(znorm_log))) + 1    
    # return values for cmap() 
    return znorm_log_norm

# ====================== WAVELENGTHS ======================

# list wavelengths of the photometry filters (in Angstroms) 
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
    ]) # Used lambda_ref from: http://svo2.cab.inta-csic.es/theory/fps3/

# returns array of observed wavelengths in Angstroms
def GetObservedWavelengths_A(print=True) : 
    if(print): 
        PrintShape(lam_A)
    return lam_A

# ====================== HELPER FUNCTIONS ======================

# returns array of ID numbers
def GetID(data, print=True) :
    # get id numbers
    ids = np.array(data['ID_COSMOS2015'])
    if(print):
        PrintShape(ids)
    return ids

# prints shape of array
def PrintShape(arr) : 
    print('Array shape:\t', np.shape(arr))

# returns absolute maximum
def MaxMagOfUneven2DList(list) : 
    max_each = []
    for row in list : 
        max_each.append(max(abs(row)))
    return max(max_each)

# returns maximum
def MaxOfUneven2DList(list) : 
    max_each = []
    for row in list : 
        max_each.append(np.nanmax(row))
    return max(max_each)

# returns maximum
def MinOfUneven2DList(list) : 
    max_each = []
    for row in list : 
        max_each.append(np.nanmin(row))
    return min(max_each)

# ====================== MATH ======================

# unit conversions
def Convert_A_um(ang) : 
    return (ang * 1E-4)
def Convert_A_cm(ang) : 
    return (ang * 1E-8)

# returns interpolated function (log scale) from x and y (not log)
def Interpolate_log(x,y) :
    # get log scale and exclude NaN
    logx = np.log10(x[~np.isnan(y)])
    logy = np.log10(y[~np.isnan(y)])
    # interpolate curve 
    f = interpolate.interp1d(logx, logy, kind='linear', fill_value='extrapolate')
    # return interpolated function 
    return f

# returns f(x) corrected for log scales 
def Flog_X(f,x) : 
    return 10**f(np.log10(x))

# ====================== SED PREP ======================

# returns array of all observed photometry and their IDs. Bad values and measurements with fracErr are set to NaN.
def GetPhotometry(data, fracErr=0.5, print=True) :

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
    if(print):
        PrintShape(flux_nu_uJy)

    # return table of observed photometry with bad values set to nan
    return flux_nu_uJy

# Returns observed wavelength array
def ConvertToRestWavelength(redshifts, lamObs=lam_A, print=True) : 
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
    if(print):
        PrintShape(lamRest)
    # return array of rest wavelengths
    return lamRest

# returns lambda*F_lambda after converting F_nu
def ConvertToEnergyDensity(lamRest_A, Fnu_uJy, print=True) : 
    # comvert units from angrstroms to centimeters
    lamRest_cm = Convert_A_cm(lamRest_A)
    # c[cm/s]/lambda^2
    c_lamR2 = 2.99792458 * 1E10  / lamRest_cm**2
    # convert Jy to cgs units: 1 Jy = 10^-23 erg/s/cm2/Hz, uJy = 10^-6 Jy
    Fnu_ergscm2Hz = Fnu_uJy * 1E-23 * 1E-6
    # multiply F_nu by c/lam^2 to get F_lam
    Flam_ergscm3 = np.multiply(Fnu_ergscm2Hz, c_lamR2)
    # lam * F_lam [erg/s/cm]
    lamFlam_ergscm2 = np.multiply(Flam_ergscm3, lamRest_cm)
    # verify same shape 
    if(print):
        PrintShape(lamFlam_ergscm2)
    # return lamFlam energy density
    return lamFlam_ergscm2

# returns normalized SED at one micron
def NormalizeSED_1um(lamRest_A, lamFlam_ergscm2, print=True) :

    # convert angstrom to microns
    lamRest_um = Convert_A_um(lamRest_A)

    # initialize list
    lamFlam_ergscm2_NORM = []

    # interpolate and normalize at 1um 
    for x,y in zip(lamRest_um,lamFlam_ergscm2) : 
        # interpolate
        f = Interpolate_log(x,y)
        # normalize
        at1um = Flog_X(f,1)
        lamFlam_ergscm2_NORM.append( y / at1um)

    # convert list to array 
    lamFlam_ergscm2_NORM = np.array(lamFlam_ergscm2_NORM)

    # print shape 
    if(print):
        PrintShape(lamFlam_ergscm2_NORM)

    # return normalized lam*F_lam
    return lamFlam_ergscm2_NORM

# ====================== SINGLE SED PLOTTING ======================

# plot SED curves from x and y arrays 
def PlotSED(
        x,                  # x-axis data   lam [A]
        y,                  # y-axis data:  lamFlam [erg/s/cm2]
        cmap='',            # colormap options: red, yel, blu, (jet otherwise)
        n_ticks=9,          # number of ticks on colorbar
        showBar=True,       # show the colorbar 
        title='',           # plot title
        save='',            # filename to save
        median=True,        # plots a median line when true
        xmin=10**-2,        # plot range 
        xmax=10**3,         #   "    "
        ymin=10**-2.5,      #   "    "
        ymax=10**2.5        #   "    "
    ) : 
    # standardize plot
    adp.SetStyle()

    # prepare x
    x_um = Convert_A_um(x)
    n = np.shape(x)[0]

    # get colormap
    cmap_use = GetCmap(cmap, n_ticks)

    # prepare colormap values 
    z = y.T [14]        # transpose y to get 24um column
    z_log = np.log10(z)      
    znorm_log_norm = NormalizeForLogCmap(z)

    # get masks for no and yes 24um measurements 
    no24 = np.where(np.isnan(z))
    no24_mask = (np.zeros(n, dtype=bool))
    no24_mask[no24] = True 
    ye24_mask = ~ no24_mask

    # apply mask 
    x_no24 = x_um[no24_mask]
    y_no24 = y[no24_mask]
    x_ye24 = x_um[ye24_mask]
    y_ye24 = y[ye24_mask]
    znorm_log_norm_ye24 = znorm_log_norm[ye24_mask]

    # plot SED curves for no24
    n_no24 = np.shape(x_no24)[0] # row
    for i in range(n_no24) : 
        plt.plot(x_no24[i],y_no24[i], color='gray')

    # plot SED curves for ye24
    n_ye24 = np.shape(x_ye24)[0] # row
    for i in range(n_ye24) : 
        plt.plot(x_ye24[i],y_ye24[i], color=cmap_use(znorm_log_norm_ye24[i]))

    # setup colorbar 
    if(showBar) :
        PlotColorbar(
            cmap=cmap_use, 
            min=np.nanmin(z_log), 
            max=np.nanmax(z_log), 
            n_ticks=n_ticks, 
            label='$Normalized \; \lambda F_{\lambda} \; at \; 24 \mu m$'
        )
    
    # plot median
    if(median) : 
        x_m, y_m = MedianCurve(x_um, y, xmin=1E-1,xmax=1E+1)
        plt.plot(x_m,y_m,c='k',linewidth=2)

    # setup plot 
    PlotSED_Settings( 
        n=n, 
        title=title,    save=save, 
        xmin=xmin,      xmax=xmax, 
        ymin=ymin,      ymax=ymax 
    ) 

    # display 
    plt.show()

    # end 
    return   

# plot settings for SED
def PlotSED_Settings(
        n=0,                # number of sources for text 
        title='',           # plot title
        save='',            # filename to save
        xmin=10**-2,        # plot range 
        xmax=10**3,         #   "    "
        ymin=10**-2.5,      #   "    "
        ymax=10**2.5,       #   "    "
        setLabels=True      # sets x and y lables when true
    ) : 
    # scale
    plt.yscale('log')
    plt.xscale('log')
    plt.grid()
    # range
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    
    # axis
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_adjustable('box')
    ax.set_xticks([1E-2,1E-1,1E0,1E1,1E2,1E3])
    ax.set_yticks([1E-2,1E-1,1E0,1E1,1E2])

    # set labels 
    if(setLabels) : 
        plt.xlabel('$\lambda_{rest} \; [\mu m]$') 
        plt.ylabel('$Normalized \; \lambda F_{\lambda} \; [erg \; s^{-1} \; cm^{-2}]$')
    # set text
    if(n>0) : 
        adp.AddText_n(n)
    # set title
    if(title) : 
        plt.title(title)
    # save plot as image 
    if(save) : 
        adp.Save(save)

    # end
    return

# returns x and y arrays for the parameter's median curve
def MedianCurve(x,y,xmin=1E-1,xmax=1E+2) : 
    
    # initialize objects
    f_all = []
    f_all_discrete = []
    x_sample = np.arange(xmin,xmax,0.1)

    # interpolate each source
    for xx,yy in zip(x,y) : 
        f = Interpolate_log(xx,yy)
        f_all.append(f)

    # get discrete points for each f(x)
    for f in f_all : 
        discrcete = Flog_X(f,x_sample)
        f_all_discrete.append(discrcete)

    # get median y value for each x
    y_median = np.nanmedian(f_all_discrete, axis=0)

    # return x and y 
    return x_sample, y_median