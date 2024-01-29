# all imports
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import csv 
from rasterize_and_save import rasterize_and_save

# color constants
c_ir        = '#D00000' # use this for IR sources                   ALT: '#CE2727' # Fire Engine Red 
c_xray      = '#0461A7' # use this for Xray sources                 ALT: '#384E77' # Y In Mn Blue    
c_ir_xray   = '#81C138' # use this for matched IR and Xray sources  ALT: '#F6CD13' # Jonquil       # OLD: '#D1F520'

# path variables 
path_cat = 'D:\\IfA REU\\Catalogs'
path_fig = 'D:\\IfA REU\\IfA-REU-Project\\Figures'
path_csv = 'D:\\IfA REU\\IfA-REU-Project\\Lists'

# text sizes
BIG = 16
MID = 12
SML = 8

##### File IO #####

def UnpackTable(fname: str) -> dict[str,np.array]:
    """Get a table dictionary like {x: xdata, y: ydata} from a CSV file containing data""" 
    # get data from the csv file
    arr = np.genfromtxt(fname, delimiter=",", dtype=str)
    # separate column names and data arrays 
    cols = np.array(arr[0 ], dtype=str)
    data = np.array(arr[1:]).T # transpose to get x and y
    # build dictionary to hold the table {x: xdata, y: ydata}
    table = {}
    for col, dat in zip(cols,data) : 
        table[col] = dat
    return table

# read from a given csv file and return a list of its contents 
def ReadFile(filename, crop=True) : 
    # initialize list
    file = []

    # open csv file
    with open(filename, newline='') as csv_file : 
        reader = csv.reader(csv_file, delimiter=',', quotechar='"')
        # output each line in the file to the list
        for row in reader :
            file.append(row)
    
    if(crop) : 
        print('Cropped: ', file[0])
        # return an array of type int with column names clipped off
        return(file[1:])

    # try using these in the future for reading files!
    # * numpy - genfromtext can remove headers
    # * pandas - dataframes good at reading csv files and removing headers
    # * astropy ascii - ascii.read()

    # else, return the full list
    return(file)

def ReadRedshifts(file, printInfo=True) : 
    # read redshifts from file 
    zAll = np.array(ReadFile(file))
    zID, zZ = zAll.T
    # get correct data type
    zID = np.array(zID, dtype=int)
    zZ  = np.array(zZ,  dtype=float)
    # apply mask to remove bad redshifts
    mask_z = (zZ >= 0) & (zZ < 99)
    zID = zID[mask_z]
    zZ  = zZ[mask_z]
    # print info
    if(printInfo):
        print('Number of redshifts:', len(zID))
    # build dict
    zdict = {
        'ID' : zID,
        'z' : zZ
    }
    # return dictionary of valid redshifts and IDs
    return(zdict)

# returns array of ( min <= z < max )
def GetRedshiftsInRange(
        zDict,   # dictionary of ID and z  
        min=0,   # minimum
        max=99,  # maximum 
        incusiveMin=False,
        printLength=True
) : 
    """
    zDict_split = {
        '$0 < z \leq 1$'   : adp.GetRedshiftsInRange(zDict_all,          max=1  ),
        '$1 < z \leq 1.5$' : adp.GetRedshiftsInRange(zDict_all, min=1,   max=1.5),
        '$1.5 < z \leq 2$' : adp.GetRedshiftsInRange(zDict_all, min=1.5, max=2  ),
        '$2 < z \leq 3$'   : adp.GetRedshiftsInRange(zDict_all, min=2,   max=3  ),
        '$3 < z \leq 6$'   : adp.GetRedshiftsInRange(zDict_all, min=3           )
    }
    """
    # get mask of indecies that are in redshift range (true)
    if(incusiveMin):
        mask = (zDict['z'] >= min) & (zDict['z'] < max)
    else: 
        mask = (zDict['z'] > min) & (zDict['z'] <= max)
    # idolate ID and z
    inRange = {
        'ID' : zDict['ID'][mask] ,
        'z'  : zDict['z'] [mask]
    }
    if(printLength) : 
        if(incusiveMin): 
            print(min, '<= z < ', max, ': ', len(inRange['ID']))
        else:
             print(min, '< z <= ', max, ': ', len(inRange['ID']))       
    return(inRange)



##### General Plotting #####

# redefine the defaults for plots
def SetStyle() :     
    # figure
    mpl.rcParams['figure.figsize']  = 6, 4  # figure size in inches
    mpl.rcParams['figure.dpi']      = 150   # dots-per-inch
    # axes
    mpl.rcParams['axes.titlesize']  = BIG   # size of title
    mpl.rcParams['axes.titlepad']   = BIG   # space between title and plot 
    mpl.rcParams['axes.labelsize']  = MID   # font size of the x and y labels
    mpl.rcParams['axes.labelpad']   = SML   # space between label and axis
    # ticks
    mpl.rcParams['xtick.labelsize'] = SML   # font size of the tick labels
    mpl.rcParams['ytick.labelsize'] = SML   # font size of the tick labels
    # legend
    mpl.rcParams['legend.fontsize'] = SML   # font size of the legend lables  
    # lines 
    mpl.rcParams['lines.linewidth'] = 1     # line width in points

# darkens a color
def darken_color(color, amount=1.4) :
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

# add 'n = #' to bottom right of plot
def AddText_n(n, pre='n = ', fontsize=8):
    ax = plt.gca()
    AddText_n_ax(ax, n, pre=pre, fontsize=fontsize)

# adds n=# to bottom right of ax
def AddText_n_ax(ax, n, pre='n = ', fontsize=SML):
    ax.text(    0.95,                           # x
                0.05,                           # y 
                pre + str(n),                   # string
                transform=ax.transAxes,         # use axis coordinants
                horizontalalignment='right',    # alignment 
                fontsize=fontsize               # font size
        )
    
def AddText_n2_ax(ax, n, pre='n = ', fontsize=SML, color='k'):
    ax.text(    0.95,                           # x
                0.10,                           # y 
                pre + str(n),                   # string
                transform=ax.transAxes,         # use axis coordinants
                horizontalalignment='right',    # alignment 
                color=color,                    # font color 
                fontsize=fontsize               # font size
        )

# adds redshift range to top left of ax
def AddText_z_ax(ax, fullText='', min=-1, max=-1, greaterEqual=False, lessEqual=True, fontsize=SML) :
    # initialize string
    text = ''
    # if full text is given
    if(fullText!='') : 
        text = fullText
    # otherwise, build text from parameters 
    else :
        if(min != -1) : # has min
            text += str(min)
            if(greaterEqual) :  text += ' $\leq$ '
            else :              text += ' < '
        text += 'z'
        if(max != -1) : # has max
            if(lessEqual):      text += ' $\leq$ '
            else:               text += ' < '
            text += str(max)
    # add text to top left on axis plot
    ax.text(    0.05,                           # x
                0.93,                           # y 
                text,                           # string
                transform=ax.transAxes,         # use axis coordinants
                horizontalalignment='left',     # alignment 
                fontsize=fontsize               # font size
        )

# draws median line for dataset x for a histogram axis
def MeanLineForHist_ax(ax,x,c='k',xtext=0.998, ytext=0.94, horizAlign='right', pre='Mean: ') :
    mean = np.array(x).mean()
    min_ylim, max_ylim = ax.set_ylim()
    ax.axvline(mean, color=c, linestyle='dashed', linewidth=2)
    ax.text(mean*xtext, max_ylim*ytext, pre+'{:.1f}'.format(mean), c=c, horizontalalignment=horizAlign)

# draws median line for dataset x for a histogram plot
def MeanLineForHist(x,c='k',xtext=0.998, ytext=0.94, horizAlign='right', pre='Mean: ') :
    ax = plt.gca()
    MeanLineForHist_ax(ax,x,c,xtext,ytext,horizAlign,pre)

# save plot 
def Save(filename, rasterize_list: list=None, dpi=1200) : 
    # get axis 
    rast: list  = rasterize_list
    if(rast == None) : rast = plt.gcf().get_axes() 
    # save
    rasterize_and_save(
        fname = filename, 
        rasterize_list = rast, 
        fig = plt.gcf(), 
        dpi = dpi,
        savefig_kw = {
            'bbox_inches'   : "tight",
            'pad_inches'    : 0.05,
            'facecolor'     : 'w',
            'edgecolor'     : 'w',
            'dpi'           : dpi
        }                   
    )

def BoldSubplot(ax, lineWidth = 2) : 
    ax.spines["top"     ].set_linewidth(lineWidth)
    ax.spines["bottom"  ].set_linewidth(lineWidth)
    ax.spines["left"    ].set_linewidth(lineWidth)
    ax.spines["right"   ].set_linewidth(lineWidth)


def PlotContours_ax(ax,x,y,colors='k',linewidths=0.4) :
    # modified from https://saturncloud.io/blog/scatter-plot-with-contour-according-to-a-parameter-using-matplotlib/ 
    from scipy.stats import kde
    # estimate density of points
    density = kde.gaussian_kde(np.vstack([x, y]))
    # create grid for the contour plot
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    Xgrid, Ygrid = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    Z = density.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))
    # plot density contour lines
    ax.contour(Xgrid, Ygrid, Z.reshape(Xgrid.shape), colors=colors,linewidths=linewidths)

##### IRAC color-color diagrams #####

# Draw the selection wedge for Donley 2012 on an axis  
def PlotDonleyWedge_ax(ax, linewidth=1) : 
    # constants
    x_min = 0.08    # x >= 0.08
    y_min = 0.15    # y >= 0.15
    max = 10        # arbritrary 

    # calculate intercepts
    x_int_ymin = (y_min + 0.27)/1.21  # intercept between y_min and (y>=1.21x-0.27)
    y_int_xmin = (1.21*x_min) + 0.27  # intercept between x_min and (y>=1.21x+0.27)

    # calculate y_low intercepts (y>=1.21x-0.27)
    y1_low = (1.21*x_int_ymin)  - 0.27   
    y2_low = (1.21*max)         - 0.27

    # calculate y_high intercepts (y>=1.21x+0.27)
    y1_high = (1.21*x_min)      + 0.27
    y2_high = (1.21*max)        + 0.27

    # plot lines between intercepts 
    ax.plot( [x_min,       x_min],         [y_min,     y_int_xmin], 'k', linewidth=linewidth)    # x >= 0.08
    ax.plot( [x_min,       x_int_ymin],    [y_min,     y_min],      'k', linewidth=linewidth)    # y >= 0.15
    ax.plot( [x_int_ymin,  max],           [y1_low,    y2_low],     'k', linewidth=linewidth)    # y >= 1.21x - 0.27
    ax.plot( [x_min,       max],           [y1_high,   y2_high],    'k', linewidth=linewidth)    # y <= 1.21x + 0.27

# Draw the selection wedge for Donley 2012 on an plot
def PlotDonleyWedge(linewidth=1) : 
    ax = plt.gca()
    PlotDonleyWedge_ax(ax,linewidth)
    plt.xlabel('$\log(f_{5.8um}/f_{3.6um})$')
    plt.ylabel('$\log(f_{8.0um}/f_{4.5um})$')

# calculate x position using IRAC channels
def IRACx(ch1, ch3):
    return(np.log10(ch3/ch1))   # x = log10 ( f_5.6um / f_3.6um )

# calculate y position using IRAC channels
def IRACy(ch2,ch4):
    return(np.log10(ch4/ch2))    # y = log10 ( f_8.0um / f_4.5 um )


##### Xray AGN selection #####

# correct Luminosity for absorbtion
def IntrinsicLuminosity(Lx,k) :
    # correct Luminosity for absorbtion: k_abs = L_abs / L_int --> L_int = L_abs / k_abs 
    luminosity = []
    for lum,abs in zip(Lx,k) : 
        if( lum <= 0) : 
            luminosity.append(-99)
        elif(abs <= 0) : # lum > 0
            luminosity.append(lum)      # assume this source is unobscured 
        else : # lum > 0 and abs > =
            luminosity.append( lum - np.log10(abs))
    return(np.array(luminosity))

# retuns a bool mask that is true if a source is an AGN (has a Lx>agnIf)
def AGN_Xray(L_xray, agnIf=43) : 
    # get number of sources (assume Lx and k are same )
    numSources = len(L_xray)
    # initialize mask of False
    AGNmask = np.zeros( numSources, dtype=bool)
    # a souce is an AGN if luminosity is greater than threshold
    for i,lum in enumerate(L_xray) : 
        if(lum >= agnIf) :
            AGNmask[i] = True 
    # return boolean mask (same size as param) that is True for AGN. 
    return(AGNmask)


##### MIR AGN selection #####

# returns true if the IRAC (x,y) position is inside the Donley et. al. 2012 wedge
def InDonleyWedge(x,y):
    if( x>=0.08 and 
        y>=0.15 and 
        y>=(1.21*x-0.27) and 
        y<=(1.21*x+0.27)
    ) :  
        return(True) 
    else :
        return(False)

# returns true if IRAC channels are monatimically rising (aka. always increacing)
def IsMonatomicallRising(
    f36,    # IRAC Ch 1
    f45,    # IRAC Ch 2
    f58,    # IRAC Ch 3
    f80     # IRAC Ch 4
):
    if( f45 > f36 and
        f58 > f45 and
        f80 > f58
    ):
        return(True)
    else :
        return(False)
    
# returns a bool mask that is true if source is AGN (meets Donley et. al. 2012 criterea)
def AGN_DonleyCriterea(
        f36,    # IRAC Ch 1
        f45,    # IRAC Ch 2
        f58,    # IRAC Ch 3
        f80     # IRAC Ch 4
    ) : 
    # get number of sources (assume all IRAC channels are same length)
    numSources = len(f36)
    # initialize mask of False
    AGNmask = np.zeros( numSources, dtype=bool)
    # look at each source 
    for i in range(numSources): 
        # must be nonzero
        if ((f36[i] <= 0) or
            (f45[i] <= 0) or
            (f58[i] <= 0) or
            (f80[i] <= 0)   ):
                continue 
        # must fall inside wedge and be monatomically rising
        x = np.log10(f58[i]/f36[i])
        y = np.log10(f80[i]/f45[i])
        if( InDonleyWedge(x,y) and 
            IsMonatomicallRising(f36[i],f45[i],f58[i],f80[i]) 
        ):
            AGNmask[i] = True 
    # return mask where True is an AGN 
    return(AGNmask)

##### By Z multipanel plotting #####

# generate figure and axis for n subplots oriented in row
def ByZ_SetupFig(n,orientVertical=True, figsizeConstVal=6):
    # standardize figure style
    SetStyle() 
    # calculate size using number of subplots
    if(orientVertical):
        figsize = (figsizeConstVal, 1+(3*n))
        nrows = n
        ncols = 1
    else : 
        figsize = (1+(3*n),figsizeConstVal)
        nrows = 1
        ncols = n
    # create figure and axis
    fig, ax = plt.subplots(nrows=nrows,ncols=ncols,sharex=True,sharey=True,figsize=figsize,layout='constrained',facecolor='w')
    return fig,ax

# finish by z plots
def ByZ_FinishPlot(
        fig,
        xaxis='',
        yaxis='',
        save='',
        xpos=(0.5, -0.01),
        ypos=(0.15,  0.5),
        fontsize=MID
):
    # name the x and y axis 
    if(xaxis) : fig.text(xpos[0], xpos[1], xaxis, ha='center', fontsize=fontsize)
    if(yaxis) : fig.text(ypos[0], ypos[1], yaxis, va='center', fontsize=fontsize, rotation='vertical')
    # finish
    if(save): Save(save)
    plt.show()
    plt.close()

def ByZ_SetupFig_Rectangle(nrow, ncol) : 
    # standardize figure style
    SetStyle() 
    # create figure and axis
    fig, ax = plt.subplots(nrows=nrow,ncols=ncol,sharex=True,sharey=True,figsize=((3*ncol),(3*nrow)),layout='constrained',facecolor='w')
    return fig, ax

##### Conversions #####

def LogErgsToLogSols(logLx_erg: float) -> float :
    """Convert log-scale Xray luminosity in ergs/s to log-scale solar luminosities 

    Args:
        logLx_erg (float): log-scale Xray luminosity in ergs/s

    Returns:
        float: log-scale solar luminosities
    """
    Lx_erg = 10**logLx_erg
    Lx_sol = Lx_erg / (3.826 * 10**33)
    logLx_sol = np.log10(Lx_sol)
    return logLx_sol

def LogSolsToLogErgs(logLx_sol: float)  -> float : 
    """Convert log-scale solar luminosities to log-scale Xray luminosity in ergs/s 

    Args:
        logLx_sol (float): log-scale solar luminosities

    Returns:
        float: log-scale Xray luminosity in ergs/s
    """
    Lx_sol = 10**logLx_sol
    Lx_erg = Lx_sol * (3.826 * 10**33)
    logLx_erg = np.log10(Lx_erg)
    return logLx_erg

def Lx_to_kx(logLx_sol: float, a: float = 15.33, b: float = 11.48, c: float = 16.20)  -> float : 
    """Calculate X-ray correction factor from X-ray luminosity band. 
    See Equation 3 in https://ui.adsabs.harvard.edu/abs/2020A%26A...636A..73D/abstract

    Args:
        logLx_sol (float): log-scale X-ray luminosity in solar luminosities.
        a (float, optional): See Table 1 row 'general' column a under Kx(Lx). Defaults to 15.33.
        b (float, optional): See Table 1 row 'general' column b under Kx(Lx). Defaults to 11.48.
        c (float, optional): See Table 1 row 'general' column c under Kx(Lx). Defaults to 16.20.

    Returns:
        float: X-ray correction factor
    """
    return ( a * (1.0 + (logLx_sol/b)**c)  )

def Lbol_to_kx(logLbol_sol: float, a: float = 10.96, b: float = 11.93, c: float = 17.79) -> float : 
    """Calculate X-ray correction factor from bolometric luminosity in solar luminosities. 
    See Equation 2 in https://ui.adsabs.harvard.edu/abs/2020A%26A...636A..73D/abstract

    Args:
        logLbol_sol (float): _description_
        a (float, optional): See Table 1 row 'general' column a under Kx(Lbol). Defaults to 10.96.
        b (float, optional): See Table 1 row 'general' column b under Kx(Lbol). Defaults to 11.93.
        c (float, optional): See Table 1 row 'general' column c under Kx(Lbol). Defaults to 17.79.

    Returns:
        float: X-ray correction factor
    """
    return ( a * (1.0 + (logLbol_sol/b)**c)  )

def Lx_To_Lbol(logLx_erg: float)  -> float : 
    """Calculate bolometric luminosity from correction factor.

    Args:
        logLx_erg (float): log-scale X-ray luminosity in erg/s

    Returns:
        float: Log-scale Bolometric luminosity in solar luminosities.
    """
    logLx_sol = LogErgsToLogSols(logLx_erg)
    kx = Lx_to_kx(logLx_sol)
    Lbol_sol = 10**(logLx_sol) * kx
    return np.log10(Lbol_sol)


def Lbol_to_Lx(logLbol_sol: float) -> float : 
    """Calculates log-scale X-ray luminosity in erg/s from log-scale bolometric luminosity in solar luminosities. 

    Args:
        logLbol_sol (float): log-scale bolometric luminosity in solar luminosity.

    Returns:
        float: log-scale X-ray luminosity in erg/s
    """
    kx = Lbol_to_kx(logLbol_sol)
    Lx_sol = 10**(logLbol_sol) / kx
    logLx_erg = LogSolsToLogErgs(np.log10(Lx_sol))
    return logLx_erg
    
    
def AddLbolAx(ax, isYax=True) :
    """Adds a secondary axis of bolomatric luminosity in solar luminosities for a log-scale X-ray luminosity in ergs/s.

    Args:
        ax: plot axis.
        isYax (bool, optional): Add to y axis when true, add to x when false. Defaults to True.

    Returns:
        _type_: secondary axis.
    """
    if(isYax) :
        secax = ax.secondary_yaxis('right', functions=(Lx_To_Lbol, Lbol_to_Lx))
        secax.set_ylabel('$\log( \; L_{BOL} \; [L_{\odot}] )$')
    else : 
        secax = ax.secondary_xaxis('right', functions=(Lx_To_Lbol, Lbol_to_Lx))
        secax.set_xlabel('$\log( \; L_{BOL} \; [L_{\odot}] )$')
    return secax