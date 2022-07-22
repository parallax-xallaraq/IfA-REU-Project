# all imports
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import csv 

# color constants
c_ir        = '#D00000' # use this for IR sources                   ALT: '#CE2727' # Fire Engine Red 
c_xray      = '#0461A7' # use this for Xray sources                 ALT: '#384E77' # Y In Mn Blue    
c_ir_xray   = '#D1F520' # use this for matched IR and Xray sources  ALT: '#F6CD13' # Jonquil         
c_xray_pur  = '#8a236b' # use this for plots relating to X-ray colorbar

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

# redefine the defaults for plots
def SetStyle() : 
    mpl.rcParams['figure.figsize']  = 6, 4  # figure size in inches
    mpl.rcParams['axes.titlesize']  = 16    # size of title
    mpl.rcParams['axes.titlepad']   = 16    # space between title and plot 
    mpl.rcParams['axes.labelsize']  = 14    # font size of the x and y labels
    mpl.rcParams['axes.labelpad']   = 10    # space between label and axis
    mpl.rcParams['lines.linewidth'] = 1     # line width in points
    mpl.rcParams['figure.dpi']      = 150   # dots-per-inch

# Draw the selection wedge for Donley 2012 on a plot 
def PlotDonleyWedge(linewidth=1) : 
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
    plt.plot( [x_min,       x_min],         [y_min,     y_int_xmin], 'k', linewidth=linewidth)    # x >= 0.08
    plt.plot( [x_min,       x_int_ymin],    [y_min,     y_min],      'k', linewidth=linewidth)    # y >= 0.15
    plt.plot( [x_int_ymin,  max],           [y1_low,    y2_low],     'k', linewidth=linewidth)    # y >= 1.21x - 0.27
    plt.plot( [x_min,       max],           [y1_high,   y2_high],    'k', linewidth=linewidth)    # y <= 1.21x + 0.27

    # name the x and y axis 
    plt.xlabel('$\log(f_{5.8um}/f_{3.6um})$') 
    plt.ylabel('$\log(f_{8.0um}/f_{4.5um})$')

# returns a boolean mask that is True for sources within the Donley 2012 wedge 
def SourcesInDonleyWedge(x,y) :
    # initialize mask with all False 
    mask_inWedge = np.zeros(x.size, dtype=bool)
    # check all datapoints
    for i in range(len(x)) :
        # inside wedge 
        if(     x[i]>=0.08 and 
                y[i]>=0.15 and 
                y[i]>=(1.21*x[i]-0.27) and 
                y[i]<=(1.21*x[i]+0.27)
            ) :  
            # set index to true
            mask_inWedge[i] = True
    # return mask
    return(mask_inWedge)

# add 'n = #' to bottom right of plot
def addtext_n(n, pre='n = '):
    plt.text(   0.95,                           # x
                0.05,                           # y 
                pre + str(n),                # string
                transform=plt.gca().transAxes,  # use axis coordinants
                horizontalalignment='right'     # alignment 
    )

# save plot 
def save(filename) :
    plt.savefig(    filename,
                    bbox_inches ="tight",
                    pad_inches=0.2,
                    facecolor='w',
                    edgecolor='w'
    )

# plot 3 types of AGN on Donley IRAC color-color plot
def PlotDonleyTricolor(
        x_xr, y_xr,             # X-ray (blue)
        x_ir, y_ir,             # infrared (red)
        x_ix, y_ix,             # overlap (yellow)
        path='', fileName='',   # save
        saveNum=2,              # 3=full,zoom,zoom w/legend. 2= zoom,zoom w/legend. 1=zoom w/legend. 0=none
        printInfo=True          # output text
    ) :

    # plot data
    plt.scatter(x_xr, y_xr, marker='s', c=c_xray,     s=3, label='X-ray luminous (n='+str(len(x_xr))+')')
    plt.scatter(x_ir, y_ir, marker='^', c=c_ir,       s=3, label='Donley selected (n='+str(len(x_ir))+')')
    plt.scatter(x_ix, y_ix, marker='o', c=c_ir_xray,  s=3, label='Donley & X-ray (n='+str(len(x_ix))+')')

    # plotting class
    PlotDonleyWedge()
    addtext_n(len(x_ir)+len(x_xr)+len(x_ix))

    # make plot square
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_adjustable('box') 

    # save
    if(path and fileName) :
        if(saveNum>=3) : 
            plt.axis([-1.5,2.5,-1.5,2.5]) 
            save(path+'\\'+fileName+'_FULL.png')
        if(saveNum>=2) : 
            plt.axis([-0.7,1.0,-0.7,1.0]) 
            save(path+'\\'+fileName+'_ZOOM.png')
        if(saveNum>=1) :
            plt.axis([-0.7,1.0,-0.7,1.0]) 
            plt.legend(markerscale=3)
            save(path+'\\'+fileName+'_ZOOM_legend.png')
        if(saveNum<1) :
            plt.axis([-0.7,1.0,-0.7,1.0]) 
            plt.legend(markerscale=3)
    else :
        plt.axis([-0.7,1.0,-0.7,1.0]) 
        plt.legend(markerscale=3) 

    # display
    plt.show()

    # print info 
    if(printInfo) :
        print('Number of IR only:\t',   len(x_ir))
        print('Number of Xray only:\t', len(x_xr))
        print('Number of matches:\t',   len(x_ix))
        print('Total Number:\t\t',      len(x_ir)+len(x_xr)+len(x_ix))
        print( len(x_xr) / (len(x_xr)+len(x_ix)) * 100., '% of Xray sources fall outside wedge.')
        print( len(x_ir) / (len(x_ir)+len(x_ix)) * 100. ,'% of IR souces have no Xray greater than 10^43')
    
    # done
    return

# plot AGN w/ and w/o Xrays on Donley IRAC color-color plot with X-ray luminosity colorbar
def PlotDonleyXray(
        x_nX, y_nX,             # no X-ray detections 
        x_yX, y_yX,             # yes X-ray detections
        Lx,                     # X-ray luminosity (colorbar)
        path='', fileName='',   # save
        saveAll=False,          # save 3 versions of plot 
        printInfo=True          # output text
    ) :
    
    # plot data16
    plt.scatter(x_nX, y_nX, marker='x', s=4, c='lightgrey', label='No X-ray (n='+str(len(x_nX))+')')
    plt.scatter(x_yX, y_yX, marker='D', s=3, c=Lx, cmap=plt.cm.get_cmap('turbo',9), label='Has X-ray (n='+str(len(x_yX))+')')

    # color bar 
    plt.clim(41.5, 46.0) # colorbar limits 
    plt.colorbar(label='$\log( \; L_{x(0.5-10keV)} \; [erg \; s^{-1}] \;)$')

    # plotting
    PlotDonleyWedge()
    addtext_n(len(x_nX)+len(x_yX))

    # make plot square
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_adjustable('box') 

    # save
    if(path and fileName) :
        if(saveAll) : 
            plt.axis([-1.5,2.5,-1.5,2.5]) 
            save(path+'\\'+fileName+'_FULL.png')
            plt.axis([-0.7,1.0,-0.7,1.0]) 
            save(path+'\\'+fileName+'_ZOOM.png')    
            plt.legend(markerscale=3) 
            save(path+'\\'+fileName+'_ZOOM_legend.png')
        else : 
            plt.axis([-0.7,1.0,-0.7,1.0]) 
            plt.legend(markerscale=3) 
            save(path+'\\'+fileName+'_ZOOM_legend.png')

    else :
        plt.axis([-0.7,1.0,-0.7,1.0]) 
        plt.legend(markerscale=3) 

    # display
    plt.show()

    #print info 
    if(printInfo) : 
        print('Number with Xray:\t',    len(x_yX))
        print('Number without Xray:\t', len(x_nX))
        print('Total Number:\t\t',      len(x_yX)+len(x_nX))
    
    # done    
    return

# plot AGN w/ Xrays on Donley IRAC color-color plot with X-ray luminosity colorbar
def PlotDonleyXray_Xonly(
        x_yX, y_yX,             # yes X-ray detections
        Lx,                     # X-ray luminosity (colorbar)
        path='', fileName='',   # save
        saveAll=False,          # save 3 versions of plot 
    ) :
    
    # plot data16
    plt.scatter(x_yX, y_yX, marker='D', s=3, c=Lx, cmap=plt.cm.get_cmap('turbo',9))

    # color bar 
    plt.clim(41.5, 46.0) # colorbar limits 
    plt.colorbar(label='$\log( \; L_{x(0.5-10keV)} \; [erg \; s^{-1}] \;)$')

    # plotting
    PlotDonleyWedge()
    addtext_n(len(x_yX))

    # make plot square
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_adjustable('box') 

    # save
    if(path and fileName) :
        if(saveAll) : 
            plt.axis([-1.5,2.5,-1.5,2.5]) 
            save(path+'\\'+fileName+'_FULL.png')
            plt.axis([-0.7,1.0,-0.7,1.0]) 
            save(path+'\\'+fileName+'_ZOOM.png')    
        else : 
            plt.axis([-0.7,1.0,-0.7,1.0]) 
            save(path+'\\'+fileName+'_ZOOM.png')
    else :
        plt.axis([-0.7,1.0,-0.7,1.0]) 

    # display
    plt.show()

# Plots one X-ray luminosty (log scale) histogram 
def PlotHistOne(x,saveStr=''):
    # Plot histogram.
    plt.hist(x, 25, edgecolor='w', color=c_xray_pur)
    # mean
    mean = np.array(x).mean()
    min_ylim, max_ylim = plt.ylim()
    plt.axvline(mean, color='k', linestyle='dashed')
    plt.text(mean*1.001, max_ylim*0.94, 'Mean: {:.2f}'.format(mean))
    # set axis lables
    plt.xlabel('$\log( \; L_{x(0.5-10keV)} \; [erg \; s^{-1}] \;)$')
    plt.ylabel('Number')
    # save
    if(saveStr) :
        save(saveStr)
    # display
    plt.show()

# Plots two X-ray luminosty (log scale) histograms
def PlotHistTwo(x1,x2,h=300,saveStr='',c=c_xray) : 
    # subplots 
    fig, (x1hist, x2hist) = plt.subplots(1,2, sharey=True)
    ## inWedge subplot
    # plot all redshift histogram
    x1hist.hist(x1, bins=np.arange(42,46,0.25), edgecolor='w', color=c)
    x1hist.set_ylim(ymin=0, ymax=h)
    # axis and titles 
    x1hist.set_xlabel('$\log( \; L_{x(0.5-10keV)} \; [erg \; s^{-1}] \;)$')
    x1hist.set_ylabel('Number')
    x1hist.set_xticks([42,43,44,45,46])
    # mean 
    mean_x1 = x1.mean()
    min_ylim_all, max_ylim_all = x1hist.get_ylim()
    x1hist.axvline(mean_x1, color='k',linestyle='dashed')
    x1hist.text(mean_x1*1.001, max_ylim_all*0.92, 'Mean: {:.2f}'.format(mean_x1))
    ## outWedge subplot
    # plot agn redshift histogram
    x2hist.hist(x2, bins=np.arange(42,46,0.25), edgecolor='w', color=c)
    # axis and titles 
    x2hist.set_xlabel('$\log( \; L_{x(0.5-10keV)} \; [erg \; s^{-1}] \;)$')
    x2hist.set_ylim(ymin=0, ymax=h)
    x2hist.set_xticks([42,43,44,45,46])
    # mean 
    mean_x2 = x2.mean()
    min_ylim_agn, max_ylim_agn = x2hist.get_ylim()
    x2hist.axvline(mean_x2, color='k',linestyle='dashed')
    x2hist.text(mean_x2*1.001, max_ylim_agn*0.92, 'Mean: {:.2f}'.format(mean_x2))
    ## end subplots 
    # formatting and save  
    plt.tight_layout()
    # save
    if(saveStr) :
        save(saveStr)
    # print info 
    print('Number of x1:\t', len(x1))
    print('Number of x2:\t', len(x2))
    print('Range x1:\t',  min(x1), '-', max(x1))
    print('Range x2:\t',  min(x2), '-', max(x2))