import matplotlib.pyplot as plt
import matplotlib as mpl


# color constants
c_ir        = '#CE2727'     # Fire Engine Red   # use this for IR sources
c_xray      = '#384E77'     # Y In Mn Blue      # use this for Xray sources
c_ir_xray   = '#F9DC5C'     # Naples Yellow     # use this for matched IR and Xray sources


# redefine the defaults for plots
def SetStyle() : 
    mpl.rcParams['figure.figsize']  = 6, 4  # figure size in inches
    mpl.rcParams['axes.titlesize']  = 16    # size of title
    mpl.rcParams['axes.titlepad']   = 16    # space between title and plot 
    mpl.rcParams['axes.labelsize']  = 14    # font size of the x and y labels
    mpl.rcParams['axes.labelpad']   = 10    # space between label and axis
    mpl.rcParams['lines.linewidth'] = 0.5   # line width in points


# Draw the selection wedge for Donley 2012 on a plot 
def PlotDonleyWedge() : 
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
    plt.plot( [x_min,       x_min],         [y_min,     y_int_xmin], 'k' )    # x >= 0.08
    plt.plot( [x_min,       x_int_ymin],    [y_min,     y_min],      'k' )    # y >= 0.15
    plt.plot( [x_int_ymin,  max],           [y1_low,    y2_low],     'k' )    # y >= 1.21x - 0.27
    plt.plot( [x_min,       max],           [y1_high,   y2_high],    'k' )    # y <= 1.21x + 0.27

    # name the x and y axis 
    plt.xlabel('$\log(f_{5.8um}/f_{3.6um})$') 
    plt.ylabel('$\log(f_{8.0um}/f_{4.5um})$') 

