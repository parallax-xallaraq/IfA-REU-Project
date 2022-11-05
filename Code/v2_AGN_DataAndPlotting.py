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

# path variables 
path_cat = 'D:\\IfA REU\\Catalogs'
path_fig = 'D:\\IfA REU\\IfA-REU-Project\\Figures'
path_csv = 'D:\\IfA REU\\IfA-REU-Project\\Lists'

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

# add 'n = #' to bottom right of plot
def AddText_n(n, pre='n = '):
    plt.text(   0.95,                           # x
                0.05,                           # y 
                pre + str(n),                # string
                transform=plt.gca().transAxes,  # use axis coordinants
                horizontalalignment='right'     # alignment 
    )

# save plot 
def Save(filename) :
    plt.savefig(    filename,
                    bbox_inches ="tight",
                    pad_inches=0.2,
                    facecolor='w',
                    edgecolor='w'
    )