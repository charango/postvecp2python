import numpy as np
import sys
import os
import re

# ---------------
# Plotting parameters
# ---------------
import matplotlib
matplotlib.rcParams.update({'font.size': 20})
#matplotlib.rc('font', family='Arial')
matplotlib.rcParams['font.family'] = 'DeJavu Serif'
matplotlib.rcParams['font.serif'] = ['Helvetica']
fontcolor='white'

# ---------------
# special color scale for 1D plots:
# ---------------
c20 = [(31, 119, 180), (255, 127, 14), (174, 199, 232), (255, 187, 120),    
       (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
       (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
       (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
       (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
for i in range(len(c20)):
    r, g, b = c20[i]    
    c20[i] = (r / 255., g / 255., b / 255.)

    
# ---------------
# read parameter file
# ---------------
params = open("paramsp2p.dat",'r')
lines_params = params.readlines()
params.close()                 

par = []                       # allocating a dictionary
var = []                       # allocating a dictionary
regex = re.compile(',')

for line in lines_params:      
    try:
        line.split()[0][0]     # check if blank line (GWF)
    except:
        continue
    if (line.split()[0][0]=='#'): # check if line starts with a # (comment)
        continue
    else:
        name, value = line.split()[0:2]
    try:
        float(value)           # first trying with float
    except ValueError:         # if it is not float
        try:
            int(value)         # we try with integer
        except ValueError:     # if it is not integer nor float
            if(regex.search(value) != None):  # we try with array with several values separated by a comma (,)
                if name == 'caract_s_z':
                    value = [float(x) for x in value.split(',')]
                if name == 'directory':
                    value = [str(x) for x in value.split(',')]
            else:
                value = '"' + value + '"'   # if none of the above tests work, we know value is a string
    par.append(value)
    var.append(name)

for i in range(len(par)):
    exec(var[i]+"="+str(par[i]))

    
# ---------------
# below we work out specific parameters and/or error messages
# ---------------
if directory == 'all':
    from glob import glob
    directory = sorted(glob("out*"))
    
# read edp.choices and check that the right options have been ticked:
if isinstance(directory, str) == True:
    directory = [directory]
if isinstance(directory, (int, float)) == True:
    directory = [str(directory)]

for i in range(len(directory)):
    with open(directory[i]+'/edp.choices', 'r') as f:
        f.readline()
        f.readline()
        f.readline()
        options = f.readline().split()
        opt_ekdiss = int(options[0])
        opt_spectr = int(options[2])
        opt_temp   = int(options[8])

    if ( (plot_zcut == 'Yes') and (field == 'ek' or field == 'dissv' or field == 'shear') and opt_ekdiss == 0):
        sys.exit("A meridional cut of the modes kinetic energy or viscous dissipation is requested, but such quantities have not been computed according to the options set in the edp.choice file. Please edit the edp.choices file accordingly (first integer at fourth line should be set to 1!) and run postvecp again.")

    if ( (plot_zcut == 'Yes') and (field == 'et' or field == 'disst') and opt_temp == 0):
        sys.exit("A meridional cut of the modes thermal energy or thermal dissipation is requested, but such quantities have not been computed according to the options set in the edp.choice file. Please edit the edp.choices file accordingly (ninth integer at fourth line should be set to 1!) and run postvecp again.")
        
    if ( plot_spectrum == 'Yes' and opt_spectr == 0):
        sys.exit("The modes spectral content is requested to be plotted but the spectra have not been computed according to the options set in the edp.choices file. Please edit the edp.choices fie accordingly (third integer at fourth line should be set to 1!) and run postvecp again.")

        
# Case of QZ calculation
if plot_qz_eigenfq == 'Yes':
    if plot_zcut == 'Yes':
        plot_zcut = 'No'
        print('When plot_qz_eigenfq, you cannot display z-cuts since no specific eigenmodes are computed for, i am therefore setting plot_zcut to No.')
    if plot_spectrum == 'Yes':
        plot_spectrum = 'No'
        print('When plot_qz_eigenfq, you cannot display the spectral content of a mode since no specific eigenmodes are computed for, i am therefore setting plot_spectrum to No.')


# Case total dissipation plotted against frequency
if plot_total_dissipation == 'Yes' and len(directory) == 1:
    sys.exit('A plot with total dissipation for several directories has been requested and yet only one directory has been specified. Please modify your paramsp2p.dat file by setting eg directory to "all".')
    
    
if not('elevationfactor' in open('paramsp2p.dat').read()): 
    elevationfactor = 1         # used by plot_zcut_3D
if not('convert_density_factor' in open('paramsp2p.dat').read()):
    convert_density_factor = 72 # used by plot_zcut_3D
if not('auto_huefactor' in open('paramsp2p.dat').read()):
    auto_huefactor = 1.         # used by plot_zcut and plot_zcut_3D
if not('pv_specular' in open('paramsp2p.dat').read()):
    pv_specular = 0.            # used by plot_zcut (3D)
if not('pv_light_intensity' in open('paramsp2p.dat').read()):
    pv_light_intensity = 0.     # used by plot_zcut (3D). If 0 we use a function of elevationfactor

# case an animation is requested
if movie == 'Yes':
    onemode = 'Yes'
    saveaspdf = 'No'
    saveaspng = 'Yes'
    if isinstance(directory, str) == True:
        sys.exit("ERROR: you requested an animation but specified only one directory...")

# color map
if not('mycolormap' in open('paramsp2p.dat').read()):
    mycolormap = 'nipy_spectral'

    
# ---------------
# Ancillary functions
# ---------------
    
def round_to_n(x, n):
    " Round x to n significant figures "
    return round(x, -int(np.floor(np.sign(x) * np.log10(abs(x)))) + n)

def str_fmt(x, n=3):
    " Format x into nice Latex rounding to n"
    if x > 0:
        power = np.floor(np.log10(round_to_n(x, 0)))
        f_SF = round_to_n(x, n) * pow(10, -power)
    else:
        power = np.floor(np.log10(round_to_n(-x, 0)))
        f_SF = round_to_n(-x, n) * pow(10, -power)
    if f_SF != 1.0:
        mystr = "$"+format(f_SF,'.2f')+r"\times 10^{"+format(power,'.0f')+"}$"
        if x < 0:
            mystr = "$"+format(-f_SF,'.2f')+r"\times 10^{"+format(power,'.0f')+"}$"
    else:
        mystr = r"$10^{"+format(power,'.0f')+"}$"
        if x < 0:
            mystr = -r"$10^{"+format(power,'.0f')+"}$"
    return mystr

def str_fmt_ek(x, n=3):
    " Format x into nice Latex rounding to n"
    power = np.floor(np.log10(round_to_n(x, 0)))
    f_SF = round_to_n(x, n) * pow(10, -power)
    if f_SF != 1.0:
        mystr = "$"+format(f_SF,'.1f')+r"\times 10^{"+format(power,'.0f')+"}$"
    else:
        mystr = r"$10^{"+format(power,'.0f')+"}$"
    return mystr
