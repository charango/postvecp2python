import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
import re

from field import *

#   =============
#   MAIN FUNCTION
#   =============
def plottotaldissipation():

    # first import global variables
    import par

    globstr = Field(field=par.field, directory=par.directory[0]).globstr

    X = np.zeros(len(par.directory))  # forcing frequency
    Y = np.zeros(len(par.directory))  # dissipation
    
    # loop over directories
    for i in range(len(par.directory)):

        # read dissipation file, which contains (i) viscous
        # dissipation, (ii) power in pressure force, (iii) power in
        # differential rotation and (iv) forcing frequency
        with open(par.directory[i]+'/dissipation', 'r') as f:
            line = f.readline()
        res = line.split()

        X[i] = res[3]  # forcing frequency
        Y[i] = res[0]  # viscous dissipation
        # if you'd like to display instead the power exerted by the pressure
        # force, set:
        # Y = res[1]        
        # if you'd like to display instead the power exerted by the
        # differential rotation force term, set:
        # Y = res[2] 
        
    # -----------------------
    # figure and axes properties
    # -----------------------
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.16, right=0.98, top=0.96, bottom=0.14)
    ax = fig.gca()

    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
    ax.tick_params(axis='x', which='minor', top=True)
    ax.tick_params(axis='y', which='minor', right=True)
    ax.set_yscale('log')
    ax.set_xlabel(r'Forcing frequency $\omega_p$ in inertial frame')
    ax.set_ylabel(r'Dissipation')
    ax.set_xlim(X.min(),X.max())
    ax.grid(axis='both', which='major', ls='-', alpha=0.8)

    # modes with negative damping rates are shown by a blue filled
    # circle, those with positive damping rates are shown by a red
    # filled circle
    ax.scatter(X,np.abs(Y),s=10,marker='o',color='tab:blue')

    # display global string with main parameters in the bottom
    # use of set_title along with negative pad allows string
    # to be automatically centred in x-position
    ax.set_title(globstr,y=0, pad=-75,fontsize=16,color='black')
        
    # ------------------
    # save in pdf or png files
    # ------------------
    outfile = 'total_dissipation'+'_'+par.directory[0]+'_'+par.directory[len(par.directory)-1]
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=80)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=80)
    plt.close(fig)
