import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
import re

from field import *

def pltcolor(lst):
    cols=[]
    for l in lst:
        if l < 0.0:
            cols.append('tab:blue')
        else:
            cols.append('tab:red')
    return cols

#   =============
#   MAIN FUNCTION
#   =============
def plotqzeigenfq():

    # first import global variables
    import par

    # loop over directories
    for i in range(len(par.directory)):

        # read data
        myfield = Field(field=par.field, directory=par.directory[i])
        omr = myfield.omr
        omi = myfield.omi
        globstr = myfield.globstr

        # -----------------------
        # figure and axes properties
        # -----------------------
        fig = plt.figure(figsize=(8.,8.))
        plt.subplots_adjust(left=0.16, right=0.98, top=0.98, bottom=0.14)
        ax = fig.gca()

        ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
        ax.tick_params(axis='x', which='minor', top=True)
        ax.tick_params(axis='y', which='minor', right=True)
        ax.set_yscale('log')
        ax.set_xlabel(r'Eigenfrequency $\omega$')
        ax.set_ylabel(r'Damping rate |$\tau$|')
        ax.set_xlim(omi.min(),omi.max())
        ax.grid(axis='both', which='major', ls='-', alpha=0.8)

        # modes with negative damping rates are shown by a blue filled
        # circle, those with positive damping rates are shown by a red
        # filled circle
        cols=pltcolor(omr)
        ax.scatter(omi,np.abs(omr),s=10,marker='o',c=cols)

        # display global string with main parameters in the bottom
        # use of set_title along with negative pad allows string
        # to be automatically centred in x-position
        ax.set_title(globstr,y=0, pad=-75,fontsize=16,color='black')
        
        # ------------------
        # save in pdf or png files
        # ------------------
        outfile = 'qz_eigenfq'+'_'+par.directory[i]
        fileout = outfile+'.pdf'
        if par.saveaspdf == 'Yes':
            plt.savefig('./'+fileout, dpi=80)
        if par.saveaspng == 'Yes':
            plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=80)
        plt.close(fig)
