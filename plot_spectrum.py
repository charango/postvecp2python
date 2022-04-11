import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
import sys
import re

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, LogFormatter)
from matplotlib.lines import Line2D

from field import *

#   =============
#   MAIN FUNCTION
#   =============
def plotspectrum():

    # first import global variables
    import par
    colored_cmap=par.mycolormap

    # loop over directories
    for i in range(len(par.directory)):

        # read data
        myfield = Field(field=par.field, directory=par.directory[i])
        omr = myfield.omr
        omi = myfield.omi
        globstr = myfield.globstr
        
        # loop over modes
        for k in range(myfield.nmodes):

            # -----------------------
            # get data to be displayed and array name
            # -----------------------
            # one domain assumed -> second index set to 0
            ckw = myfield.ck[:,0,0,k]
            cku = myfield.ck[:,0,1,k]
            # means a priori that temperature was computed as a third variable in the equations
            if len(myfield.ck[0,0,:,0]) == 3:
                ckt = myfield.ck[:,0,2,k]
            clw = myfield.cl[:,0,0,k]
            clu = myfield.cl[:,0,1,k]
            # means a priori that temperature was computed as a third variable in the equations
            if len(myfield.cl[0,0,:,0]) == 3:
                clt = myfield.cl[:,0,2,k]
            
            # -----------------------
            # figure and axes properties
            # -----------------------
            fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(16., 8.))

            # display global string with main parameters in the bottom
            fig.suptitle(globstr,y=0.04,fontsize=16,color='black')
            
            # left plot: ck spectrum
            ax[0].tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
            ax[0].tick_params(axis='x', which='minor', top=True)
            ax[0].tick_params(axis='y', which='minor', right=True)
            ax[0].set_yscale('log')
            ax[0].set_xlabel('k')
            ax[0].set_ylabel(r'$C_k$')
            X = np.arange(len(cku))
            ax[0].set_xlim(X.min(),X.max())
            ax[0].grid(axis='both', which='major', ls='-', alpha=0.8)
            ax[0].plot(X, cku, color=par.c20[0], lw=2., linestyle = 'solid', label=r'$u_l$')
            ax[0].plot(X, ckw, color=par.c20[1], lw=2., linestyle = 'solid', label=r'$w_l$')
            if len(myfield.ck[0,0,:,0]) == 3:
                ax[0].plot(X, ckt, color=par.c20[2], lw=2., linestyle = 'solid', label=r'$t_l$')

            # legend put only once
            legend = fig.legend(loc='lower left',fontsize=16,facecolor='white',edgecolor='white',framealpha=1.0,numpoints=1,bbox_to_anchor=(0.02,0.02),bbox_transform=plt.gcf().transFigure)
            for line, text in zip(legend.get_lines(), legend.get_texts()):
                text.set_color(line.get_color())
                        
            # right plot: cl spectrum
            ax[1].tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
            ax[1].tick_params(axis='x', which='minor', top=True)
            ax[1].tick_params(axis='y', which='minor', right=True)
            ax[1].set_yscale('log')
            ax[1].set_xlabel('l')
            ax[1].set_ylabel(r'$C_l$')
            X = np.arange(len(clu))
            ax[1].set_xlim(X.min(),X.max())
            ax[1].grid(axis='both', which='major', ls='-', alpha=0.8)
            index0 = myfield.lmin[0]+2*np.arange(1+(myfield.lmax-myfield.lmin[0])/2.0)
            index1 = myfield.lmin[1]+2*np.arange(1+(myfield.lmax-myfield.lmin[1])/2.0)
            ax[1].plot(index0, clu[0:len(index0)], color=par.c20[0], lw=2., linestyle = 'solid', label=r'$u_l$')
            ax[1].plot(index1, clw[0:len(index1)], color=par.c20[1], lw=2., linestyle = 'solid', label=r'$w_l$')
            if len(myfield.cl[0,0,:,0]) == 3:
                ax[1].plot(index1, clt[0:len(index1)], color=par.c20[2], lw=2., linestyle = 'solid', label=r'$t_l$')
            
            plt.subplots_adjust(left=0.10, bottom=0.15, right=0.97, top=0.95, wspace=0.3)

            # ------------------
            # save in pdf or png files
            # ------------------
            outfile = 'spectrum_mode'+str(k)+'_'+par.directory[i]
            fileout = outfile+'.pdf'
            if par.saveaspdf == 'Yes':
                plt.savefig('./'+fileout, dpi=80)
            if par.saveaspng == 'Yes':
                plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=80)
            plt.close(fig)

            if par.onemode == 'Yes':
                break
