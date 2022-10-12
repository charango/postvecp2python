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

        globstr += r' $\gamma=$'+format(myfield.gamma,'.4f')
        

        # get number of domains
        ndom = len(myfield.ck[0,:,0,0])

        # get number of variables
        nvar = len(myfield.ck[0,0,:,0])

        # different linestyles for different domains (dimension should be enough!...)
        linestyle_dom = ['solid','dotted','dashed','dashdot','loosely dotted','loosely dashed','loosely dashdotted']
        
        # loop over modes
        for k in range(myfield.nmodes):

            # -----------------------
            # get data to be displayed and array name
            # -----------------------

            # first: one domain assumed -> second index set to 0
            ck_var0 = myfield.ck[:,0,0,k]
            ck_var1 = myfield.ck[:,0,1,k]
            cl_var0 = myfield.cl[:,0,0,k]
            cl_var1 = myfield.cl[:,0,1,k]

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
            ax[0].grid(axis='both', which='major', ls='-', alpha=0.8)
            
            # loop over variables
            for v in range(nvar):
                # loop over domains
                for d in range(ndom):
                    ck_var = myfield.ck[:,d,v,k]
                    X = np.arange(len(ck_var[ck_var!=0.0]))
                    if v == 0 and d == 0:
                        ax[0].set_xlim(X.min(),X.max())
                    if d == 0:
                        ax[0].plot(X, ck_var[ck_var!=0.0], color=par.c20[3*v], lw=2., linestyle=linestyle_dom[d], label=myfield.variables[v])
                    else:
                        ax[0].plot(X, ck_var[ck_var!=0.0], color=par.c20[3*v], lw=2., linestyle=linestyle_dom[d])
                        
            # legend put only once
            legend = fig.legend(loc='lower left',fontsize=16,facecolor='white',edgecolor='white',framealpha=0.7,numpoints=1,bbox_to_anchor=(0.001,0.001),bbox_transform=plt.gcf().transFigure)
            for line, text in zip(legend.get_lines(), legend.get_texts()):
                text.set_color(line.get_color())
                        
            # right plot: cl spectrum
            ax[1].tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
            ax[1].tick_params(axis='x', which='minor', top=True)
            ax[1].tick_params(axis='y', which='minor', right=True)
            ax[1].set_yscale('log')
            ax[1].set_xlabel('l')
            ax[1].set_ylabel(r'$C_l$')
            ax[1].grid(axis='both', which='major', ls='-', alpha=0.8)
            index0 = myfield.lmin[0]+2*np.arange(1+(myfield.lmax-myfield.lmin[0])/2.0)
            index1 = myfield.lmin[1]+2*np.arange(1+(myfield.lmax-myfield.lmin[1])/2.0)

            # loop over variables
            for v in range(nvar):
                # loop over domains
                for d in range(ndom):
                    cl_var = myfield.cl[:,d,v,k]
                    X = np.arange(len(cl_var))
                    if v == 0 and d == 0:
                        ax[1].set_xlim(X.min(),X.max())
                    ax[1].plot(index1, cl_var[0:len(index1)], color=par.c20[3*v], lw=2., linestyle=linestyle_dom[d], label=myfield.variables[v])
            
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


    # ------------------
    # finally concatenate png if movie requested
    # ------------------
    if par.movie == 'Yes':
        # png files that have been created above
        allpngfiles = ['spectrum_mode0_out'+str(x).zfill(3)+'.png' for x in range(len(par.directory))]
        str_on_start_number = str(0)
        # input files for ffpmeg
        input_files = 'spectrum_mode0_out%03d.png'
        # output file for ffmpeg
        filempg = 'spectrum_mode0_'+str(par.directory[0])+'_'+str(par.directory[len(par.directory)-1])+'.mpg'
        # call to python-ffmpeg
        import ffmpeg
        (
            ffmpeg            
            .input(input_files, framerate=10, start_number=str_on_start_number)
            # framerate=10 means the video will play at 10 of the original images per second
            .output(filempg, r=30, pix_fmt='yuv420p', **{'qscale:v': 3})
            # r=30 means the video will play at 30 frames per second
            .overwrite_output()
            .run()
        )
        # erase png files
        allfiles = ' '.join(allpngfiles)
        os.system('rm -f '+allfiles)
