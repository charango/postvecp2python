import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
import sys
import re

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, LogFormatter)

from field import *

#   =============
#   MAIN FUNCTION
#   =============
def plotzcut():

    # first import global variables
    import par
    colored_cmap=par.mycolormap
    
    #print('par.directory =  ', par.directory)
    
    # loop over directories
    for i in range(len(par.directory)):

        if par.movie == 'Yes':
            print('directory ',str(i),' over ',len(par.directory), ' (',par.directory[i],')')
            
       # read data
        myfield = Field(field=par.field, directory=par.directory[i])
        X = myfield.x
        Y = myfield.y
        omr = myfield.omr
        omi = myfield.omi
        globstr = myfield.globstr
        
        # loop over modes
        for k in range(myfield.nmodes):

            # -----------------------
            # get data to be displayed and array name
            # -----------------------
            if par.field == 'ek':
                array = myfield.data[:,:,0,k]
                strfield = 'Kinetic energy'
            if par.field == 'dissv':
                array = myfield.data[:,:,1,k]
                strfield = 'Viscous dissipation'
            if par.field == 'shear':
                array = myfield.data[:,:,2,k]
                strfield = 'Shear power'
            if par.field == 'et':
                array = myfield.data[:,:,0,k]
                strfield = 'Thermal energy'
            if par.field == 'disst':
                array = myfield.data[:,:,1,k]
                strfield = 'Thermal dissipation'

            # Option to multiply field by distance to rotation axis,
            # to enhance contrast:
            if par.multbyaxisdist == 'Yes':
                for ii in range(len(X)):
                    for jj in range(len(X[ii])):
                        array[ii][jj]=array[ii][jj]*X[ii][jj]
                strfield += r' $\times$ s'
                
            # -----------------------
            # work out min/max colorbar
            # -----------------------
            if par.normalizetomax == 'Yes':
                array /= array.max()
                myfieldmin = 1e-6
                myfieldmax = 1.0
                strfield += ' (normalized to maximum)'
            else:
                myfieldmax = array.max()
                myfieldmin = 1e-6 * myfieldmax
                strfield += ' (code units)'
                
            if par.fieldmin != '#':
                myfieldmin = par.fieldmin
            if par.fieldmax != '#':
                myfieldmax = par.fieldmax
                
            mynorm = matplotlib.colors.LogNorm(vmin=myfieldmin,vmax=myfieldmax)

            # -----------------------
            # figure and axes properties
            # -----------------------
            fig = plt.figure(figsize=(8.,8.))
            plt.subplots_adjust(left=0.17, right=0.92, top=0.90, bottom=0.12)
            ax = plt.gca()
            ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
            ax.tick_params(axis='x', which='minor', top=True)
            ax.tick_params(axis='y', which='minor', right=True)
            ax.xaxis.set_minor_locator(MultipleLocator(0.1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.1))
            ax.set_xlabel('s')
            ax.set_ylabel('z')
            ax.set_xlim(0.0,1.0)
            ax.set_ylim(0.0,1.0)

            # -----------------------
            # display contour field
            # -----------------------        
            CF = ax.pcolormesh(X,Y,array,cmap=colored_cmap,norm=mynorm,rasterized=True)

            # ------------------
            # overplot characteristics
            # ------------------
            if par.plot_caract == 'Yes':
                nb_iterations = 5e4 # 5e4
                if myfield.gamma == 0.0:
                    sc,zc = myfield.compute_characteristics(niterations=nb_iterations,omega=omi[k])
                else:
                    sc,zc = myfield.compute_characteristics(niterations=nb_iterations,omega=myfield.gamma)
                ax.scatter(sc,zc,s=3,marker='.',alpha=0.5,color='white')

            # ------------------
            # overplot critical latitudes
            # ------------------
            if par.plot_critical_latitudes == 'Yes' and (myfield.rotation == 'solid' or myfield.rotation == 'shellular'):
                # only for solid-body rotation or shellular
                # differential rotation, as there is no simple
                # analytical expression for the critical latitudes
                # otherwise
                from matplotlib.lines import Line2D
                if myfield.gamma == 0.0:
                    sinthetain,sinthetaout = myfield.compute_critical_latitudes(omega=omi[k])
                else:
                    sinthetain,sinthetaout = myfield.compute_critical_latitudes(omega=myfield.gamma)
                xct = np.asarray([myfield.eta-0.02,myfield.eta+0.02])*np.sqrt(1.0-sinthetain**2.0)
                yct = np.asarray([myfield.eta-0.02,myfield.eta+0.02])*sinthetain
                line = Line2D(xct, yct, lw=2., color=par.c20[1], alpha=0.7)
                ax.add_line(line)
                xct = np.asarray([1.0-0.02,1.0+0.02])*np.sqrt(1.0-sinthetaout**2.0)
                yct = np.asarray([1.0-0.02,1.0+0.02])*sinthetaout
                line = Line2D(xct, yct, lw=2., color=par.c20[1], alpha=0.7)
                ax.add_line(line)

            # ------------------
            # overplot turning surfaces and/or critical layers, if any
            # ------------------
            if par.plot_turning_surfaces == 'Yes' or par.plot_critical_layers == 'Yes':
                if myfield.gamma == 0.0:
                    (xi,buf,buf,omegatilde) = myfield.compute_dzds_dsdz_caract(omi[k],myfield.n,1.0,X,Y)
                else:
                    (xi,buf,buf,omegatilde) = myfield.compute_dzds_dsdz_caract(myfield.gamma,myfield.n,1.0,X,Y)
                if par.plot_turning_surfaces == 'Yes':
                    plt.contour(X,Y,xi,levels=[0],colors='tab:brown',alpha=1.0,linewidths=2)
                if par.plot_critical_layers == 'Yes':
                    # case there is no differential rotation, in which
                    # case omegatilde is not an array but a single float:
                    if (not omegatilde == 'False'):
                        omegatilde = omegatilde*np.ones(xi.shape[0]*xi.shape[1]).reshape(xi.shape[0],xi.shape[1])
                    plt.contour(X,Y,omegatilde,levels=[0],colors='tab:red',alpha=1.0,linewidths=2)
                 
            # ---------------
            # display strings
            # ---------------
            # display real and imaginary parts of mode's
            # eigenfrequencies in top-right corner except if tidal
            # forcing is applied, in which case we simply display the
            # forcing frequency
            if myfield.gamma == 0.0:
                if omi[k] != 0.0:
                    strfq  = r'$|\omega|=$'+format(omi[k],'.4f')
                    ax.text(0.99,0.99,strfq,fontsize=16,color='black',horizontalalignment='right',verticalalignment='top')
                if omr[k] != 0.0:
                    strtau = r'$\tau=$'+par.str_fmt(omr[k])
                    ax.text(0.99,0.94,strtau,fontsize=16,color='black',horizontalalignment='right',verticalalignment='top')
            else:
                strfq  = r'$|\gamma|=$'+format(myfield.gamma,'.4f')
                ax.text(0.99,0.99,strfq,fontsize=16,color='black',horizontalalignment='right',verticalalignment='top')
            # display global string with main parameters in the bottom
            # use of set_title along with negative pad allows string
            # to be automatically centred in x-position
            ax.set_title(globstr,y=0, pad=-65,fontsize=16,color='black')
        
            # ----------------
            # plot color-bars
            # ----------------
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("top", size="2.5%", pad=0.12)
            cb = plt.colorbar(CF, cax=cax, orientation='horizontal')
            cax.xaxis.tick_top()
            cax.xaxis.set_tick_params(direction='out')
            cax.xaxis.set_label_position('top')
            cax.set_xlabel(strfield)
            cax.xaxis.labelpad = 8
            cax.xaxis.set_major_locator(ticker.LogLocator(base=10.0,numticks=8))
                
            # ------------------
            # save in pdf or png files
            # ------------------
            outfile = par.field+'_mode'+str(k)+'_'+par.directory[i]
            if par.movie == 'Yes':
                outfile = par.field+str(i).zfill(4)
            fileout = outfile+'.pdf'
            if par.saveaspdf == 'Yes':
                plt.savefig('./'+fileout, dpi=80)
            if par.saveaspng == 'Yes':
                plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=100)
            plt.close(fig)

            if par.onemode == 'Yes':
                break

    # ------------------
    # finally concatenate png if movie requested
    # ------------------
    if par.movie == 'Yes':
        # png files that have been created above
        allpngfiles = [par.field+str(x).zfill(4)+'.png' for x in range(len(par.directory))]
        str_on_start_number = str(0)
        # input files for ffpmeg
        input_files = par.field+'%04d.png'
        # output file for ffmpeg
        filempg = par.field+'_'+str(par.directory[0])+'_'+str(par.directory[len(par.directory)-1])+'.mpg'
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