from field import *             # python library for postvecp
from pyevtk.hl import gridToVTK # https://pypi.org/project/pyevtk/
import colorsys                 # hsv2rgb
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import re

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, LogFormatter)

# to insert image:
from PIL import Image
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

convert_density = 72. # density for conversion between pdf and png.

#   =============
#   MAIN FUNCTION
#   =============
def plotzcut3D():

    # first import global variables
    import par
    
    # loop over directories
    for i in range(len(par.directory)):
    
    # read data
        myfield = Field(field=par.field, directory=par.directory[i])
        omr = myfield.omr
        omi = myfield.omi
        globstr = myfield.globstr

        for k in range(myfield.nmodes):

            if (par.directory[i] == '.'):
                fileout_prefix = par.field+'3D_mode'+str(k)
                pv_fileout_prefix = par.field+'_mode'+str(k)+'_pv'
                file_colormap = "colormap_"+par.field+"_mode"+str(k)+".xml"
            else:
                fileout_prefix = par.field+'3D_mode'+str(k)+'_'+par.directory[i]
                pv_fileout_prefix = par.field+'_mode'+str(k)+'_pv_'+par.directory[i]
                file_colormap = "colormap_"+par.field+"_mode"+str(k)+'_'+par.directory[i]+".xml"

            # -----------------------
            # Prepare image
            # -----------------------
            strfield = prepare_image(myfield, k, fileout_prefix, pv_fileout_prefix, file_colormap)

            # -----------------------
            # figure and axes properties
            # -----------------------
            fig = plt.figure(figsize=(8.,8.))
            plt.subplots_adjust(left=0.17, right=0.92, top=0.92, bottom=0.17)
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
            im = Image.open(pv_fileout_prefix+'.png')
            im = OffsetImage(im, zoom=0.5635*(1+par.elevationfactor/100)*72/convert_density)
            ab = AnnotationBbox(im, (0.5, 0.499), frameon=False)
            ab.set(zorder=-100.)
            ax.add_artist(ab)

            # ---------------
            # display strings
            # ---------------
            # display real and imaginary parts of mode's
            # eigenfrequencies in top-right corner except if tidal
            # forcing is applied, in which case we simply display the
            # forcing frequency
            if myfield.gamma == 0.0:
                strsuptitle=strfield
                if omi[k] != 0.0:
                    strfq  = r'$|\omega|=$'+format(omi[k],'.4f')
                    ax.text(0.99,0.82,strfq,fontsize=16,color='black',horizontalalignment='right',verticalalignment='top', clip_on=True) # it doesn't show, so we use suptitle
                if omr[k] != 0.0:
                    strtau = r'$\tau=$'+par.str_fmt(omr[k])
                    ax.text(0.99,0.87,strtau,fontsize=16,color='black',horizontalalignment='right',verticalalignment='top')
            else:
                strfq  = r'$|\gamma|=$'+format(myfield.gamma,'.4f')
                ax.text(0.99,0.90,strfq,fontsize=16,color='black',horizontalalignment='right',verticalalignment='top')
            # display global string with main parameters in the bottom
            # use of set_title along with negative pad allows string
            # to be automatically centred in x-position
            ax.set_title(globstr,y=0, pad=-65,fontsize=16,color='black')
            plt.suptitle(strsuptitle, y=0.98, x=0.5, fontsize=20)


            # ------------------
            # save in pdf or png files
            # ------------------
            outfile = par.field+'3D_mode'+str(k)+'_'+par.directory[i]
            if par.movie == 'Yes':
                outfile = par.field+str(i).zfill(4)
            fileout = fileout_prefix+'.pdf'
            if par.saveaspdf == 'Yes':
                plt.savefig('./'+fileout, dpi=80)
            if par.saveaspng == 'Yes':
                plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=100)
            plt.close(fig)

# cleaning:
            cmd = 'rm -f '+file_colormap+' '+pv_fileout_prefix+'*'
#           print(cmd)
#           os.system(cmd)
            
            if par.onemode == 'Yes':
                break

            
# -----------------------
def prepare_image(myfield, k, fileout_prefix, pv_fileout_prefix, file_colormap):
# -----------------------
    import par
    nmode = str(k)    # mode number. Start from 0
    VTK_name = pv_fileout_prefix
    X = myfield.x
    Y = myfield.y
    # Dimensions
    nr  = len(X)-1
    nth = len(X[0])-1
    nz = 1
    npoints = nr*nth
    nbins = max(int(np.sqrt(npoints)),1000)


    array = np.empty_like(myfield.data[:,:,0,k])
    if par.field == 'ek':
        strfield = 'Kinetic energy'
        if par.multbyaxisdist == 'Yes':
            for ii in range(len(X)):
                for jj in range(len(X[ii])):
                    array[ii][jj]=abs(myfield.data[ii, jj, 0, k]*X[ii][jj])
        else:
            for ii in range(len(X)):
                for jj in range(len(X[ii])):
                    array[ii][jj]=abs(myfield.data[ii, jj, 0, k])

    if par.field == 'dissv':
        strfield = '|Viscous dissipation|'
        if par.multbyaxisdist == 'Yes':
            for ii in range(len(X)):
                for jj in range(len(X[ii])):
                    array[ii][jj]=abs(myfield.data[ii, jj, 1, k]*X[ii][jj])
        else:
            for ii in range(len(X)):
                for jj in range(len(X[ii])):
                    array[ii][jj]=abs(myfield.data[ii, jj, 1, k])

    if par.multbyaxisdist == 'Yes':
        strfield = r'log$_{10}$('+strfield+r' $\times$ s)'
    else:
        strfield = 'log10('+strfield+')'
    
    if par.normalizetomax == 'Yes':
      # arraymax=array.max()
        array = np.log10(array/array.max()+1.e-11)
        strfield += ' (normalized to max.)'
    else:
        array = np.log10(array+1.e-32)
    myfieldmin = array.min()
    myfieldmax = array.max()
    print ('minmax:',myfieldmin,myfieldmax)

#   ncells = nr * nth * nz
    
    x  = np.zeros((nr, nth, nz))
    y  = np.zeros((nr, nth, nz))
    z  = np.zeros((nr, nth, nz))
    scalar = np.zeros((nr, nth, nz))
    
    for k in range(nz):
        for j in range(nth):
            for i in range(nr):
                x     [i, j, k] = X[i,j]
                y     [i, j, k] = Y[i,j]
                z     [i, j, k] = array[i,j]/(myfieldmax-myfieldmin) * 0.06 * par.elevationfactor
                scalar[i, j, k] = array[i,j]
    
#-------------------------------------------------------------------------------
# create color palette:

    bins = []
    scalarflat = scalar.reshape(-1)
    scalarsorted = np.sort(scalarflat)
    data_points_per_bin = len(scalarsorted) // nbins
    bins = [scalarsorted[_ * data_points_per_bin: (_+1)*data_points_per_bin] for _ in range(nbins)]

    hue  = np.zeros(nbins)
    sat  = np.zeros(nbins)
    val  = np.zeros(nbins)
    ndata = (nr)*(nth)

    for k in range(nbins):
        hue[k] = (float(k)/nbins)**3.0*0.14*par.auto_huefactor
        val[k] = (float(k)/nbins)**.5
        sat[k] = 1

    fh = open (file_colormap, "w")
    fh.write ('<ColorMaps>' + "\n" + '<ColorMap name="inertial" space="RGB">' + "\n")
    for ii in range(nbins):
       red, green, blue = colorsys.hsv_to_rgb(hue[ii],1,val[ii])
       scalarvalue=scalarsorted[int(ii*ndata/nbins)]
       fh.write ("<Point x=\"%-13.6f\" o=\"1\" r=\"%-13.6f\" g=\"%-13.6f\" b=\"%-13.6f\"/>\n" % (scalarvalue, red, green, blue))

    fh.write ('</ColorMap>' + "\n" + '</ColorMaps>')
    fh.close()
#-------------------------------------------------------------------------------

    gridToVTK(
        VTK_name,
        x,
        y,
        z,
        pointData={"scalar": scalar},
    )
    HOME = os.getenv("HOME")
    cmd = 'rm -f '+HOME+'/.config/ParaView/ParaView-UserSettings.json; pv_zcut3D.py '+file_colormap+' '+pv_fileout_prefix+'; convert -density '+str(convert_density)+' '+pv_fileout_prefix+'.pdf '+pv_fileout_prefix+'.png ' # we need to remove ParaView-UserSettings.json otherwise color palette is taken from that file.

#   print ('cmd=',cmd)
    os.system(cmd)
    return strfield
