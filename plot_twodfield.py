import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import fnmatch
import os
import subprocess
import sys
import re
import math

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, LogFormatter)

from mesh import *
from field import *


def round_to_n(x, n):
    " Round x to n significant figures "
    return round(x, -int(np.floor(np.sign(x) * np.log10(abs(x)))) + n)

def str_fmt(x, n=2):
    " Format x into nice Latex rounding to n"
    power = np.floor(np.log10(round_to_n(x, 0)))
    f_SF = round_to_n(x, n) * pow(10, -power)
    if f_SF != 1.0:
        mystr = "$"+format(f_SF,'.1f')+r"\times 10^{"+format(power,'.0f')+"}$"
        return mystr
    else:
        mystr = r"$10^{"+format(power,'.0f')+"}$"
        return mystr

def fmt_colorbar(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


#   =============
#   MAIN FUNCTION
#   =============
def plottwodfield():

    # first import global variables
    import par
    colored_cmap=par.mycolormap
    
    # how many color bars
    two_color_bars = False
    if par.showdust == 'Yes':
        two_color_bars = True

    # several output numbers
    if isinstance(par.on, int) == True:
        on = [par.on]
    else:
        on = par.on
    # if several directories specified in parameter file, only deal with the first one
    if isinstance(par.directory, str) == False:
        directory = par.directory[0]
    else:
        directory = par.directory
    # if movie requested, then compute fields for all output numbers between range
    # specified in parameter file via on (onmin,onmax)
    if par.movie == 'Yes':
        on = range(par.on[0],par.on[1]+1,par.take_one_point_every)

    for k in range(len(on)):
        if par.movie == 'Yes':
            print('animation: output number '+str(k)+' / '+str(len(on)-1),end='\r')

        if par.allfluids == 'Yes':
            # prepare figure for each output number
            if par.nbfluids > 6:  # could be refined...
                fig, axes = plt.subplots(2,int(np.ceil(par.nbfluids/2)), figsize=(24.,12.))
            else:
                fig, axes = plt.subplots(2,int(np.ceil(par.nbfluids/2)), figsize=(18.,12.))
            if par.fieldofview == 'cart':  
                plt.subplots_adjust(left=0.25, right=0.75, top=0.90, bottom=0.10, wspace=0.6, hspace=0.6)
            else:
                plt.subplots_adjust(left=0.05, bottom=0.08, right=0.95, top=0.92, wspace=0.4, hspace=0.4)
            
        # loop over fluids
        for f in range(len(par.fluids)):

            if par.allfluids == 'Yes':
                if f < int(np.ceil(par.nbfluids/2)):
                    rowindex = 0
                    colindex = f
                else:
                    rowindex = 1
                    colindex = f-int(np.ceil(par.nbfluids/2))
                ax = axes[rowindex,colindex]
                # case we have an odd number of 'windows' to be displayed
                if len(par.fluids) % 2 != 0 and f == len(par.fluids)-1:
                    fig.delaxes(axes[1,int(np.ceil(par.nbfluids/2))-1])

            myfield = Field(field=par.whatfield, fluid=par.fluids[f], on=on[k], directory=directory, physical_units=par.physical_units, nodiff=par.nodiff, fieldofview=par.fieldofview, slice=par.slice, onedprofile='No', override_units=par.override_units)
            array = myfield.data
            
            # plot relative difference wrt initial field 
            if par.nodiff == 'No':
                myfield0 = Field(field=par.whatfield, fluid=par.fluids[f], on=0, directory=directory, physical_units=par.physical_units, nodiff=par.nodiff, fieldofview=par.fieldofview, slice=par.slice, onedprofile='No', override_units=par.override_units)
                array0 = myfield0.data
                for i in range(myfield0.nrad):
                    axisym = np.sum(array0[i,:]) / myfield0.nsec
                    array0[i,:] = axisym
                array = (myfield.data-array0)/array0
            else:
                # conversion in physical units
                if par.physical_units == 'Yes':
                    array = myfield.data * myfield.unit
            
            # rotate field by user-defined angle specified in degree in .dat file
            if ('rotate_angle' in open('paramsf2p.dat').read()) and (par.rotate_angle != '#'):
                shift = int(myfield.nsec*par.rotate_angle/360.0)
                array = np.roll(array,shift=shift,axis=1)

            # roll array by nsec/2 along azimuthal axis if grid's
            # azimuthal extent is smaller than 2pi:
            if abs(myfield.pedge.max()-2.0*np.pi) > 0.1:
                array = np.roll(array,shift=myfield.nsec//2,axis=1)
                
            # get field's full name
            strfield = myfield.strname

            # -------------------
            # Stuff we do only at first output number
            # -------------------
            if k == 0 and f==0:
                # get computational grid coordinates
                # start with cylindrical radius (Fargo2D) or spherical radius (Fargo3D)
                R = myfield.redge # myfield.rmed
            
                # radius conversion in physical units
                if par.physical_units == 'Yes':
                    R *= (myfield.culength / 1.5e11) # in au
            
                # define visualisation params
                if ('rbox' in open('paramsf2p.dat').read()) and (par.rbox != '#'):
                    # start by finding planet's orbita radius
                    if myfield.fargo3d == 'Yes':
                        f1, xpla, ypla, f4, f5, f6, f7, f8, date, omega = np.loadtxt(directory+"/planet0.dat",unpack=True)
                    else:
                        f1, xpla, ypla, f4, f5, f6, f7, date, omega, f10, f11 = np.loadtxt(directory+"/planet0.dat",unpack=True)                    
                else:
                    if (par.myrmin != '#'):
                        myrmin = par.myrmin
                    else:
                        myrmin = R.min()
                    if (par.myrmax != '#'):
                        myrmax = par.myrmax
                    else:
                        myrmax = R.max()
                    imin = np.argmin(np.abs(R-myrmin))
                    imax = np.argmin(np.abs(R-myrmax))   

                # VISUALISATION IN MIDPLANE: get azimuth 
                if par.fargo3d == 'No' or (par.fieldofview != 'latitudinal' and par.fieldofview != 'vertical'):
                    T = myfield.pedge # myfield.pmed
                    myphimin = par.myphimin
                    if (par.myphimin == '#'):
                        myphimin = T.min()
                    jmin = np.argmin(np.abs(T-myphimin))
                    myphimax = par.myphimax
                    if (myphimax == '#'):
                        myphimax = T.max()
                    jmax = np.argmin(np.abs(T-myphimax))
                    if ('flip_xaxis' in open('paramsf2p.dat').read()) and (par.flip_xaxis == 'Yes'):
                        bufmin = myphimin
                        bufmax = myphimax
                        myphimax = bufmin
                        myphimin = bufmax
                
                # VISUALISATION IN VERTICAL (LATITUDINAL) PLANE
                if par.fargo3d == 'Yes' and par.fieldofview == 'latitudinal':
                    #T = myfield.tmed       # latitude
                    T = myfield.tedge       # latitude
                    # number of grid cells in the radial and azimuthal directions
                    jmin = np.argmin(np.abs(T-T.min())) 
                    jmax = np.argmin(np.abs(T-T.max()))
                if par.fargo3d == 'Yes' and par.fieldofview == 'vertical':
                    T = myfield.zedge       # altitude above midplane
                    # number of grid cells in the radial and azimuthal directions
                    jmin = np.argmin(np.abs(T-T.min())) 
                    jmax = np.argmin(np.abs(T-T.max()))

                # color bar for background field
                mycolormap = colored_cmap
                '''
                if par.showdust == 'No':
                    mycolormap = colored_cmap
                else:
                    mycolormap = 'gray'
                '''
                
                # end of stuff done only at first output number
            
            # -------------------
            # case where fields are displayed with a fixed radial range about planet's orbital radius
            # -------------------
            if ('rbox' in open('paramsf2p.dat').read()) and (par.rbox != '#'):
                if par.take_one_point_every == '#':
                    take_one_point_every = 1
                else:
                    take_one_point_every = par.take_one_point_every
                rpla = np.sqrt( xpla[int(on[k])]*xpla[int(on[k])] + ypla[int(on[k])]*ypla[int(on[k])] )
                myrmin = rpla-par.rbox
                myrmax = rpla+par.rbox
                imin = np.argmin(np.abs(R-myrmin))
                imax = np.argmin(np.abs(R-myrmax)) 


            # -------------------
            # read information on the dust particles
            # -------------------
            if par.showdust == 'Yes':
                (rd, td, vrd, vtd, Stokes, sizedust) = np.loadtxt(directory+'/dustsystat'+str(on[k])+'.dat', unpack=True)
                if par.physical_units == 'Yes':
                    rd *= (myfield.culength / 1.5e11)  # in au
                sizemin = par.sizemin
                if (par.sizemin == '#'):
                    sizemin = sizedust.min()
                sizemax = par.sizemax
                if (par.sizemax == '#'):
                    sizemax = sizedust.max()

                # rotate particles azimuth by user-defined angle specified in degree in .dat file
                if ('rotate_angle' in open('paramsf2p.dat').read()) and (par.rotate_angle != '#'):
                    pc_shift_angle = np.pi*par.rotate_angle/180  # in radian
                    for i in range(len(td)):
                        td[i] += pc_shift_angle                  
                        if td[i] > 2.0*np.pi:
                            td[i] -= 2.0*np.pi
                        if td[i] < 0.0:
                            td[i] += 2.0*np.pi
                
                # Only select particles with size between sizemin and sizemax
                mydustsize = sizedust.compress((sizedust<sizemax)&(sizedust>sizemin).flat)
                mydustsize2 = mydustsize[mydustsize.ravel().argsort()]
                myrd = rd.compress((sizedust<sizemax)&(sizedust>sizemin).flat)
                myrd2 = myrd[mydustsize.ravel().argsort()]
                mytd = td.compress((sizedust<sizemax)&(sizedust>sizemin).flat)
                mytd2 = mytd[mydustsize.ravel().argsort()]
                rd = myrd2
                td = mytd2
                sizedust = mydustsize2

            # -------------------
            # POLAR FIELD OF VIEW
            # -------------------
            if par.fieldofview == 'polar':
                # shift by nsec/2 along azimuthal direction, same for dust
                if par.fargo3d == 'No':
                    array = np.roll(array, shift=int(myfield.nsec//2), axis=1)
                array_orig = array
                if par.showdust == 'Yes':
                    td += np.pi
                    for i in range(len(td)):
                        if td[i] > 2.0*np.pi:
                            td[i] -= 2.0*np.pi
                # plot radius in y-axis, azimuth in x-axis
                if par.rvsphi == 'Yes':
                    X = T
                    Y = R
                    if par.showdust == 'Yes':   # particles
                        xdust = td
                        ydust = rd
                else:
                # plot azimuth in y-axis, radius in x-axis
                    array = np.transpose(array)
                    X = R
                    Y = T
                    if par.showdust == 'Yes':   # particles
                        xdust = rd
                        ydust = td
                        
                # figure
                if par.allfluids == 'No':
                    if two_color_bars == True:
                        fig = plt.figure(figsize=(8.5,8.))        
                        plt.subplots_adjust(left=0.12, right=0.86, top=0.88, bottom=0.11)
                    else:
                        fig = plt.figure(figsize=(8.,8.))
                        plt.subplots_adjust(left=0.16, right=0.94, top=0.88, bottom=0.11)
                    ax = plt.gca()
                    
                if par.physical_units == 'No':
                    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
                    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
                if par.rvsphi == 'Yes':
                    ax.set_xlabel('Azimuth [rad]')
                    if par.physical_units == 'Yes':
                        ax.set_ylabel('Radius [au]')
                    else:
                        ax.set_ylabel('Radius [code units]')
                    xlim_min = myphimin
                    xlim_max = myphimax
                    ylim_min = myrmin
                    ylim_max = myrmax
                    #ax.set_xlim(myphimin,myphimax)
                    #ax.set_ylim(myrmin,myrmax)
                else:
                    ax.set_ylabel('Azimuth [rad]')
                    if par.physical_units == 'Yes':
                        ax.set_xlabel('Radius [au]')
                    else:
                        ax.set_xlabel('Radius [code units]')
                    xlim_min = myrmin
                    xlim_max = myrmax
                    ylim_min = myphimin
                    ylim_max = myphimax
                    #ax.set_ylim(myphimin,myphimax)
                    #ax.set_xlim(myrmin,myrmax)
            #
            # -------------------------
            # LATITUDINAL FIELD OF VIEW
            # -------------------------
            if par.fieldofview == 'latitudinal':
                X = R
                Y = T
                array_orig = array
                array = np.transpose(array)

                if par.allfluids == 'No':
                    if two_color_bars == True:
                        fig = plt.figure(figsize=(8.5,8.))        
                        plt.subplots_adjust(left=0.19, right=0.86, top=0.88, bottom=0.11)
                    else:
                        fig = plt.figure(figsize=(8.,8.))
                        plt.subplots_adjust(left=0.19, right=0.94, top=0.88, bottom=0.11)
                    ax = plt.gca()
                    
                if par.physical_units == 'No':
                    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
                    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
                    ax.set_ylabel('Latitude [rad]')
                    ax.set_xlabel('Radius [code units]')                                
                else:
                    ax.set_xlabel('Radius [au]')
                    ax.set_ylabel('Latitude [rad]')
                xlim_min = myrmin
                xlim_max = myrmax
                ylim_min = T.min()
                ylim_max = T.max()
                #ax.set_ylim(T.min(),T.max())
                #ax.set_xlim(myrmin,myrmax)
            #
            # -----------------------
            # CARTESIAN FIELD OF VIEW
            # -----------------------
            if par.fieldofview == 'cart':            
                if par.fargo3d == 'Yes' and myfield.cartesian_grid == 'No':
                    array = np.roll(array, shift=int(myfield.nsec//2), axis=1)
                array_orig = array
                
                if myfield.cartesian_grid == 'No':
                    radius_matrix, theta_matrix = np.meshgrid(R,myfield.pedge)
                    X = radius_matrix * np.cos(theta_matrix)
                    if ('flip_xaxis' in open('paramsf2p.dat').read()) and (par.flip_xaxis == 'Yes'):
                        X = -X
                    Y = radius_matrix * np.sin(theta_matrix)
                    # roll array by nsec/2 along azimuthal axis if grid's
                    # azimuthal extent is smaller than 2pi:
                    if abs(myfield.pedge.max()-2.0*np.pi) > 0.1:
                        array = np.roll(array,shift=myfield.nsec//2,axis=1)
                    #
                    array = np.transpose(array)
                    #
                    if par.showdust == 'Yes':    # particles
                        xdust = rd*np.cos(td)
                        if ('flip_xaxis' in open('paramsf2p.dat').read()) and (par.flip_xaxis == 'Yes'):
                            xdust = -xdust
                        ydust = rd*np.sin(td)

                if myfield.cartesian_grid == 'Yes':
                    X = myfield.xmed
                    Y = myfield.ymed
                    
                # figure
                if par.allfluids == 'No':
                    if two_color_bars == True:
                        fig = plt.figure(figsize=(8.8,8.))
                        # do not edit subplot position below!
                        plt.subplots_adjust(left=0.14, right=0.85, top=0.88, bottom=0.1)
                    else:
                        fig = plt.figure(figsize=(8.,8.))
                        # do not edit subplot position below!
                        plt.subplots_adjust(left=0.17, right=0.92, top=0.88, bottom=0.1)
                    ax = plt.gca()

                xlim_min = -myrmax
                xlim_max = myrmax
                ylim_min = -myrmax
                ylim_max = myrmax
                #ax.set_xlim(-myrmax,myrmax)
                #ax.set_ylim(-myrmax,myrmax)
                if par.physical_units == 'Yes':
                    ax.set_xlabel('x [au]')
                    ax.set_ylabel('y [au]')
                else:
                    ax.set_xlabel('x [code units]')
                    ax.set_ylabel('y [code units]')
            #
            # -----------------------
            # VERTICAL FIELD OF VIEW
            # -----------------------
            if par.fieldofview == 'vertical':
                array_orig = array
                array = np.transpose(array)
                radius_matrix, theta_matrix = np.meshgrid(R,myfield.tedge)
                #radius_matrix, theta_matrix = np.meshgrid(myfield.rmed,myfield.tmed)
                X = radius_matrix * np.cos(theta_matrix)
                Y = radius_matrix * np.sin(theta_matrix)

                # figure
                if par.allfluids == 'No':
                    if two_color_bars == True:
                        fig = plt.figure(figsize=(8.5,8.))        
                        plt.subplots_adjust(left=0.15, right=0.86, top=0.88, bottom=0.11)
                    else:
                        fig = plt.figure(figsize=(8.,8.))
                        plt.subplots_adjust(left=0.15, right=0.94, top=0.88, bottom=0.11)
                    ax = plt.gca()
                    
                xlim_min = myrmin #X.min()
                xlim_max = myrmax #X.max()
                ylim_min = Y.min()
                ylim_max = Y.max()
                #ax.set_ylim(Y.min(),Y.max())
                #ax.set_xlim(X.min(),X.max())
                if par.physical_units == 'No':
                    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
                    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
                    ax.set_ylabel('Altitude [code units]')
                    ax.set_xlabel('Radius [code units]')                                
                else:
                    ax.set_xlabel('Radius [au]')
                    ax.set_ylabel('Altitude [au]')

            # common display
            ax.set_xlim(xlim_min,xlim_max)
            ax.set_ylim(ylim_min,ylim_max)
            ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
            ax.tick_params(axis='x', which='minor', top=True)
            ax.tick_params(axis='y', which='minor', right=True)
            ax.xaxis.set_major_locator(plt.MaxNLocator(5))
            ax.yaxis.set_major_locator(plt.MaxNLocator(5))

            # -------------------
            # adjust min and max in color bar according to field of view
            # -------------------
            #print(array_orig.shape,imin,imax,jmin,jmax)
            myfieldmin = par.fieldmin
            if (par.fieldmin == '#'):
                myfieldmin = array_orig[imin:imax+1,jmin:jmax+1].min()
                if par.fieldofview == 'latitudinal':
                    myfieldmin = array_orig[imin:imax+1,jmin:jmax+1].min()
            myfieldmax = par.fieldmax
            if (par.fieldmax == '#'):
                myfieldmax = array_orig[imin:imax+1,jmin:jmax+1].max()
                if par.fieldofview == 'latitudinal':
                    myfieldmax = array_orig[imin:imax+1,jmin:jmax+1].max()

            print('imin, imax = ', imin, imax)
            print('jmax, jmax = ', jmin, jmax)
            print('fieldmin, fieldmax = ', myfieldmin, myfieldmax)
            print('array_orig shape = ', array_orig.shape)
            if par.log_colorscale == 'Yes':
                if (par.fieldmin == 'auto' or par.fieldmax == 'auto'):
                    minarray = array.min() #1e-3*array.max()?
                    maxarray = array.max()
                    array = np.log(array/minarray)/np.log(maxarray/minarray)
                    #print(par.fluids[f],minarray,maxarray,array.min(),array.max())
                    #print(array)
                    myfieldmin = 1e-5
                    myfieldmax = 1.0
                    mynorm = matplotlib.colors.Normalize(vmin=myfieldmin,vmax=myfieldmax)
                    strfield = r' $\log(\rho/\rho_{min})~/~\log(\rho_{max}/\rho_{min})$'+' at '+myfield.strtime
                else:
                    if (myfieldmax/myfieldmin > 1e3 and (par.fieldmin == '#')):
                        myfieldmin = 1e-3*myfieldmax

            # Normalization for colorbar: linear or logarithmic scale
            if par.log_colorscale == 'Yes':
                if (par.fieldmin != 'auto' and par.fieldmax != 'auto'):
                    mynorm = matplotlib.colors.LogNorm(vmin=myfieldmin,vmax=myfieldmax)
            else:
                mynorm = matplotlib.colors.Normalize(vmin=myfieldmin,vmax=myfieldmax)
                
            # -----------------------
            # display contour field
            # -----------------------
            print('shapes = ', X.shape, Y.shape, array.shape)
            CF = ax.pcolormesh(X,Y,array,cmap=mycolormap,norm=mynorm,rasterized=True)
            #CF = ax.imshow(array, origin='lower', cmap=mycolormap, interpolation='bilinear', vmin=myfieldmin, vmax=myfieldmax, aspect='auto', extent=[X.min(),X.max(),Y.min(),Y.max()])

            # ------------------
            # overlay streamlines
            # ------------------
            if par.streamlines == 'Yes':
                for s in range(par.nstreamlines):
                    # initial radius and azimuth of streamlines
                    myR0 = myrmin   + np.random.rand()*(myrmax-myrmin)
                    myT0 = myphimin + np.random.rand()*(myphimax-myphimin)
                    # forward integration of streamlines
                    xstr,ystr = myfield.compute_streamline(niterations=10000,R0=myR0,T0=myT0,rmin=myrmin,rmax=myrmax,pmin=myphimin,pmax=myphimax,forward=True,fieldofview=par.fieldofview,slice=par.slice)
                    ax.scatter(xstr,ystr,s=3,marker='o',color='white')
                    # backward integration of streamlines
                    xstr,ystr = myfield.compute_streamline(niterations=10000,R0=myR0,T0=myT0,rmin=myrmin,rmax=myrmax,pmin=myphimin,pmax=myphimax,forward=False,fieldofview=par.fieldofview,slice=par.slice)
                    ax.scatter(xstr,ystr,s=3,marker='o',color='white')
            
            # ------------------
            # overlay particles
            # ------------------
            if par.showdust == 'Yes':
                mynorm = matplotlib.colors.LogNorm(vmin=sizemin,vmax=sizemax)
                CD = ax.scatter(xdust,ydust,s=1,c=sizedust,cmap='nipy_spectral',alpha=0.3,norm=mynorm)
                #CD = ax.scatter(xdust,ydust,s=1,c=sizedust,cmap='nipy_spectral',alpha=0.3,vmin=sizemin,vmax=sizemax,norm=matplotlib.colors.LogNorm())

            # ------------------
            # overlay planets
            # ------------------
            if par.showplanet == 'Yes':

                # Find how many 'planets' there are
                nbplanets = len(fnmatch.filter(os.listdir(directory), 'planet*.dat'))
                xp = np.zeros(nbplanets)
                yp = np.zeros(nbplanets)

                for l in range(nbplanets):
                    # read information on the planets (inner one so far)
                    if par.fargo3d == 'Yes':
                        f1, xpla, ypla, f4, f5, f6, f7, f8, date, omega = np.loadtxt(directory+"/planet"+str(l)+".dat",unpack=True)
                    else:
                        if par.fargo_orig == 'No':
                            f1, xpla, ypla, f4, f5, f6, f7, date, omega, f10, f11 = np.loadtxt(directory+"/planet"+str(l)+".dat",unpack=True)
                        else:
                            f1, xpla, ypla, f4, f5, f6, f7, date, omega = np.loadtxt(directory+"/planet"+str(l)+".dat",unpack=True)
                    xp[l] = xpla[on[k]]
                    yp[l] = ypla[on[k]]
                    if par.physical_units == 'Yes':
                        xp[l] *= (myfield.culength / 1.5e11)  # in au
                        yp[l] *= (myfield.culength / 1.5e11)  # in au
                    rp = np.sqrt(xp[l]*xp[l] + yp[l]*yp[l])  # planet's orbital radius
                    tp = math.atan2(yp[l],xp[l])       # planet's azimuthal angle

                    # rotate particles azimuth by user-defined angle specified in degree in .dat file
                    if ('rotate_angle' in open('paramsf2p.dat').read()) and (par.rotate_angle != '#'):
                        planet_shift_angle = np.pi*par.rotate_angle/180  # in radian
                        tp += planet_shift_angle                  
                        if tp > 2.0*np.pi:
                            tp -= 2.0*np.pi
                        if tp < 0.0:
                            tp += 2.0*np.pi
                        xp[l] = rp*np.cos(tp)
                        yp[l] = rp*np.sin(tp)

                    if par.fieldofview == 'polar':
                        tp += np.pi
                        if tp > 2.0*np.pi:
                            tp -= 2.0*np.pi
                        if par.rvsphi == 'Yes':
                            xp[l] = tp
                            yp[l] = rp
                        else:
                            xp[l] = rp
                            yp[l] = tp

                    if par.fieldofview == 'cart' and myfield.cartesian_grid == 'No':
                        if ('flip_xaxis' in open('paramsf2p.dat').read()) and (par.flip_xaxis == 'Yes'):
                            xp[l] = -xp[l]
                    
                # add planets via scatter plot
                if par.verbose == 'Yes':
                    print('xp = ', xp, ' , yp = ', yp)
                CP = ax.scatter(xp,yp,s=10,c='lightpink',cmap=colored_cmap,alpha=1)


            # ------------------
            # overlay CPU
            # ------------------
            if ('showcpus' in open('paramsf2p.dat').read()) and par.showcpus == 'Yes':
                cpunb, cpurmin, cpurmax = np.loadtxt(directory+"/minmaxradii.dat",unpack=True)
                if par.physical_units == 'Yes':
                    cpurmin *= (myfield.culength / 1.5e11) # in au  
                    cpurmax *= (myfield.culength / 1.5e11) # in au              
                for i in range(len(cpurmin)):
                    if par.fieldofview == 'polar':
                        if par.rvsphi == 'No':
                            ax.plot([cpurmin[i],cpurmin[i]],[Y.min(),Y.max()],'-',linewidth=1,color='grey')
                        else:
                            ax.plot([X.min(),X.max()],[cpurmin[i],cpurmin[i]],'-',linewidth=1,color='grey')

            # ------------------        
            # Add user-defined string in top right corner
            # ------------------
            if ( ('display_label' in open('paramsf2p.dat').read()) and (par.display_label != '#') ):
                mylabel = par.display_label.replace('_',' ')
                xlabel = xlim_max - 0.05*(xlim_max-xlim_min)
                ylabel = ylim_max - 0.05*(ylim_max-ylim_min)
                #print('xlabel=',xlabel, ' ylabel=', ylabel)
                ax.text(xlabel, ylabel, mylabel, fontsize=15, color = 'white',weight='bold', horizontalalignment='right')

            # ----------------
            # special case all fluids 
            # ----------------
            if par.allfluids == 'Yes':
                colorstr = 'lightpink' # 'white'
                xmin,xmax = ax.get_xlim()
                ymin,ymax = ax.get_ylim()
                if f > 0:
                    strsize = str_fmt(par.dust_size[f-1])
                    strsize += ' m'
                    ax.text(xmax,ymin,strsize,fontsize=20,color=colorstr,horizontalalignment='right',verticalalignment='bottom')
                else:
                    ax.text(xmax,ymin,'gas',fontsize=20,color=colorstr,horizontalalignment='right',verticalalignment='bottom')

            # ----------------
            # plot color-bars
            # ----------------
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("top", size="2.5%", pad=0.12)
            cb = plt.colorbar(CF, cax=cax, orientation='horizontal')
            cax.xaxis.tick_top()
            cax.xaxis.set_tick_params(direction='out')
            cax.xaxis.set_major_locator(plt.MaxNLocator(4))
            # title on top
            cax.xaxis.set_label_position('top')
            cax.set_xlabel(strfield)
            cax.xaxis.labelpad = 8
            if par.log_colorscale == 'Yes' and par.fieldmin != 'auto' and par.fieldmax != 'auto':
                cax.xaxis.set_major_locator(ticker.LogLocator(base=10.0,numticks=8))
            if par.log_colorscale == 'Yes' and (par.fieldmin == 'auto' or par.fieldmax == 'auto'):
                cax.xaxis.set_tick_params(direction='out')
            if two_color_bars == True:
                cax2 = divider.append_axes("right", size="2.5%", pad=0.12)
                cb2 = plt.colorbar(CD, cax=cax2, orientation='vertical')
                cax2.yaxis.tick_right()        
                cb2.set_alpha(0.7)
                from packaging import version
                # if matplotlib version >= 3.6, draw_all no longer
                # exists and becomes _draw_all
                if version.parse(matplotlib.__version__) < version.parse("3.6"):
                    cb2.draw_all()
                else:
                    cb2._draw_all()
                cax2.yaxis.set_tick_params(direction='out')
                # title on right-hand side
                cax2.set_ylabel('dust size [meter]')
                cax2.yaxis.labelpad = 8

        # ------------------
        # save in pdf or png files
        # ------------------

        if ('filename' in open('paramsf2p.dat').read()) and (par.filename != '#'):
            outfile = par.filename
        else:
            outfile = par.fluid+'_'+par.whatfield+'_'+directory+'_'+par.fieldofview+'_'+par.slice+'_'+str(on[k]).zfill(4)
        if par.movie == 'Yes' and par.take_one_point_every != 1:
            outfile = par.fluid+'_'+par.whatfield+'_'+directory+'_'+par.fieldofview+'_'+par.slice+'_'+str(k).zfill(4)
        if par.showdust == 'Yes':
            outfile += '_dust'
        fileout = outfile+'.pdf'
        if par.saveaspdf == 'Yes':
            plt.savefig('./'+fileout, dpi=80)
        if par.saveaspng == 'Yes':
            plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=80)  #80 is optimum
        if par.movie == 'Yes':
            plt.close(fig)  # close figure as we reopen figure at every output number

    
    # ------------------
    # finally concatenate png if movie requested
    # ------------------
    if par.movie == 'Yes':
        # png files that have been created above
        allpngfiles = [par.fluid+'_'+par.whatfield+'_'+directory+'_'+par.fieldofview+'_'+par.slice+'_'+str(on[x]).zfill(4)+'.png' for x in range(len(on))]
        if par.take_one_point_every != 1:
            allpngfiles = [par.fluid+'_'+par.whatfield+'_'+directory+'_'+par.fieldofview+'_'+par.slice+'_'+str(x).zfill(4)+'.png' for x in range(len(on))]
            str_on_start_number = str(0)
        else:
            str_on_start_number = str(on[0])
        # input files for ffpmeg
        input_files = par.fluid+'_'+par.whatfield+'_'+directory+'_'+par.fieldofview+'_'+par.slice+'_%04d.png'
        # output file for ffmpeg
        filempg = par.fluid+'_'+par.whatfield+'_'+directory+'_'+par.fieldofview+'_'+par.slice+'_'+str(on[0])+'_'+str(on[len(on)-1])+'.mpg'
        # options
        if par.showdust == 'Yes':
            allpngfiles = [par.fluid+'_'+par.whatfield+'_'+directory+'_'+par.fieldofview+'_'+par.slice+'_'+str(on[x]).zfill(4)+'_dust.png' for x in range(len(on))]
            if par.take_one_point_every != 1:
                allpngfiles = [par.fluid+'_'+par.whatfield+'_'+directory+'_'+par.fieldofview+'_'+par.slice+'_'+str(x).zfill(4)+'_dust.png' for x in range(len(on))]
            input_files = par.fluid+'_'+par.whatfield+'_'+directory+'_'+par.fieldofview+'_'+par.slice+'_%04d_dust.png'
            filempg = re.sub('.mpg', '_dust.mpg', filempg)
        if par.nodiff == 'Yes':
            filempg = re.sub('.mpg', '_nodiff.mpg', filempg)
        # call to ffmpeg-python (you also need to install ffmpeg on your local environement!)
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
        
