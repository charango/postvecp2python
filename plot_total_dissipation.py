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

    X  = np.zeros(len(par.directory))  # forcing frequency
    Yv = np.zeros(len(par.directory))  # viscous dissipation
    Yt = np.zeros(len(par.directory))  # thermal dissipation
    Yc = np.zeros(len(par.directory))  # chemical dissipation
    Yr = np.zeros(len(par.directory))  # power in differential rotation
    Yp = np.zeros(len(par.directory))  # pressure work
    Ye = np.zeros(len(par.directory))  # relative erreur
    
    # loop over directories
    for i in range(len(par.directory)):

        # read dissipation file, which contains (i) viscous
        # dissipation, (ii) power in pressure force, (iii) power in
        # differential rotation and (iv) forcing frequency
        with open(par.directory[i]+'/dissipation', 'r') as f:
            line = f.readline()
        res = line.split()

        Yv[i] = res[0]  # viscous dissipation
        Yt[i] = res[1]  # thermal dissipation
        Yp[i] = res[2]  # pressure work
        Yr[i] = res[3]  # power in differential rotation
        X[i]  = res[4]  # forcing frequency
        Yc[i] = res[5]  # chemical dissipation
        Ye[i] = (Yp[i]+Yv[i]+Yt[i]+Yc[i])/Yp[i]
        
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
    #ax.set_xlim(X.min(),X.max())
    ax.grid(axis='both', which='major', ls='-', alpha=0.8)

    # modes with negative damping rates are shown by a blue filled
    # circle, those with positive damping rates are shown by a red
    # filled circle
    ax.scatter(X,np.abs(Yp),s=10,marker='o',color='tab:orange',label='Pressure work')
    ax.scatter(X,np.abs(Yv),s=10,marker='o',color='tab:blue',label='viscous dissipation')
    ax.scatter(X,np.abs(Yt),s=10,marker='o',color='tab:red',label='|thermal dissipation|')
    ax.scatter(X,np.abs(Yc),s=10,marker='o',color='tab:green',label='chemical dissipation')

    # display global string with main parameters in the bottom
    # use of set_title along with negative pad allows string
    # to be automatically centred in x-position
    ax.set_title(globstr,y=0, pad=-75,fontsize=16,color='black')

    # legend
    from matplotlib.lines import Line2D
    legend = plt.legend(loc='upper left',fontsize=13,facecolor='white',edgecolor='white',framealpha=1.0,numpoints=1,bbox_to_anchor=(0.6,0.93),bbox_transform=plt.gcf().transFigure)
    for line, text in zip(legend.get_lines(), legend.get_texts()):
        text.set_color(line.get_color())
    
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


    # ===================
    # second figure: relative error #
    # ===================
    fig = plt.figure(figsize=(8.,8.))
    plt.subplots_adjust(left=0.16, right=0.98, top=0.96, bottom=0.14)
    ax = fig.gca()
    
    ax.tick_params(top='on', right='on', length = 5, width=1.0, direction='out')
    ax.tick_params(axis='x', which='minor', top=True)
    ax.tick_params(axis='y', which='minor', right=True)

    ax.set_xlabel(r'Forcing frequency $\omega_p$ in inertial frame')
    ax.set_ylabel(r'Relative error')
    ax.grid(axis='both', which='major', ls='-', alpha=0.8)

    ax.scatter(X,Ye,s=10,marker='o',color='tab:blue')

    outfile = 'error_dissipation'+'_'+par.directory[0]+'_'+par.directory[len(par.directory)-1]
    fileout = outfile+'.pdf'
    if par.saveaspdf == 'Yes':
        plt.savefig('./'+fileout, dpi=80)
    if par.saveaspng == 'Yes':
        plt.savefig('./'+re.sub('.pdf', '.png', fileout), dpi=80)
    plt.close(fig)
