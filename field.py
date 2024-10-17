import numpy as np
import os
import sys
import subprocess

import math
import itertools

from mesh import *

# ----------------------------------------
# reading fields output by Fargo (density, energy, velocities, etc) or
# constructing fields from outputs
# ----------------------------------------

class Field(Mesh):
    # based on P. Benitez Llambay routine
    """
    Field class, it stores all the mesh, parameters and scalar data 
    for a scalar field.
    Input: field [string] -> name of the field
           staggered='c' [string] -> staggered direction of the field. 
                                      Possible values: 'x', 'y', 'xy', 'yx'
           directory='' [string] -> where filename is
           on='0' [integer] -> output number
           fluid = '' [string] -> gas or dust
           dtype='float64' (numpy dtype) -> 'float64', 'float32', 
                                             depends if FARGO_OPT+=-DFLOAT is activated
    """
    def __init__(self, field, fluid='gas', staggered='c', directory='', on=0, dtype='float64', physical_units='Yes', nodiff='Yes', fieldofview='polar', slice='midplane', onedprofile='No', override_units='No'):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'

        # first import global variables
        import par
    
        # check if simulation has been carried out with fargo3d or
        # with fargo2d / original fargo code (somewhat redundant with
        # what is done in par.py...)
        if isinstance(directory, str) == False:
            summary0_file = directory[0]+'/summary0.dat'
            usedazi_file  = directory[0]+'/used_azi.dat'
        else:
            summary0_file = directory+'/summary0.dat'
            usedazi_file  = directory+'/used_azi.dat'
        if os.path.isfile(summary0_file) == True:
            # Simulations were carried out with Fargo3D
            self.fargo3d = 'Yes'
        else:
            # Simulations were carried out with Fargo2D
            self.fargo3d = 'No'
            if os.path.isfile(usedazi_file) == True:
                # Simulations were carried out with Dusty FARGO-ADSG
                self.fargo_orig = 'No'
            else:
                # Simulations were carried out with original FARGO code
                self.fargo_orig = 'Yes'
            
        self.cartesian_grid = 'No'
        if self.fargo3d == 'Yes':
            command = par.awk_command+' " /^COORDINATES/ " '+directory+'/variables.par'
            if sys.version_info[0] < 3:   # python 2.X
                buf = subprocess.check_output(command, shell=True)
            else:                         # python 3.X
                buf = subprocess.getoutput(command)
            if buf.split()[1] == 'cartesian':
                self.cartesian_grid = 'Yes'
                #print('A cartesian grid has been used in the simulation')
            if field == 'vrad' and self.cartesian_grid == 'No':
                field = 'vy'
            if field == 'vtheta' and self.cartesian_grid == 'No':
                field = 'vx'
            if field == 'vcol':
                field = 'vz'
            command = par.awk_command+' " /^ZMAX/ " '+directory+'/variables.par'
            if sys.version_info[0] < 3:   # python 2.X
                buf = subprocess.check_output(command, shell=True)
            else:                         # python 3.X
                buf = subprocess.getoutput(command)
            self.zmax = float(buf.split()[1])
                                
        # get nrad and nsec (number of cells in radial and azimuthal directions)
        buf, buf, buf, buf, buf, buf, nrad, nsec = np.loadtxt(directory+"dims.dat",unpack=True)
        self.nrad = int(nrad)
        self.nsec = int(nsec)

        # get number of cells in (co)latitude, if any
        if self.fargo3d == 'Yes':
            command = par.awk_command+' " /^NZ/ " '+directory+'variables.par'
            # check which version of python we're using
            if sys.version_info[0] < 3:   # python 2.X
                buf = subprocess.check_output(command, shell=True)
            else:                         # python 3.X
                buf = subprocess.getoutput(command)
            self.ncol = int(buf.split()[1])
        else:
            self.ncol = 1
            
        # all Mesh attributes
        Mesh.__init__(self, directory)    

        if physical_units == 'Yes':
            if self.fargo3d == 'No' and self.fargo_orig == 'No':
                # units.dat contains physical units of mass [kg], length [m], time [s], and temperature [k] 
                cumass, culength, cutime, cutemp = np.loadtxt(directory+"units.dat",unpack=True)
                self.cumass = cumass
                self.culength = culength
                self.cutime = cutime 
                self.cutemp = cutemp
            if self.fargo3d == 'Yes':
                # get units via variable.par file
                command = par.awk_command+' " /^UNITOFLENGTHAU/ " '+directory+'variables.par'
                # check which version of python we're using
                if sys.version_info[0] < 3:   # python 2.X
                    buf = subprocess.check_output(command, shell=True)
                else:                         # python 3.X
                    buf = subprocess.getoutput(command)
                self.culength = float(buf.split()[1])*1.5e11  #from au to meters
                command = par.awk_command+' " /^UNITOFMASSMSUN/ " '+directory+'variables.par'
                # check which version of python we're using
                if sys.version_info[0] < 3:   # python 2.X
                    buf = subprocess.check_output(command, shell=True)
                else:                         # python 3.X
                    buf = subprocess.getoutput(command)
                self.cumass = float(buf.split()[1])*2e30  #from Msol to kg
                # unit of time = sqrt( pow(L,3.) / 6.673e-11 / M );
                self.cutime = np.sqrt( self.culength**3.0 / 6.673e-11 / self.cumass)
                # unit of temperature = mean molecular weight * 8.0841643e-15 * M / L;
                self.cutemp = 2.35 * 8.0841643e-15 * self.cumass / self.culength

            if override_units == 'Yes':
                if par.new_unit_length == 0.0:
                    sys.exit('override_units set to yes but new_unit_length is not defined in params.dat, I must exit!')
                #else:
                    #print('new unit of length in meters : ', par.new_unit_length)
                if par.new_unit_mass == 0.0:
                    sys.exit('override_units set to yes but new_unit_mass is not defined in params.dat, I must exit!')
                #else:
                    #print('new unit of mass in kg : ', par.new_unit_mass)
                self.cumass = par.new_unit_mass
                self.culength = par.new_unit_length
                # Deduce new units of time and temperature:
                # T = sqrt( pow(L,3.) / 6.673e-11 / M )
                # U = mmw * 8.0841643e-15 * M / L;
                self.cutime = np.sqrt( self.culength**3 / 6.673e-11 / self.cumass)
                self.cutemp = 2.35 * 8.0841643e-15 * self.cumass / self.culength
                '''
                print('### NEW UNITS SPECIFIED: ###')
                print('new unit of length [m] = ',self.culength)
                print('new unit of mass [kg]  = ',self.cumass)
                print('new unit of time [s] = ',self.cutime)
                print('new unit of temperature [K] = ',self.cutemp)
                '''

        # now, staggering:
        if staggered.count('r')>0:
            self.r = self.redge[:-1] # do not dump last element
        else:
            self.r = self.rmed

        # get time in orbital periods at initial location of inner planet
        # also get omegaframe for vtheta
        if self.fargo3d == 'Yes':
            f1, xpla, ypla, f4, f5, f6, f7, f8, date, omega = np.loadtxt(directory+"/planet0.dat",unpack=True)
        else:
            if self.fargo_orig == 'Yes':
                f1, xpla, ypla, f4, f5, f6, f7, date, omega = np.loadtxt(directory+"planet0.dat",unpack=True)
            else:
                f1, xpla, ypla, f4, f5, f6, f7, date, omega, f10, f11 = np.loadtxt(directory+"planet0.dat",unpack=True)

        # read only first line of file orbit0.dat to get initial semi-major axis:
        # procedure is independent of how many columns there is in the file
        if os.path.isfile(directory+"/planet0.dat") == True:
            with open(directory+"/planet0.dat") as f_in:
                firstline_orbitfile = np.genfromtxt(itertools.islice(f_in, 0, 1, None), dtype=float)
            apla = firstline_orbitfile[1]
        else:
            with open(directory+"/orbit0.dat") as f_in:
                firstline_orbitfile = np.genfromtxt(itertools.islice(f_in, 0, 1, None), dtype=float)
            apla = firstline_orbitfile[2]
        if (apla <= 1e-5):
            apla = 1.0
            
        # check if planet0.dat file has only one line or more!
        if isinstance(xpla, (list, tuple, np.ndarray)) == True:
            omegaframe = omega[on]
            #rpla_0 = np.sqrt( xpla[0]*xpla[0] + ypla[0]*ypla[0] )        
            #time_in_code_units = round(date[on]/2./np.pi/rpla_0/np.sqrt(rpla_0),1)
            time_in_code_units = round(date[on]/2./np.pi/apla/np.sqrt(apla),1)
        else:
            omegaframe = omega
            #rpla_0 = np.sqrt( xpla*xpla + ypla*ypla )        
            #time_in_code_units = round(date/2./np.pi/rpla_0/np.sqrt(rpla_0),1)
            time_in_code_units = round(date/2./np.pi/apla/np.sqrt(apla),1)
        self.strtime = str(time_in_code_units)+' orbits'
        
        # --------- NB: it can take WAY more time to read orbit0.dat
        # than planet0.dat since the former file can contain many more
        # lines! Also, if the planet has no eccentricity to start
        # with, then r = a initially and it does not change anything
        # for the calculation of the orbital time. If not, it should
        # be slightly edited... -------
        
        if nodiff == 'No':
            self.strname = 'perturbed '+fluid
        else:
            self.strname = fluid

        # check if input file exists: if so, read it, otherwise build array of desired quantity
        input_file = directory+fluid+field+str(on)+'.dat'

        if field == 'temp' and self.fargo3d == 'No':
            input_file = directory+'Temperature'+str(on)+'.dat'

        if os.path.isfile(input_file) == True:
            self.data = self.__open_field(input_file,dtype,fieldofview,slice)
            if (field == 'vtheta' or field == 'vx') and self.cartesian_grid == 'No':
                for i in range(self.nrad):
                    self.data[i,:] += (self.rmed)[i]*omegaframe
                print(self.data.min(), self.data.max())
            '''
            if field == 'dens' and 'cavity_gas' in open('paramsf2p.dat').read() and cavity_gas == 'Yes':
                imin = np.argmin(np.abs(self.rmed-1.3))
                for i in range(self.nrad):
                    if i < imin:
                        for j in range(self.nsec):
                            self.data[i,j] *= ((self.rmed[i]/self.rmed[imin])**(6.0))
            '''
            # ----
            # print out disc mass if density field is computed
            # ----
            if field == 'dens' and par.verbose == 'Yes' and ( par.fieldofview == 'cartesian' or par.fieldofview == 'polar'):
                mass    = np.zeros((self.nrad,self.nsec))
                surface = np.zeros((self.nrad,self.nsec))

                Rinf = self.redge[0:len(self.redge)-1]
                Rsup = self.redge[1:len(self.redge)]
                surf = np.pi * (Rsup*Rsup - Rinf*Rinf) / self.nsec
                for th in range(self.nsec):
                    surface[:,th] = surf

                # mass of each grid cell
                mass = self.data*surface
                # total disc mass 
                print('disc mass / star mass = ', np.sum(mass)) 

            # ----
            # PASSIVE SCALAR
            # ----
            if field == 'label':
                self.strname = 'concentration'
                    
        else:
            #print('input file ',input_file,' does not exist!')
            self.data = np.zeros((self.nrad,self.nsec))  # default
            self.unit = 1.0                              # default
            # ----
            # TEMPERATURE or PRESSURE or TOOMRE Q-parameter = c_s Omega / pi G Sigma
            # checked by CB on January 2022
            # ----
            if field == 'temp' or field == 'pressure' or field == 'toomre':
                if field == 'pressure':
                    self.strname += ' pressure'
                if field == 'toomre':
                    self.strname += ' Toomre parameter'
                    
                # check that no energy equation was employed
                if self.fargo3d == 'No':
                    command = par.awk_command+' " /^EnergyEquation/ " '+directory+'*.par'
                    # check which version of python we're using
                    if sys.version_info[0] < 3:   # python 2.X
                        buf = subprocess.check_output(command, shell=True)
                    else:                         # python 3.X
                        buf = subprocess.getoutput(command)
                    energyequation = str(buf.split()[1])
                    
                    if energyequation == 'No':
                        # get the aspect ratio and flaring index used in the numerical simulation
                        command = par.awk_command+' " /^AspectRatio/ " '+directory+'*.par'
                        # check which version of python we're using
                        if sys.version_info[0] < 3:   # python 2.X
                            buf = subprocess.check_output(command, shell=True)
                        else:                         # python 3.X
                            buf = subprocess.getoutput(command)
                        aspectratio = float(buf.split()[1])
                        # get the flaring index used in the numerical simulation
                        command = par.awk_command+' " /^FlaringIndex/ " '+directory+'*.par'
                        if sys.version_info[0] < 3:
                            buf = subprocess.check_output(command, shell=True)
                        else:
                            buf = subprocess.getoutput(command)
                        flaringindex = float(buf.split()[1])
                    
                    # get the adiabatic index index used in the numerical simulation
                    command = par.awk_command+' " /^AdiabaticIndex/ " '+directory+'*.par'
                    # check which version of python we're using
                    if sys.version_info[0] < 3:   # python 2.X
                        buf = subprocess.check_output(command, shell=True)
                    else:                         # python 3.X
                        buf = subprocess.getoutput(command)
                    gamma = float(buf.split()[1])

                    if field == 'toomre':
                        vphi = self.__open_field(directory+fluid+'vtheta'+str(on)+'.dat',dtype,fieldofview,slice)
                    
                # case we're running with Fargo3D    
                else:
                    if "ISOTHERMAL" in open(directory+'summary'+str(on)+'.dat',"r").read():
                        energyequation = "No"
                        if energyequation == 'No':
                            command = par.awk_command+' " /^ASPECTRATIO/ " '+directory+'*.par'
                            if sys.version_info[0] < 3:   # python 2.X
                                buf = subprocess.check_output(command, shell=True)
                            else:                         # python 3.X
                                buf = subprocess.getoutput(command)
                            aspectratio = float(buf.split()[1])                
                            # then get the flaring index used in the numerical simulation
                            command = par.awk_command+' " /^FLARINGINDEX/ " '+directory+'*.par'
                            if sys.version_info[0] < 3:
                                buf = subprocess.check_output(command, shell=True)
                            else:
                                buf = subprocess.getoutput(command)
                            flaringindex = float(buf.split()[1])                        
                    else:
                        energyequation = "Yes"
                     # then get the adiabatic index 'gamma' in the numerical simulation
                    command = par.awk_command+' " /^GAMMA/ " '+directory+'*.par'
                    if sys.version_info[0] < 3:
                        buf = subprocess.check_output(command, shell=True)
                    else:
                        buf = subprocess.getoutput(command)
                    gamma = float(buf.split()[1])          

                    if field == 'toomre':
                        vphi = self.__open_field(directory+fluid+'vx'+str(on)+'.dat',dtype,fieldofview,slice)
                        
                # work out temperature first: self.data contains the gas temperature
                if energyequation == 'No':
                    gamma = 1.0  # reset gamma to 1
                    # expression below valid both for Fargo-3D and Dusty FARGO-ADSG:
                    for i in range(self.nrad):
                        self.data[i,:] = aspectratio*aspectratio*(((self.rmed)[i])**(-1.0+2.0*flaringindex))
                else:
                    if self.fargo3d == 'No':
                        self.data = self.__open_field(directory+'Temperature'+str(on)+'.dat',dtype,fieldofview,slice)
                    else:
                        cs = self.__open_field(directory+fluid+'energy'+str(on)+'.dat',dtype,fieldofview,slice)
                        self.data = np.sqrt(cs)

                # work out pressure then
                if field == 'pressure':
                    dens = self.__open_field(directory+fluid+'dens'+str(on)+'.dat',dtype,fieldofview,slice)
                    self.data *= dens

                # finally work out Toomre Q-parameter
                if field == 'toomre':
                    dens = self.__open_field(directory+fluid+'dens'+str(on)+'.dat',dtype,fieldofview,slice)
                    cs = np.sqrt(gamma*self.data)
                    self.data = cs/np.pi/dens  # (nrad,nsec)
                    
                    omega = np.zeros((self.nrad,self.nsec))
                    # vphi, read above, is in the corotating frame!
                    for i in range(self.nrad):
                        vphi[i,:] += (self.rmed)[i]*omegaframe
                    axivphi = np.sum(vphi,axis=1)/self.nsec  # just in case...
                    
                    for i in range(self.nrad):
                        omega[i,:] = vphi[i,:] / self.rmed[i]

                    self.data *= omega

            #
            # ----
            # SPECIFIC ENTROPY Sigma^{1-gamma} x Temp
            # ----
            if field == 'entropy':
                # get the adiabatic index index used in the numerical simulation
                command = par.awk_command+' " /^AdiabaticIndex/ " '+directory+'*.par'
                # check which version of python we're using
                if sys.version_info[0] < 3:   # python 2.X
                    buf = subprocess.check_output(command, shell=True)
                else:                         # python 3.X
                    buf = subprocess.getoutput(command)
                gamma = float(buf.split()[1])
                dens = self.__open_field(directory+fluid+'dens'+str(on)+'.dat',dtype,fieldofview,slice)
                temp = self.__open_field(directory+'Temperature'+str(on)+'.dat',dtype,fieldofview,slice)
                self.data = temp*(dens**(1.0-gamma))
                self.strname += ' specific entropy'
            # ----
            # MASS ACCRETION RATE Mdot = abs(2pi R v_R Sigma)
            # ----
            if field == 'mdot':
                dens = self.__open_field(directory+fluid+'dens'+str(on)+'.dat',dtype,fieldofview,slice)
                if self.fargo3d == 'No':
                    vrad = self.__open_field(directory+fluid+'vrad'+str(on)+'.dat',dtype,fieldofview,slice)
                else:
                    vrad = self.__open_field(directory+fluid+'vy'+str(on)+'.dat',dtype,fieldofview,slice)
                for j in range(self.nsec):
                    for i in range(self.nrad):
                        self.data[i,j] = -2.0*np.pi*self.rmed[i]*vrad[i,j]*dens[i,j]
                self.strname += r' $\dot{M}$'

            #
            # ----
            # DISC ECCENTRICITY
            # ----
            if field == 'ecc':
                dens = self.__open_field(directory+fluid+'dens'+str(on)+'.dat',dtype,fieldofview,slice)
                if self.fargo3d == 'No':
                    vrad = self.__open_field(directory+fluid+'vrad'+str(on)+'.dat',dtype,fieldofview,slice)
                    vphi = self.__open_field(directory+fluid+'vtheta'+str(on)+'.dat',dtype,fieldofview,slice)
                else:
                    vrad = self.__open_field(directory+fluid+'vy'+str(on)+'.dat',dtype,fieldofview,slice)
                    vphi = self.__open_field(directory+fluid+'vx'+str(on)+'.dat',dtype,fieldofview,slice)

                # vphi is in the corotating frame!
                for i in range(self.nrad):
                    vphi[i,:] += (self.rmed)[i]*omegaframe
                    
                mass    = np.zeros((self.nrad,self.nsec))
                vrcent  = np.zeros((self.nrad,self.nsec))
                vtcent  = np.zeros((self.nrad,self.nsec))
                surface = np.zeros((self.nrad,self.nsec))
                rmed2D  = np.zeros((self.nrad,self.nsec))
                
                Rinf = self.redge[0:len(self.redge)-1]
                Rsup = self.redge[1:len(self.redge)]
                surf = np.pi * (Rsup*Rsup - Rinf*Rinf) / self.nsec
                for th in range(self.nsec):
                    surface[:,th] = surf
                    rmed2D[:,th]  = self.rmed

                # mass of each grid cell
                mass = dens*surface
                #print('disc total mass is: ', np.sum(mass))

                # loop below could probably be made more concise / optimized
                for i in range(self.nrad):
                    for j in range(self.nsec):
                        if i < self.nrad-1:
                            vrcent[i,j] = (self.rmed[i] - self.redge[i])*vrad[i+1,j] + (self.redge[i+1] - self.rmed[i])*vrad[i,j]
                            vrcent[i,j] /= (self.redge[i+1] - self.redge[i])
                        else:
                            vrcent[i,j] = vrad[i,j]
                        jm = j
                        jp = j+1
                        if (jp > self.nsec-1):
                            jp -= self.nsec
                        vtcent[i,j] = 0.5*(vphi[i,jm] + vphi[i,jp])

                Ar  = rmed2D*vtcent*vtcent/(1.0+mass)-1.0
                At = -rmed2D*vrcent*vtcent/(1.0+mass)
                self.data = np.sqrt(Ar*Ar + At*At)
                
                self.strname += ' eccentricity'
            #
            #
            # ----
            # VORTICITY or VORTENSITY
            # ----
            if (field == 'vorticity' or field == 'drl' or field == 'vortensity' or field == 'invvortensity'):

                if self.fargo3d == 'No':
                    vrad = self.__open_field(directory+fluid+'vrad'+str(on)+'.dat',dtype,fieldofview,slice)
                    vphi = self.__open_field(directory+fluid+'vtheta'+str(on)+'.dat',dtype,fieldofview,slice)
                else:
                    vrad = self.__open_field(directory+fluid+'vy'+str(on)+'.dat',dtype,fieldofview,slice)
                    vphi = self.__open_field(directory+fluid+'vx'+str(on)+'.dat',dtype,fieldofview,slice)
                
                # vphi is in the corotating frame!
                for i in range(self.nrad):
                    vphi[i,:] += (self.rmed)[i]*omegaframe
                # we first calculate drrvphi
                drrvphi = np.zeros((self.nrad,self.nsec))
                for j in range(self.nsec):
                    for i in range(1,self.nrad):
                        drrvphi[i,j] = ( (self.rmed)[i]*vphi[i,j] - (self.rmed)[i-1]*vphi[i-1,j] ) / ((self.rmed)[i] - (self.rmed)[i-1] )
                    drrvphi[0,j] = drrvphi[1,j]
                # then we calculate dphivr
                dphivr = np.zeros((self.nrad,self.nsec))
                for j in range(self.nsec):
                    if j==0:
                        jm1 = self.nsec-1
                    else:
                        jm1 = j-1
                    for i in range(self.nrad):
                        dphivr[i,j] = (vrad[i,j]-vrad[i,jm1])/2.0/np.pi*self.nsec
                # we deduce the vorticity or vortensity
                for j in range(self.nsec):
                    for i in range(self.nrad):
                        self.data[i,j] = (drrvphi[i,j] - dphivr[i,j]) / (self.redge)[i]

                # this is the radial derivative of the specific angular momentum
                if field == 'drl':
                    for j in range(self.nsec):
                        for i in range(self.nrad):
                            self.data[i,j] *= (self.redge)[i]
                    self.strname += r' $\partial_r \ell$'
                    
                # vortensity or inverse vortensity
                if field == 'vortensity' or field == 'invvortensity':
                    dens = self.__open_field(directory+fluid+'dens'+str(on)+'.dat',dtype,fieldofview,slice)
                    self.data /= dens
                    if field == 'invvortensity':
                        self.data = (1.0/self.data)
                        self.strname += ' inverse vortensity'
                    else:
                        self.strname += ' vortensity'
                if field == 'vorticity':
                    self.strname += ' vorticity'
            
            # ----
            # GAS self-gravitating accelerations
            # ----
            if (field == 'sgacctheta'):
                input_file = directory+'sgacctheta'+str(on)+'.dat'
                self.data = self.__open_field(input_file,dtype,fieldofview,slice)
                self.strname += r' SG $a_{\varphi}$'
                if physical_units == 'Yes' and nodiff == 'Yes':
                    self.unit = 1e-3*(self.culength)/(self.cutime)/(self.cutime)
                    self.strname += r' [km s$^{-2}$]'
            if (field == 'sgaccr'):
                input_file = directory+'sgaccr'+str(on)+'.dat'
                self.data = self.__open_field(input_file,dtype,fieldofview,slice)
                if par.log_xyplots_y == 'Yes':
                    self.data = np.abs(self.data)
                self.strname += r' SG $a_{r}$'
                if physical_units == 'Yes' and nodiff == 'Yes':
                    self.unit = 1e-3*(self.culength)/(self.cutime)/(self.cutime)
                    self.strname += r' [km s$^{-2}$]'


            # ----
            # Torque in every cell
            # ----
            if (field == 'torquesg'):
                input_file = directory+'torquesg'+str(on)+'.dat'
                self.data = self.__open_field(input_file,dtype,fieldofview,slice)
                self.strname += r' SG spec. torque'
                if physical_units == 'Yes' and nodiff == 'Yes':
                    self.unit = (self.culength)*(self.culength)/(self.cutime)
                    self.strname += r' [m$^{2}$ s$^{-1}$]'
            if (field == 'torquesumdisc'):
                input_file = directory+'torquesumdisc'+str(on)+'.dat'
                self.data = self.__open_field(input_file,dtype,fieldofview,slice)
                self.strname += r' spec. torque via summation'
                if physical_units == 'Yes' and nodiff == 'Yes':
                    self.unit = (self.culength)*(self.culength)/(self.cutime)
                    self.strname += r' [m$^{2}$ s$^{-1}$]'


            # ----
            # VERTICALLY-INTEGRATED (=surface) DENSITY
            # ----
            if self.fargo3d == 'Yes' and field == 'surfacedens':
                field = np.fromfile(directory+fluid+'dens'+str(on)+'.dat', dtype)
                datacube = field.reshape(self.ncol,self.nrad,self.nsec)
                datacube_cyl = np.zeros((self.nver,self.nrad,self.nsec))
                # sweep through the 3D cylindrical grid:
                for k in range(self.nver):
                    for i in range(self.nrad):
                        r     = np.sqrt( self.rmed[i]*self.rmed[i] + self.zmed[k]*self.zmed[k] )  # spherical radius
                        theta = math.atan( self.zmed[k]/self.rmed[i] )  # latitude ~ 0
                        isph = np.argmin(np.abs(self.rmed-r))
                        if r < self.rmed[isph] and isph > 0:
                            isph-=1
                        ksph = np.argmin(np.abs(self.tmed-theta))
                        if theta < self.tmed[ksph] and ksph > 0:
                            ksph-=1
                        if (isph < self.nrad-1 and ksph < self.ncol-1):
                            datacube_cyl[k,i,:] = ( datacube[ksph,isph,:]*(self.rmed[isph+1]-r)*(self.tmed[ksph+1]-theta) + datacube[ksph+1,isph,:]*(self.rmed[isph+1]-r)*(theta-self.tmed[ksph]) + datacube[ksph,isph+1,:]*(r-self.rmed[isph])*(self.tmed[ksph+1]-theta) + datacube[ksph+1,isph+1,:]*(r-self.rmed[isph])*(theta-self.tmed[ksph]) ) / ( (self.rmed[isph+1]-self.rmed[isph]) * (self.tmed[ksph+1]-self.tmed[ksph]) )
                        else:
                            # simple nearest-grid point interpolation...
                            datacube_cyl[k,i,:] = datacube[ksph,isph,:]
                # vertically-integrated density: \int_0^zmax rhoxdz
                # dz = Rxdtheta with dtheta uniform
                buf = np.sum(datacube_cyl,axis=0)  # (nrad,nsec)
                for i in range(self.nrad):
                    self.data[i,:] = buf[i,:]*np.abs(self.tmed[1]-self.tmed[0])*self.rmed[i]
                self.strname += ' surface density'
                if physical_units == 'Yes' and nodiff == 'Yes':
                    self.unit = (self.cumass*1e3)/((self.culength*1e2)**2.)
                    self.strname += r' [g cm$^{-2}$]'

            # ----
            # DUST STOKES NUMBER St = sqrt(pi/8) x (s rho_dust_int) / (H rho_gas)
            # ----
            if field == 'stokes':
                # gas mass surface (2D run) or volume density (3D run)
                if self.ncol > 1:  # 3D
                    field = np.fromfile(directory+'gasdens'+str(on)+'.dat', dtype)
                    rho_gas = field.reshape(self.ncol,self.nrad,self.nsec)
                else:
                    sigma_gas = self.__open_field(directory+'gasdens'+str(on)+'.dat',dtype,fieldofview,slice)
                
                # get dust internal density
                if self.fargo3d == 'Yes':
                    command = par.awk_command+' " /^DUSTINTERNALRHO/ " '+directory+'variables.par'
                else:
                    command = par.awk_command+' " /^Rhopart/ " '+directory+'*.par'

                if sys.version_info[0] < 3:   # python 2.X
                    buf = subprocess.check_output(command, shell=True)
                else:                         # python 3.X
                    buf = subprocess.getoutput(command)

                rho_dust_int = float(buf.split()[1])   # in g/cm^3
                rho_dust_int *= 1e3 # in kg/m^3
                rho_dust_int /= (self.cumass)
                rho_dust_int *= (self.culength**3.)  # in code units

                # get dust size
                if self.fargo3d == 'Yes':
                    dust_id, dust_size, dust_gas_ratio = np.loadtxt(directory+'/dustsizes.dat',unpack=True)
                    if fluid != 'gas':
                        s = dust_size[int(fluid[-1])-1] # fluid[-1] = index of dust fluid
                    else:
                        s = 1e-10  # arbitrarily small
                else:
                    command = par.awk_command+' " /^Sizepart/ " '+directory+'*.par'
                    if sys.version_info[0] < 3:   # python 2.X
                        buf = subprocess.check_output(command, shell=True)
                    else:                         # python 3.X
                        buf = subprocess.getoutput(command)
                    s = float(buf.split()[1])     # in meters

                s /= self.culength  # in code units

                # get angular frequency and sound speed, assuming
                # locally isothermal equation of state: first get the
                # aspect ratio and flaring index used in the numerical
                # simulation
                if self.ncol > 1:  # 3D
                    command = par.awk_command+' " /^ASPECTRATIO/ " '+directory+'*.par'
                    if sys.version_info[0] < 3:   # python 2.X
                        buf = subprocess.check_output(command, shell=True)
                    else:                         # python 3.X
                        buf = subprocess.getoutput(command)
                    aspectratio = float(buf.split()[1])                
                    # then get the flaring index used in the numerical simulation
                    command = par.awk_command+' " /^FLARINGINDEX/ " '+directory+'*.par'
                    if sys.version_info[0] < 3:
                        buf = subprocess.check_output(command, shell=True)
                    else:
                        buf = subprocess.getoutput(command)
                    flaringindex = float(buf.split()[1])
                    cs = np.zeros((self.ncol,self.nrad,self.nsec))
                    omega = np.zeros((self.ncol,self.nrad,self.nsec))
                    buf = np.zeros((self.ncol,self.nrad,self.nsec))
                    for i in range(self.nrad):
                        cs[:,i,:] = aspectratio*(((self.rmed)[i])**(-0.5+flaringindex))
                        omega[:,i,:] = self.rmed[i]**(-1.5)
                    buf = np.sqrt(np.pi/8.0) * (s*rho_dust_int) * omega / (cs*rho_gas) # 3D cube
                    if fieldofview == 'latitudinal' or fieldofview == 'vertical':
                        myfield = np.sum(buf,axis=2)/self.nsec  # azimuthally-averaged field (R vs. latitude)
                        self.data = myfield [::-1,:]
                    else:
                        self.data = buf[-1,:,:]   # midplane field
                else:  # 2D
                    self.data = 0.5*np.pi*s*rho_dust_int/sigma_gas
                self.strname += ' Stokes number'
                
            # ----
            # VRAD AND VTHETA for 2D CARTESIAN RUNS WITH FARGO3D
            # ----
            if self.fargo3d == 'Yes' and self.cartesian_grid == 'Yes':
                vx = self.__open_field(directory+fluid+'vx'+str(on)+'.dat',dtype,fieldofview,slice)
                vy = self.__open_field(directory+fluid+'vy'+str(on)+'.dat',dtype,fieldofview,slice)
                vrad_cart = np.zeros((self.nx,self.ny))
                vphi_cart = np.zeros((self.nx,self.ny))
                for i in range(self.nx-1):
                    for j in range(self.ny-1):
                        rmed = np.sqrt( self.xmed[i]*self.xmed[i] + self.ymed[j]*self.ymed[j] )
                        vxmed = 0.5*(vx[j,i]+vx[j,i+1])
                        vymed = 0.5*(vy[j,i]+vy[j+1,i])
                        #vxmed = vx[j,i]
                        #vymed = vy[j,i]
                        vrad_cart[i,j] = (self.xmed[i]*vxmed + self.ymed[j]*vymed)/rmed
                        vphi_cart[i,j] = (self.xmed[i]*vymed - self.ymed[j]*vxmed)/rmed
                vrad_cart[self.nx-1,:] = vrad_cart[self.nx-2,:]
                vphi_cart[self.nx-1,:] = vphi_cart[self.nx-2,:]
                vrad_cart[:,self.ny-1] = vrad_cart[:,self.ny-2]
                vphi_cart[:,self.ny-1] = vphi_cart[:,self.ny-2]
                if field == 'vrad':
                    #print(vrad_cart.min(),vrad_cart.max())
                    self.data = vrad_cart
                if field == 'vtheta':
                    #print(vphi_cart.min(),vphi_cart.max())
                    self.data = vphi_cart

            # ----
            # Non-axisymmetric part of gas density
            # ----        
            if field == 'naodens':
                dens = self.__open_field(directory+fluid+'dens'+str(on)+'.dat',dtype,fieldofview,slice)
                axidens = np.sum(dens,axis=1)/self.nsec
                self.data = dens-axidens.repeat(self.nsec).reshape(self.nrad,self.nsec)
                if par.verbose == 'Yes':
                    print('#### NAO DENSITY ###')
                    print(self.data.min(),self.data.max(),end='\r')
                self.strname = r'$\Sigma - \langle\Sigma\rangle_\varphi$'

            # ----
            # time-averaged particle density: means that we read
            # dustdensX.dat files for X from 0 to current output
            # number 'on' and we time-average arrays
            # ----        
            if fluid == 'pc' and field == 'rtadens':
                # case where on in paramsf2p.dat is specified as X,Y like for an animation: 
                # we then compute the r.t.a. density from on=X to on=Y
                if isinstance(par.on, int) == False:
                    on = range(par.on[0],par.on[1]+1,par.take_one_point_every)
                    for z in on:
                        print('reading pcdens'+str(z)+'.dat file',end='\r')
                        self.data += self.__open_field(directory+fluid+'dens'+str(z)+'.dat',dtype,fieldofview,slice)
                    self.data /= len(on)
                else:
                    for z in np.arange(on+1):
                        print('reading pcdens'+str(z)+'.dat file',end='\r')
                        self.data += self.__open_field(directory+fluid+'dens'+str(z)+'.dat',dtype,fieldofview,slice)
                    self.data /= len(np.arange(on))
                self.strname = 'r.t.a. particle density'
                if physical_units == 'Yes' and nodiff == 'Yes':
                    self.unit = (self.cumass*1e3)/((self.culength*1e2)**2.)
                    self.strname += r' [g cm$^{-2}$]'

            # ----
            # time-averaged particle density: means that we read
            # dustdensX.dat files for X from 0 to current output
            # number 'on' and we time-average arrays
            # ----        
            if field == 'densoveraxi':
                dens = self.__open_field(directory+fluid+'dens'+str(on)+'.dat',dtype,fieldofview,slice)
                axidens = np.sum(dens,axis=1)/self.nsec
                self.data = dens/(axidens.repeat(self.nsec).reshape(self.nrad,self.nsec))
                self.strname = r'$\Sigma / \langle\Sigma\rangle_\varphi$'
                
        # field name and units
        if field == 'dens':
            self.strname += ' density'
            #self.strname += r' $\Sigma$'
            if self.fargo3d == 'No':  # 2D
                if physical_units == 'Yes' and nodiff == 'Yes':
                    self.unit = (self.cumass*1e3)/((self.culength*1e2)**2.)
                    self.strname += r' [g cm$^{-2}$]'
            else:
                #self.strname += ' midplane density'
                if physical_units == 'Yes' and nodiff == 'Yes':
                    if self.ncol > 1:  # 3D
                        self.unit = (self.cumass*1e3)/((self.culength*1e2)**3.)
                        self.strname += r' [g cm$^{-3}$]'
                    else:  # 2D
                        self.unit = (self.cumass*1e3)/((self.culength*1e2)**2.)
                        self.strname += r' [g cm$^{-2}$]'
        if field == 'energy' and  self.fargo3d == 'Yes':
            self.strname += ' sound speed'
            if physical_units == 'Yes' and nodiff == 'Yes':
                self.unit = 1e-3*(self.culength)/(self.cutime)
                self.strname += r' [km s$^{-1}$]'
        if field == 'vrad' or field == 'vy':
            if field == 'vy' and self.cartesian_grid == 'Yes':
                self.strname += r' $v_{y}$'
            else:
                self.strname += r' $v_{r}$'
            if physical_units == 'Yes' and nodiff == 'Yes':
                self.unit = 1e-3*(self.culength)/(self.cutime)
                self.strname += r' [km s$^{-1}$]'
        if field == 'vtheta' or field == 'vx':
            if field == 'vx' and self.cartesian_grid == 'Yes':
                self.strname += r' $v_{x}$'
            else:
                self.strname += r' $v_{\varphi}$'
            if physical_units == 'Yes' and nodiff == 'Yes':
                self.unit = 1e-3*(self.culength)/(self.cutime)
                self.strname += r' [km s$^{-1}$]'
        if field == 'vcol' or field == 'vz':
            self.strname += r' $v_{\theta}$'
            if physical_units == 'Yes' and nodiff == 'Yes':
                self.unit = 1e-3*(self.culength)/(self.cutime)
                self.strname += r' [km s$^{-1}$]'
        if field == 'temp':
            self.strname += ' temperature'
            if physical_units == 'Yes' and nodiff == 'Yes':
                self.unit = self.cutemp
                self.strname += ' [K]'
        if field == 'mdot':
            if physical_units == 'Yes' and nodiff == 'Yes':
                self.unit = (self.cumass)/(self.cutime)
                self.unit /= 2e30  
                self.unit *= 3.15e7
                self.strname += r' [$M_{\odot}$ yr$^{-1}$]'

        #
        if physical_units == 'No' and nodiff == 'Yes':
            self.strname += r' [code units]'            

        if onedprofile == 'No':
            self.strname += ' at '+self.strtime
       

    def __open_field(self, f, dtype, fieldofview, slice):
        """
        Reading the data
        """
        field = np.fromfile(f, dtype=dtype)

        # cuidadin: replace NaNs by 0 
        where_are_NaNs = np.isnan(field)
        field[where_are_NaNs] = 0.0
        
        if self.ncol == 1: # 2D
            return field.reshape(self.nrad,self.nsec)
        else: # 3D
            datacube = field.reshape(self.ncol,self.nrad,self.nsec)
            if fieldofview == 'latitudinal' or fieldofview == 'vertical':
                if slice == 'average':
                    myfield = np.sum(datacube,axis=2)/self.nsec  # azimuthally-averaged field (R vs. latitude)
                    return np.transpose(myfield [::-1,:])
                else:
                    return np.transpose(datacube[:,:,self.nsec//2])  # azimuthal cut at planet's location
            else:  # polar or cartesian fields of view
                if np.abs(self.zmax-1.57) < 0.01:
                    if slice == 'midplane':
                        return datacube[-1,:,:]   # midplane field only if "half-a-disc" is simulated in latitudinal direction!
                    if (slice == 'upper' or slice == '#'):
                        return datacube[1,:,:]   # field at disc surface only if "half-a-disc" is simulated in latitudinal direction!
                    if slice == 'intermediate':
                        return datacube[self.ncol//2,:,:]   # field at disc surface only if "half-a-disc" is simulated in latitudinal direction!
                else:
                    if slice == 'midplane':
                        return datacube[self.ncol//2,:,:]   # midplane field only if "full" disc is simulated in latitudinal direction!
                    if slice == 'lower':
                        return datacube[1,:,:]   # field at lower surface only if "half-a-disc" is simulated in latitudinal direction!
                    if slice == 'upper':
                        return datacube[-1,:,:]   # field at upper surface only if "half-a-disc" is simulated in latitudinal direction!


    def compute_streamline(self, dtype='float64',niterations=100000,R0=0,T0=0,rmin=0,rmax=1e4,pmin=0,pmax=6.28,forward=True,fieldofview='polar',slice='midplane'):
        
        # first import global variables
        import par
        
        mydirectory = par.directory
        if mydirectory[-1] != '/':
            mydirectory += '/'

        # several output numbers
        if isinstance(par.on, int) == False:
            on = par.on[0]
        else:
            on = par.on
            
        # vrad should be subtracted by the instantaneous planet's migration rate da/dt ...
        
        if self.fargo3d == 'No':
            vrad = self.__open_field(mydirectory+par.fluid+'vrad'+str(on)+'.dat',dtype,fieldofview,slice)
            vrad = np.roll(vrad, shift=int(self.nsec/2), axis=1)
            vphi = self.__open_field(mydirectory+par.fluid+'vtheta'+str(on)+'.dat',dtype,fieldofview,slice)
            vphi = np.roll(vphi, shift=int(self.nsec/2), axis=1)
        else:
            vrad = self.__open_field(mydirectory+par.fluid+'vy'+str(on)+'.dat',dtype,fieldofview,slice)
            vphi = self.__open_field(mydirectory+par.fluid+'vx'+str(on)+'.dat',dtype,fieldofview,slice)
           
        # forward or backward integration of streamlines
        if forward == False:
            vrad = -vrad
            vphi = -vphi

        # Get radius and azimuth (code units)
        R = self.redge
        T = self.pmed

        if par.physical_units == 'Yes': # back to code units...
            R0 /= (self.culength/1.5e11)  
            R = self.redge/((self.culength/1.5e11))
        
        # Euler-first order integration step
        stepstrl = (pmax-pmin)/self.nsec  # could be refined?
        
        # XS and YS are lists with the X- and Y-coordinates of points along a given streamline
        XS = []; YS = []
        myR = R0; myT = T0
        nloop = 0

        while( myR > R.min() and myR < R.max() and myT > T.min() and myT < T.max() and nloop < niterations ):
            imin = np.argmin(np.abs(R-myR))
            jmin = np.argmin(np.abs(T-myT))
            if jmin >= self.nsec:
                jmin = self.nsec-1
            if jmin < 0:
                jmin = 0
            vrnew = vrad[imin,jmin]  # a bilinear interpolate would be better...
            vpnew = vphi[imin,jmin]  # a bilinear interpolate would be better...
            dr = stepstrl * myR*vrnew / np.sqrt ( vpnew*vpnew + myR*myR*vrnew*vrnew )
            dt = stepstrl * vpnew     / np.sqrt ( vpnew*vpnew + myR*myR*vrnew*vrnew )
            myR = myR + dr
            myT = myT + dt
            if par.fieldofview == 'polar' and par.rvsphi == 'No':
                if par.physical_units == 'Yes':
                    XS.append(myR*self.culength/1.5e11)
                else:
                    XS.append(myR)
                YS.append(myT)
            if par.fieldofview == 'polar' and par.rvsphi == 'Yes':
                XS.append(myT)
                if par.physical_units == 'Yes':
                    YS.append(myR*self.culength/1.5e11)
                else:
                    YS.append(myR)
            if par.fieldofview == 'cart':
                if par.physical_units == 'Yes':
                    # this perhaps requires editing if grid's azimuthal extent is smaller than 2pi...
                    XS.append(myR*np.cos(myT+np.pi)*self.culength/1.5e11)  # add pi since field is shifted by nsec/2 cells in azimuth
                    YS.append(myR*np.sin(myT+np.pi)*self.culength/1.5e11)
                else:
                    XS.append(myR*np.cos(myT+np.pi))  # add pi since field is shifted by nsec/2 cells in azimuth
                    YS.append(myR*np.sin(myT+np.pi))   
            nloop += 1

        return XS,YS
