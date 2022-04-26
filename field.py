import numpy as np
import os
import sys
import subprocess
import math
import itertools

# first import global variables
import par

# inherit from Mesh class
from mesh import *

class Field(Mesh):
    def __init__(self, field, directory='', dtype='float64'):
        if len(directory) > 0:
            if directory[-1] != '/':
                directory += '/'

        self.directory = directory
        self.verbose = par.verbose
        
        # read valp/vecp header
        self.read_valp_header()
        
        # get all Mesh attributes
        if par.plot_qz_eigenfq == 'No':
            Mesh.__init__(self, directory)

        # work out global string with main run parameters to be
        # displayed in the bottom part of the plots
        self.globstr  = r'$N_r$='+str(self.nr)
        self.globstr += r'$\;\;L$='+str(self.nth-1)
        if 'm' in dir(self):
            self.globstr += r'$\;\;m$='+str(self.m)
        if 'e' in dir(self):
            self.globstr += r'$\;\;E$='+par.str_fmt_ek(self.e)
        if 'eta' in dir(self):
            self.globstr += r'$\;\;\eta$='+str(self.eta)
        if 'omslope' in dir(self) and self.omslope != 0.0:
            self.globstr += r'$\;\;\sigma$='+str(self.omslope)
        if 'epsilon' in dir(self) and self.epsilon != 0.0:
            self.globstr += r'$\;\;\epsilon$='+str(self.epsilon)
        if 'n' in dir(self) and self.n != 0.0:
            self.globstr += r'$\;\;N$='+str(self.n)
        if 'pr' in dir(self) and self.pr != 0.0:
            self.globstr += r'$\;\;Pr$='+str(self.pr)

        # special case where a radial BVF profile is used, in which
        # case we read the file and assume N is no longer a constant
        # but a 1D array:
        if os.path.isfile(directory+'/BVFprofile'):
            print('BVFprofile file found: I understand that you have assumed a radially varying Brunt-Vaisala frequency')
            self.n = np.zeros(self.nr+1, dtype=float)
            with open(directory+'/BVFprofile', 'r') as f:
                for i in range(self.nr+1):
                    self.n[i] = f.readline().split()[1]
                    
        # case a meridional cut is requested
        if (par.plot_zcut == 'Yes' or par.plot_zcut_3D == 'Yes'):
            if self.verbose == 'Yes':
                print('starting calculation of meridional cuts...')
            if (par.field == 'ek' or par.field == 'dissv' or par.field == 'shear'):
                (self.data,self.field_names) = self.read_field('dissec')
            if (par.field == 'et' or par.field == 'disst'):
                (self.data,self.field_names) = self.read_field('temper')
            if self.verbose == 'Yes':
                print('done!')
             
        # case we want to display the spectral content of the modes
        # velocity field
        if par.plot_spectrum == 'Yes':
            if self.verbose == 'Yes':
                print('now plotting modes spectral content...')
            (self.ck,self.cl) = self.read_spectra()
            if self.verbose == 'Yes':
                print('done!')

                
    def read_valp_header(self):
        
        # read valp first time to get number of eigenmodes
        with open(self.directory+'valp', 'r') as f:
            nbvp = len(f.readlines())-2
            print('number of eigenmodes in valp: ', nbvp)

        # read it a second time
        with open(self.directory+'valp', 'r') as f:
            i = 2
            line = f.readline()
            if 'cyl' in line:
                self.rotation = 'cylindrical'
            if 'coni' in line:
                self.rotation = 'conical'
            if line.startswith(' st'):
                line = f.readline()
                i += 1
            while line.find('nr=') == -1:
                line = f.readline()
                i += 1

            line = line[1:]
            res = line.replace('=',' ').split()
            self.nth = 2
            self.lmin = 0
            self.omslope = 0.0
            self.epsilon = 0.0
            self.gamma = 0.0
            self.eta = 0.0
            self.n = 0.0
            self.pr = 0.0
            for (j, value) in enumerate(res):
                if not j%2:
                    name = value.lower()
                else:
                    try:
                        nvalue = int(value)
                    except:
                        nvalue = float(value.replace('D','e'))
                    setattr(self, name, nvalue)
                    # secondary parameters
                    if name == 'lmax':
                        self.nth = min(self.lmax+1, 5000)
                    elif name == 'nnmax':
                        self.nth = min(self.nnmax+1, 5000)
                        self.lmin = 1
            if self.lmin == 0:
                self.lmin = max(0, self.m)
                
            # we allocate array of eigenfrequencies:
            nc = 0
            self.omr = np.zeros(nbvp, dtype=float)
            self.omi = np.zeros_like(self.omr)
            for line in f:
                res = line.replace('=',' ').split()
                self.omr[nc] = res[0]
                self.omi[nc] = res[1]    
                nc += 1
 
            self.nmodes = nc

            # infer type of differential rotation based on reading
            # parameters omslope and epsilon via vecp header:
            if self.omslope == 0.0 and self.epsilon == 0.0:
                self.rotation = 'solid'    
            if self.omslope != 0.0:
                self.rotation = 'shellular'

            if self.verbose == 'Yes':
                print('number of eigenmodes found: ',self.nmodes)
                print('eigenfrequencies read in valp:')
                for i in range(self.nmodes):
                    print(self.omi[i],self.omr[i])
                print('rotation = ', self.rotation)

                
    def read_array(self, f, n1, n2):
        i1 = 0
        i2 = 0
        res = np.zeros((n1, n2), dtype=float)
        while i2 < n2:
            line = f.readline()
            for number in line.split():
                res[i1, i2] = float(number)
                i1 += 1
                if i1 == n1:
                    i2 += 1
                    i1 = 0
        return res

    
    def read_field(self, filename):
        with open(self.directory + filename, 'r') as f:
            f.readline()
            nc = int(f.readline())
            fields = np.zeros((self.nr+1, self.nth, nc, self.nmodes), dtype=float)
            field_names = []
            for im in range(self.nmodes):
                f.readline()
                for ic in range(nc):
                    field_names.append(f.readline().strip().replace('!d','_{').replace('!n', '}'))
                    fields[:,:, ic, im] = self.read_array(f, self.nr+1, self.nth)
        return (fields,field_names)

    
    def read_scalar(self, f, n):
        i = 0
        res = np.zeros(n, dtype=float)
        while i < n:
            values = f.readline().split()
            for value in values:
                res[i] = float(value)
                i += 1
        return res

                
    def read_spectra(self):
        with open(self.directory+'spectre_2D', 'r') as f:
            f.readline()
            nc = int(f.readline())
            
            ck = np.zeros((self.nr+1, self.ndomains, nc, self.nmodes), dtype=float)
            cl = np.zeros((self.lmax - self.lmin + 1, self.ndomains, nc, self.nmodes), dtype=float)
            self.lmin = np.zeros(nc, dtype=int)

            # loop over modes
            for k in range(self.nmodes):
                f.readline()  # (omr,omi)
                # loop on variables
                for n in range(nc):
                    f.readline()  # variable name
                    # loop on domains
                    for domain in range(self.ndomains):
                        nzd = self.nzd[domain]
                        ck[:nzd, domain, n, k] = self.read_scalar(f, nzd)

                    if self.lmax != 0:
                        self.lmin[n] = int(f.readline())
                        nl = (self.lmax - self.lmin[n]) // 2 + 1
                        for domain in range(self.ndomains):
                            cl[:nl, domain, n, k] = self.read_scalar(f, nl)

            return (ck,cl)
    

    def compute_characteristics(self,omega,niterations=10000,sstart=0.4,zstart=0.4):
        
        # first import global variables
        import par
        
        # SC and ZC are lists of the s- and z-coordinates along paths
        # of characteristic
        SC = []; ZC = []
        sarray = np.zeros(int(niterations))
        zarray = np.zeros(int(niterations))

        if par.caract_s_z[0] != '##':
            sstart = par.caract_s_z[0]
        else:
            print('starting s-coordinate for paths of characteristics is undefined, I am setting it to 0.5.')
            sstart = 0.5
        if par.caract_s_z[1] != '##':
            zstart = par.caract_s_z[1]
        else:
            print('starting z-coordinate for paths of characteristics is undefined, I am setting it to 0.5.')
            zstart = 0.5
        s = sstart
        z = zstart
        
        eps = 1e-3      # initial integration step
        ds = eps
        dz = eps
        slope = 1.0     # +1 or -1 (changes sign if a reflexion occurs)
        step_in_ds = False
        step_in_dz = False
        
        # mode's frequency in inertial frame
        omegap = omega  

        # check that initial condition is in hyperbolique domain (xi > 0)
        (xi,dzds0,dsdz0,buf) = self.compute_dzds_dsdz_caract(omegap,slope,s,z)
        if (xi < 0.0):
            sys.exit('wrong choice of initial condition for the calculation of characteristics, as xi < 0: ', xi)

        # ==============================
        # we integrate the equation of characteristics by iteration,
        # starting from s=sstart and z=zstart. We integrate either z
        # from s or s from z depending on the absolute value of the
        # slopes dz/ds and ds/dz (as in Mirouh + 2016):
        for k in range(int(niterations)):
        # ==============================
            sarray[k] = s
            zarray[k] = z

            if k != 0:
                (xi_prev,dzds_prev,dsdz_prev,buf) = self.compute_dzds_dsdz_caract(omegap,slope,sarray[k-1],zarray[k-1])
            (xi,dzds,dsdz,buf) = self.compute_dzds_dsdz_caract(omegap,slope,s,z)

            # this 'if condition' is required to avoid newz be a NaN in
            # case xi would be < 0:
            if (xi >= 0.0):
                if np.abs(dzds) < np.abs(dsdz):
                    # case we were integrating wrt z at previous iterations:
                    if step_in_dz == True:
                        if dz*dzds*ds < 0.0:
                            ds = -ds
                    # step in ds
                    news = s+ds
                    newz = z+dzds*ds
                    step_in_ds = True
                    step_in_dz = False
                else:
                    # case we were integrating wrt s at previous iterations:
                    if step_in_ds == True:
                        if dz*dsdz*ds < 0.0:
                            dz = -dz
                    # step in dz
                    newz = z+dz
                    news = s+dsdz*dz
                    step_in_ds = False
                    step_in_dz = True
            else: # if xi < 0 we get back to s and z at beginning of previous iteration:
                news = sarray[k-1] # or s?
                newz = zarray[k-1] # or z?
            newr = np.sqrt(news*news+newz*newz)

            # ------------------------------
            # case where a reflection occurs
            # ------------------------------
            if (newz < 0 or news < 0 or newr < self.r.min() or newr > self.r.max() or xi < 0.0):

                # slope changes sign whatever kind of reflexion:
                slope = -slope
                
                # reflexion along equatorial axis: we need only change
                # sign of dz if we integrate wrt z
                if (newz < 0 and step_in_dz == True):
                    dz = -dz

                # reflexion along rotation axis: we need only change
                # sign of ds if we integrate wrt s
                if (news < 0 and step_in_ds == True):
                    ds = -ds
                        
                # reflexion along inner or outer radius: dr needs to
                # change sign. Since rdr = sds+zdz, it means that
                # dsx(s + zdz/ds) needs to change sign if we integrate
                # wrt s, or dzx(z + sds/dz) needs to change sign if we
                # integrate wrt z
                if (newr < self.r.min() or newr > self.r.max()):
                    (bufxi,bufdzds,bufdsdz,buf) = self.compute_dzds_dsdz_caract(omegap,slope,s,z)
                    if step_in_ds == True:
                        bufs = s+ds
                        bufz = z+bufdzds*ds
                        if (s+z*dzds)*(bufs+bufz*bufdzds) > 0.0:
                            ds = -ds
                    if step_in_dz == True:
                        bufz = z+dz
                        bufs = s+bufdsdz*dz
                        if (z+s*dsdz)*(bufz+bufs*bufdsdz) > 0.0:
                            dz = -dz

                # reflexion along xi=0 surface: 'dxi' needs to change sign
                if (xi < 0.0):
                    if step_in_ds == True:
                        ds = -ds
                        s = news+ds
                        z = newz+dzds_prev*ds
                    if step_in_dz == True:
                        dz = -dz
                        z = newz+dz
                        s = news+dsdz_prev*dz                    

                (xi,dzds,dsdz,buf) = self.compute_dzds_dsdz_caract(omegap,slope,s,z)
                if step_in_ds == True:
                    news = s+ds
                    newz = z+dzds*ds
                if step_in_dz == True:
                    newz = z+dz
                    news = s+dsdz*dz
                
            # ------- end various cases of reflexions
            
            # We finally update s and z:
            s = news
            z = newz
            
            if par.last_only == 'Yes':
                if k > int(0.8*niterations):
                    SC.append(s)
                    ZC.append(z)
            else:
                SC.append(s)
                ZC.append(z)

            #k+=1
                
        return SC,ZC


    # ====================================
    # function that computes the slope dz/ds for solving the equation
    # of characteristics, which is valid for gravito-inertial modes
    # with differential rotation and a constant BVF.
    #
    # input parameters:
    # omegap = mode's frequency in inertial frame
    # N = Brunt-Vaisala frequency (constant so far)
    # slope = +1 or -1
    # s,z: local coordinates inside the shell
    def compute_dzds_dsdz_caract(self,omegap,slope,s,z):
    # ====================================
        # we consider the different cases for differential rotation: solid-body
        # rotation, differential rotation either shellular, cylindrical or
        # conical

        # case 1: solid-body rotation, for which Omega = 1 throughout the shell
        if self.rotation == 'solid':
            Omega = 1.0
            Az = 0.0
            As = 4.0
            
        # case 2: shellular differential rotation, with Omega(r) =
        # Omega(Rmax) x (r/Rmax)^omslope, Rmax=1, Omega(Rmax)=1 in our
        # set of code units. Parameter omslope is inhereted via the
        # Field class
        if self.rotation == 'shellular':
            Omega = (s*s + z*z)**(0.5*self.omslope)
            Az = 2.0*(Omega**2.0)*s*z*self.omslope/(s*s+z*z)
            As = 4.0*(Omega**2.0)*(1.0 + 0.5*self.omslope*s*s/(s*s+z*z))

        # case 3: cylindrical differential rotation, with Omega(s) = 1
        # + epsilon*s*s in our set of code units. Parameter epsilon is
        # also inhereted via the Field class
        if self.rotation == 'cylindrical':
            Omega = 1.0 + self.epsilon*s*s
            Az = 0.0
            As = 4.0*(Omega**2.0)*(1.0 + self.epsilon*s*s/(1.0+self.epsilon*s*s))

        # case 4: conical differential rotation, with Omega(s) = 1 +
        # epsilon*sin^2(theta) in our set of code units. Parameter
        # epsilon is also inhereted via the Field class
        if self.rotation == 'conical':
            Omega = 1.0 + self.epsilon*s*s/(s*s+z*z)
            Az = -4.0*Omega*self.epsilon*(s**3.0)*z*((s*s+z*z)**(-2.0))
            As = 4.0*(Omega**2.0)*(1.0 + self.epsilon*s*s*z*z*((s*s+z*z)**(-2.0))/4.0/Omega/Omega)
        
        # case where mode has been computed with one of Michel's eq files!
        if par.eq_file_michel == 'Yes':
            omegap *= 2.0
            omegap -= self.m*Omega
            
        # Doppler-shifted frequency and its square:
        Omegatilde = omegap + self.m*Omega
        omega2 = Omegatilde**2.0

        # Brunt-Vaisala frequency N:
        # case where N is constant throughout the domain
        if isinstance(self.n, (list, tuple, np.ndarray)) == False:
            N = self.n        
        else:
        # otherwise a radial profile for N has been used
            r = np.sqrt(s*s+z*z)
            if isinstance(r, (list, tuple, np.ndarray)) == True:
                N = np.zeros_like(r)
                for i in range(self.nr+1):
                    for j in range(self.nth):
                        index = np.argmin(np.abs(self.r-r[i,j]))
                        N[i,j] = self.n[index]
            else:
                index = np.argmin(np.abs(self.r-r))
                N = self.n[index]
        
        # xi-parameter (discriminant of the reduced equation), see
        # expression in Mirouh+ 2016 (JFM):        
        xi = Az*(Az/4.0 + 4.0*N*N*s*z) - As*(4.0*N*N*z*z - omega2) - omega2*(omega2 - 4.0*N*N*(s*s+z*z))
        
        # slopes dz/ds and ds/dz of the caracteristics:
        num = 4.0*N*N*s*z + 0.5*Az + slope*np.sqrt(xi)
        den = omega2 - 4.0*N*N*z*z
        dzds = num/den

        num = 4.0*N*N*s*z + 0.5*Az - slope*np.sqrt(xi)
        den = omega2 - As - 4.0*N*N*s*s
        dsdz = num/den

        omegatilde_crit = Omegatilde - 2.0*N*z
        return (xi,dzds,dsdz,omegatilde_crit) # instead of omegatilde as last argument
    

    # ====================================
    # function that computes the critical latitudes at the inner and outer
    # radial edges of the shell
    def compute_critical_latitudes(self,omega):
    # ====================================

        # We compute Omega at inner and outer radii for solid-body 
        # or shellular differential rotation
        if self.rotation == 'solid':
            omega_in  = 1.0
            omega_out = 1.0

        if self.rotation == 'shellular':
            omega_in = self.eta**(self.omslope)
            omega_out = 1.0

        # case where mode has been computed with one of Michel's eq files!
        if par.eq_file_michel == 'Yes':
            omega *= 2.0
            omega -= self.m
            
        # Doppler-shifted frequencies at inner and outer radii:
        omegatilde_in  = omega + self.m*omega_in
        omegatilde_out = omega + self.m*omega_out

        sintheta_in  = 0.5*omegatilde_in/omega_in
        sintheta_out = 0.5*omegatilde_out/omega_out
        return (sintheta_in,sintheta_out)
