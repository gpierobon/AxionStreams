import os, sys
import h5py as h5
import numpy as np
import scipy as scp

path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(path)

from galpy import potential
from astropy import units as u
from AxionStreams import density 
from AxionStreams import perturbations 
from AxionStreams import sample
from AxionStreams import model

class Stream():
    def __init__(self,stream_id,verb=False):
        self.ID = stream_id 
        if verb:
            print("Created stream with ID", stream_id)
    
    def print0(self,message):
        if self.ID == 0:
            print(message)

    def run(self, fin, fout, ts):
        '''
        '''
        self.print0("Reading orbits from %s ... "%fin.split('/')[-1])

        self.read_orbit(fin)
        self.set_parameters() 
        self.set_encounters(ts)
        self.set_perturb_data()
        self.set_perturb_array()
        self.set_stream_energies()
        self.perturb(saveE=False,debug=False)
        self.get_local_M(ts,debug=False)
        self.dump(fout,hdf5=False)

    # --------------------------------------------------

    def read_orbit(self, fname):
        '''
        '''
        with h5.File(fname,'r') as f:
            #self.e = f['Orbit_%03d/Eccentricity'%self.ID]  
            #self.R = np.array(f['Orbit_%03d/R'%self.ID])
            self.x = np.array(f['Orbit_%03d/x'%self.ID])
            self.y = np.array(f['Orbit_%03d/y'%self.ID])
            self.z = np.array(f['Orbit_%03d/z'%self.ID])
            self.vx = np.array(f['Orbit_%03d/vx'%self.ID])
            self.vy = np.array(f['Orbit_%03d/vy'%self.ID])
            self.vz = np.array(f['Orbit_%03d/vz'%self.ID])
    
    def set_parameters(self):
        '''
        Units: Msun, kpc, km/s
        '''
        self.isstream = 0
        self.G_N = 4.3e-6        # kpc*km^2/Msun/s^2
        self.c2 = (3e5)**2
        
        pl = 0  # To implement better
        ic = 0  # To implement better

        if pl == 0:
            self.prof = 'NFW'
        else:
            self.prof = 'PL'
        
        self.alpha2 = model.get_alpha2(pl)
        self.beta = model.get_beta(pl)
        self.kappa = model.get_kappa(pl)

        # Physical parameters 
        self.Mass = sample.draw_random_Mass(ic)
        self.IMass = self.Mass      
        self.Rad = model.get_radius(ic,self.Mass)
        #print(self.IMass,self.Rad)
        
        self.bmax = 2e-5         
        self.bmin = 1e-5         

        self.vin = np.sqrt(self.vx[0]**2+self.vx[1]**2+self.vx[2]**2)    # km/s
        self.vdisp = (self.kappa*self.G_N*self.Mass/self.Rad)**(1/2)     # km/s
        
        # To implement
        # Transition radius, Eq. (7) of 2207.11276, f_b=6
        #self.bs = 6*np.sqrt(2*np.sqrt(self.alpha2)/(3*self.beta))*self.Rad 
        

    def set_encounters(self,ts):
        '''
        '''
        kms_to_kpcGyr = 1.02269 
        star_density = density.stellar_density(self.x,self.y,self.z)
        star_density /= np.sum(star_density)

        vel = np.abs(np.sqrt(self.vx**2+self.vy**2+self.vz**2))*kms_to_kpcGyr
        inte = 1e9*star_density*vel*np.pi*self.bmax**2 # Units are kpc,Gyr
        finterp = scp.interpolate.InterpolatedUnivariateSpline(ts, inte, k=1)
        xx = np.linspace(ts[0], ts[-1], 5*len(ts))
        qq = [finterp.integral(0, t) for t in xx]        
        
        self.N_encounters  = int(round(qq[-1]))

        if self.N_encounters > 1e6:
            self.N_encounters = 1e6
        
        new_enc,t_enc_ind,enc_deg = density.get_ts_encounters(ts,star_density,self.N_encounters)
        self.N_encounters = new_enc
        self.t_enc_ind = t_enc_ind
        self.enc_deg = enc_deg
        self.print0("Stream %d has %d encounters"%(self.ID,self.N_encounters))                    
        
        del xx, qq, star_density

    def set_stream_energies(self):
        '''
        '''
        self.Ebind = self.beta*self.G_N*self.Mass**2/self.Rad*self.c2
        self.Etot = (0.5*self.kappa/self.beta - 1)*self.Ebind*self.c2

    def set_perturb_data(self):
        '''
        '''
        pdata = perturbations.Perturb(self.prof)
        pdata.load_perturb_data()
        self.dM_I = pdata.dM_interp
        self.fej_I = pdata.fej_interp
        self.fub_I = pdata.fub_interp
    
    def set_perturb_array(self):
        '''
        '''
        # Fix the relative velocity in the encounters
        kpc = u.kpc
        kms = u.km/u.s
        R_list = [np.sqrt(self.x[ti]**2+self.y[ti]**2) for ti in self.t_enc_ind]
        phi_list = [np.arctan2(-self.y[ti], self.x[ti]) for ti in self.t_enc_ind]
        vcirc_list = [potential.vcirc(potential.MWPotential2014, R=Ri,vo=220*kms, ro=8*kpc) for Ri in R_list]
        
        vstar_x = [-vcirc_list[i]*np.sin(phi_list[i]) for i in range(self.N_encounters)]
        vstar_y = [ vcirc_list[i]*np.cos(phi_list[i]) for i in range(self.N_encounters)]       
        venc_x = [self.vx[ti] - vstar_x[i] for ti,i in zip(self.t_enc_ind,range(self.N_encounters))]
        venc_y = [self.vx[ti] - vstar_x[i] for ti,i in zip(self.t_enc_ind,range(self.N_encounters))]

        self.vlist = [np.sqrt(venc_x[i]**2+venc_y[i]**2+self.vz[ti]**2) for ti,i in zip(self.t_enc_ind,range(self.N_encounters))] 
        self.vlist = np.array(self.vlist)
        
        del vstar_x, vstar_y, venc_x, venc_y, vcirc_list, R_list, phi_list

        # Draw a random impact parameter
        self.blist = sample.draw_random_b(self.bmin,self.bmax,self.N_encounters) # kpc
        
        # Delta E imparted, there is an implicit Mstar**2=1 Msun**2 in the numerator to fix units
        self.DelE = 4./3.*self.G_N**2/(self.blist**4*self.vlist**2)*self.Mass*self.alpha2*self.Rad**2*self.c2

    def perturb(self,saveE=False,debug=False):
        '''
        '''
        counter = 0; dE_max = 0; 
        totEi = 0; fs = 0
        self.t_disr_ind = -1

        Mloss_list = []
        vdisp_list = []
        if saveE == True:
            Eb_list = []; R_list = []; M_list = []
        
        # Encounter loop
        for i in range(self.N_encounters):
            
            # Fractional energy imparted, with encounter degeneracy
            dE = self.DelE[i]/self.Ebind*self.enc_deg 
            totEi += self.DelE[i]*self.enc_deg      
            
            # Identify dE max
            if dE > dE_max:
                dE_max = dE 
            
            if (dE < 1e-4):
                dM = 0.0
                dE_remain = dE*self.Ebind
            else:
                dM = self.dM_I(dE)*self.Mass
                dE_remain = dE*(1 - self.fej_I(dE)) - self.fub_I(dE)*self.Etot
                dE_remain *= self.Ebind
           
            # Mass loss modifieds mass
            newMass = self.Mass - dM
            newRad = self.Rad
            newVel = self.vdisp*(1+dE)
            
            # Energy conservation, needed to establish the new radius 
            E_f = self.Etot + dE_remain
           
            # Here we check if MC is fully disrupted, and set disruption time in case
            if totEi > 0.9995*self.Ebind:
                if fs == 0:
                    self.isstream = 1
                    self.t_disr_ind = self.t_enc_ind[counter]
                    fs = 1    
            else:
               # We updated radius as long as we have a minicluster
               newRad = (0.5*self.kappa - self.beta)*self.G_N*newMass**2/E_f*self.c2**2

            if debug == True:
                self.print0("dM %g dE_remain %g Etot %g E_f %g"%(dM,dE_remain,self.Etot,E_f))

            # Update mass, radius, vel. disp.
            self.Rad = newRad
            self.Mass = newMass
            self.vdisp = newVel
            
            # Append to save for the analysis
            Mloss_list.append(dM)
            vdisp_list.append(newVel)

            if saveE == True:
                Eb_list.append(self.Ebind)
                M_list.append(self.Mass)
                R_list.append(self.Rad)
            
            # Recompute energies 
            self.set_stream_energies()
            counter += 1
            

        # To Numpy array for the analysis
        self.MassLoss = np.array(Mloss_list)
        self.Vdisp = np.array(vdisp_list)
        if saveE == True:
            self.Ebl = np.array(Eb_list)
            self.Rl = np.array(R_list)
            self.Ml = np.array(M_list)
        
        # Example print for first stream
        self.print0("For stream %d max dE is %.8f"%(self.ID,dE_max))

    def get_local_M(self,ts,debug=False):
        '''
        '''
        dL = 1e-9               # 1 mpc in kpc
        mloc = 0.0
        Mremain = self.Mass     
        
        for i in range(self.N_encounters):
            # Exclude i = -1, which is the last encounter
            if self.t_enc_ind[i] in [-3,-2,-1]:
                continue
            mloc += self.MassLoss[i]*dL/(self.Vdisp[i]*(ts[-1]-ts[self.t_enc_ind[i]]))
        
        # Convert to Msun and add remaining MC mass
        phys_conv = 3.0857e16/3.1557e16 
        self.Slocal = mloc*phys_conv
        self.Mlocal = self.Slocal + Mremain
        
        self.Mlocal2 = mloc*phys_conv
        if debug == True:
            print(self.IMass,self.Mass,self.Mlocal,self.Slocal)

    def dump(self,fileout,hdf5=False):
        '''
        This function saves the Class data to file
        and prints what it is saving
        '''
        if hdf5 == True:
            pass
        else:
            with open(fileout,'a') as f:
                f.write("%g %g %g %g %g \n"%(self.isstream,self.vin,self.IMass,self.Mass,self.Slocal))

        
def get_ts(fname):  
    with h5.File(fname,'r') as f:
        ts = np.array(f['TimeSeries'])
    return ts


