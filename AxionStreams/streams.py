import os, sys
import h5py as h5
import numpy as np
import scipy as scp

path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(path)

from galpy import potential
from AxionStreams import density 
from AxionStreams import perturbations 
from AxionStreams import sample

class Stream():
    def __init__(self,stream_id,verb=False):
        self.ID = stream_id 
        if verb:
            print("Created stream with ID", stream_id)

    def print0(self,message):
        if self.ID == 0:
            print(message)

    def run(self, fname, ts):
    
        self.print0("Reading orbits from %s ... "%fname.split('/')[-1])
   
        self.read_orbit(fname)

        self.set_parameters() 
        
        self.set_encounters(ts)
        
        self.set_perturb_data()
        
        self.set_perturb_array()

        self.set_stream_energies()

        #print("Stream %d before: %d "%(self.ID,self.isstream))   
        self.perturb(debug=True)
        #print("Stream %d after: %d "%(self.ID,self.isstream))  
    # --------------------------------------------------

    def read_orbit(self, fname):
        with h5.File(fname,'r') as f:
            #self.e = f['Orbit_%03d/Eccentricity'%self.ID]  
            self.R = np.array(f['Orbit_%03d/R'%self.ID])
            self.x = np.array(f['Orbit_%03d/x'%self.ID])
            self.y = np.array(f['Orbit_%03d/y'%self.ID])
            self.z = np.array(f['Orbit_%03d/z'%self.ID])
            self.vx = np.array(f['Orbit_%03d/vx'%self.ID])
            self.vy = np.array(f['Orbit_%03d/vy'%self.ID])
            self.vz = np.array(f['Orbit_%03d/vz'%self.ID])
    
    def set_parameters(self):
        self.isstream = 0    # 0 for MC, 1 for stream
        
        # These will be drawn from simulations
        self.Mass = 1e-10    # Msun
        self.Rad = 3e-9      # kpc
        self.vdisp = 3e-5    # km/s    
        self.bmax = 1e-5     # kpc 
        self.bmin = 1e-6     # kpc 
        
        self.alpha2 = 0.27 
        self.beta = 1.5
        self.kappa = 1.73  
        self.prof = 'PL'     # 'PL', 'NFW' 

        # Transition radius, Eq. (7) of 2207.11276, f_b=6
        #self.bs = 6*np.sqrt(2*np.sqrt(self.alpha2)/(3*self.beta))*self.Rad 
        
        self.rho = self.Mass/(4*np.pi/3.*self.Rad**3) # SolarMass/kpc**3 
        self.G_N = 1 # Fix the units

    def set_encounters(self,ts):
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
        self.print0("[ENCOUNTERS] Stream %d will undergo %d encounters"%(self.ID,self.N_encounters))                    

        del xx, qq, star_density

    def set_stream_energies(self):
        self.Ebind = self.beta*self.G_N*self.Mass**2/self.Rad
        self.Etot = 0.5*self.Mass*self.vdisp**2-self.Ebind
        #self.Etot = (0.5*self.kappa/self.beta - 1)*self.Ebind

    def set_perturb_data(self):
        pdata = perturbations.Perturb(self.prof)
        pdata.load_perturb_data()
        self.dM_I = pdata.dM_interp
        self.fej_I = pdata.fej_interp
        self.fub_I = pdata.fub_interp

    def set_perturb_array(self):
        self.blist = sample.draw_random_b(self.bmin,self.bmax,self.N_encounters) # kpc
        self.vlist = 500*np.ones(self.N_encounters) # km/s
        self.DelE = 4./3.*self.G_N**2/(self.blist**4*self.vlist**2)*self.Mass*self.alpha2*self.Rad**2
        
        ''' Here we do vlist properly
        R_list = [np.sqrt(self.x[ti]**2+self.y[ti]**2) for ti in self.t_enc_ind]
        phi_list = [np.sqrt(self.x[ti]**2+self.y[ti]**2) for ti in self.t_enc_ind]
        vcirc_list = [potential.vcirc(potential.MWPotential2014, R=Ri,vo=220*kms, ro=8*kpc) for Ri in R_list]
        vx_list = [vcirc_list[j] for j in len(vcirc_list)]
        vy_list = [vcirc_list[j] for j in len(vcirc_list)]
        self.vlist = []
        '''

    def perturb(self,debug=False):
        counter = 0
        dE_max = 0
        Mloss_list = [] 
        for i in range(self.N_encounters):
            dE = self.DelE[i]/self.Ebind*self.enc_deg # Encounter degeneracy multiplies the encounte
            if dE > dE_max:
                dE_max = dE
            
            if (dE < 1e-4):
                dM = 0.0
                dE_remain = dE
            else:
                dM = self.dM_I(dE)*self.Mass
                dE_remain = dE*(1 - self.fej_I(dE)) - self.fub_I(dE)*self.Etot
            
            newMass = self.Mass - dM
            newRad = self.Rad

            E_f = self.Etot + dE_remain*self.Mass
            if debug == True:
                self.print0("dM %g dE_remain %g Etot %g E_f %g"%(dM,dE_remain,self.Etot,E_f))
            
            if self.isstream == 1:
                continue 
            
            if (E_f >= 0):
                self.isstream = 1
                self.t_disr_ind = self.t_enc_ind[counter]
            else:
                newRad = (0.5*self.kappa - self.beta)*self.G_N*newMass**2/E_f

            self.Rad = newRad
            self.Mass = newMass
            Mloss_list.append(dM)
            counter += 1
        
        self.set_stream_energies()
        self.MassLoss = np.array(Mloss_list)
        self.print0("For stream %d max dE is %.8f"%(self.ID,dE_max))

def get_ts(fname):
    with h5.File(fname,'r') as f:
        ts = np.array(f['TimeSeries'])
    return ts


