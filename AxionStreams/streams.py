import os, sys
import h5py as h5
import numpy as np
import scipy as scp

path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(path)

from AxionStreams import density as dens

class Stream():
    def __init__(self,stream_id,verb=False):
        self.ID = stream_id 
        if verb:
            print("Created stream with ID", stream_id)

    def print0(self,message):
        if self.ID == 0:
            print(message)

    def build(self, fname, ts):
    
        self.print0("[BUILD] Reading orbits from %s ... "%fname.split('/')[-1])
        self.read_orbit(fname)

        self.print0("[BUILD] Setting physical parameters ... ")
        self.set_parameters() 
        
        self.print0("[BUILD] Finding encounters ... ")
        self.set_encounters(ts)
    
    def perturb(self):
        # Here we pass the number of encounters and calculate the E_List and disr time
        # At the end recompute the energy 
        #self.get_stream_energies(self)
        pass

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
        self.bmax = 1e-5    # kpc 
        self.bmin = 1e-7    # kpc 
        self.alpha2 = 0.27 
        self.beta = 1.5
        
        self.rho = self.Mass/(4*np.pi/3.*self.Rad**3) # SolarMass/kpc**3 

    def set_encounters(self,ts):
        kms_to_kpcGyr = 1.02269 
        star_density = dens.stellar_density(self.x,self.y,self.z) 
        vel = np.abs(np.sqrt(self.vx**2+self.vy**2+self.vz**2))*kms_to_kpcGyr
        inte = 1e9*star_density*vel*np.pi*self.bmax**2 # Units are kpc,Gyr
        finterp = scp.interpolate.InterpolatedUnivariateSpline(ts, inte, k=1)
        xx = np.linspace(ts[0], ts[-1], 5*len(ts))
        qq = [finterp.integral(0, t) for t in xx]        
        
        self.N_encounters  = int(round(qq[-1]))
        new_enc,ts_enc,enc_deg = dens.get_ts_encounters(ts,star_density,self.N_encounters)
        self.N_encounters = new_enc
        self.ts_enc = ts_enc
        self.enc_deg = enc_deg
        self.print0("[ENCOUNTERS] Stream %d will undergo %d encounters"%(self.ID,self.N_encounters))                    
        
        del xx, qq, star_density

    def get_stream_energies(self):
        #self.Ebind = self.beta*G_N*self.M**2/self.R
        #self.Etot = 0.5*self.Mass*self.vdisp**2-self.Ebind
        pass 


def get_ts(fname):
    with h5.File(fname,'r') as f:
        ts = np.array(f['TimeSeries'])
    return ts


