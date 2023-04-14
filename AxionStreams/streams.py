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
    
    def build(self, fname, profile, verb=True):
        if verb == True:
            print("[BUILD] Reading orbits from disk ... "%fname)
        self.read_orbit(fname)

        if verb == True:
            print("[BUILD] Setting physical parameters ... "%fname)
        
        self.set_parameters(profile) # Set structural parameters
        self.set_size() # Set mass, radius
        
        if verb == True:
            print("[BUILD] Sampling encounters ... "%fname)
        self.set_encounters()

    def read_orbit(self, fname):
        with h5.File(fname,'r') as f:
            self.R = np.array(f['Orbit_%03d/R'%self.ID])
            self.x = np.array(f['Orbit_%03d/x'%self.ID])
            self.y = np.array(f['Orbit_%03d/y'%self.ID])
            self.z = np.array(f['Orbit_%03d/z'%self.ID])
            self.vx = np.array(f['Orbit_%03d/vx'%self.ID])
            self.vy = np.array(f['Orbit_%03d/vy'%self.ID])
            self.vz = np.array(f['Orbit_%03d/vz'%self.ID])
    
    def set_parameters(self,profile):
        if profile == "PL":
            self.alpha = 0.5196 
            self.beta = 1.5 
            self.b_max = 0.1e-3 # ~0.1 pc, might fix
        elif profile == "NFW":
            self.alpha = 0.3606 
            self.beta = 3.47
            self.b_max = 0.01e-3 # ~0.01 pc might fix
        elif profile == "SIM":
            pass
            # To implement sampling from simulation data
        else:
            raise ValueError("Profile not known, select bewteen PL,NFW and SIM")

    def set_size(self):
        #self.Mass = 0.0
        #self.Rad = 0.0
        #self.conc = 0.0
        pass

    def set_encounters(self,ts):
        star_density = dens.stellar_density(self.x,self.y,self.z)
        inte = star_density*np.pi # Times the velocity of the mc/stream and bmax^2
        finterp = scp.interpolate.InterpolatedUnivariateSpline(ts, inte, k=1)
        xx = np.linspace(ts[0], ts[-1], 5*len(ts))
        qq = [finterp.integral(0, t) for t in xx]
        self.N_encounters  = int(round(qq[-1]))
        # Add them to a time series 
        # self.t_encounters = 
        del xx, qq, star_density


    def get_stellar_density(self):
        self.star_dens = dens.stellar_density(self.x,self.y,self.z)


def get_ts(fname):
    with h5.File(fname,'r') as f:
        ts = np.array(f['TimeSeries'])
    return ts


