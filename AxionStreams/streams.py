import h5py as h5
import numpy as np

class Stream():
    def __init__(self,stream_id,verb=False):
        self.ID = stream_id 
        if verb:
            print("Created stream with ID", stream_id)
    
    def get_orbit(self, fname):
        with h5.File(fname,'r') as f:
            self.R = np.array(f['Orbit_%03d/R'%self.ID])
            self.x = np.array(f['Orbit_%03d/x'%self.ID])
            self.y = np.array(f['Orbit_%03d/y'%self.ID])
            self.z = np.array(f['Orbit_%03d/z'%self.ID])
            self.vx = np.array(f['Orbit_%03d/vx'%self.ID])
            self.vy = np.array(f['Orbit_%03d/vy'%self.ID])
            self.vz = np.array(f['Orbit_%03d/vz'%self.ID])
        

def get_ts(fname):
    with h5.File(fname,'r') as f:
        ts = np.array(f['TimeSeries'])
    return ts


