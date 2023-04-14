import h5py as h5
import numpy as np
from astropy import units as u
from galpy import potential
from galpy.orbit import Orbit
from AxionStreams import orbit as orb

kpc = u.kpc
kms = u.km/u.s
deg = u.deg
Gyr = u.Gyr

N       = 1000
T_Gyr   = 13.5
nframes = 1000
fname   = 'orbits_%d_01.hdf5'%np.log10(N)

def dump_orbits(fname,N,T_Gyr,nframes):

    samples = orb.orbit_sampling(N,rmin=0.0,rmax=200)
    ts = np.linspace(0.0,T_Gyr,nframes)

    with h5.File('orbit_data/'+fname,'w') as f:
        f.create_dataset('TimeSeries',data=ts)

        for i in range(N):

            R   = samples[i].R()
            z   = samples[i].z()
            vR  = samples[i].vR()
            vT  = samples[i].vT()
            vz  = samples[i].vz()
            phi = samples[i].phi()*180/np.pi
            o = Orbit(vxvv=[R*kpc,vR*kms,vT*kms,z*kpc,vz*kms,phi*deg])
            o.integrate(ts,potential.MWPotential2014)

            fg = f.create_group('Orbit_%.3d'%i)
#            fg.attrs['ID'] = i 
            fg.create_dataset('x',data=o.x(ts))
            fg.create_dataset('y',data=o.y(ts))
            fg.create_dataset('z',data=o.z(ts))
            fg.create_dataset('vx',data=o.vx(ts))
            fg.create_dataset('vy',data=o.vy(ts))
            fg.create_dataset('vz',data=o.vz(ts))

dump_orbits(fname,N,T_Gyr,nframes)
