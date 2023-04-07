import numpy as np
from galpy import potential,df
from galpy.orbit import Orbit
from astropy import units as u

kpc = u.kpc
kms = u.km/u.s
deg = u.deg
Gyr = u.Gyr

def orbit_sampling(N_samples):
    pot  = potential.NFWPotential()
    distr_funct = df.isotropicNFWdf(pot=pot,rmax=200*kpc,vo=230*kms,ro=8.2*kpc)
    samples = distr_funct.sample(n=N_samples,return_orbit=True)
    return samples 

def orbit_coordinates(samples,T_Gyr,nframes,ctype='Cartesian'):
    N_samples = len(samples)
    ts = np.linspace(0.0,T_Gyr*u.Gyr,nframes)
    Coords = np.zeros(shape=(nframes,samples.dim(),N_samples))

    for i in range(N_samples):
        R   = samples[i].R()
        z   = samples[i].z()
        vR  = samples[i].vR()
        vT  = samples[i].vT()
        vz  = samples[i].vz()
        phi = samples[i].phi()
        o = Orbit(vxvv=[R*kpc,vR*kms,vT*kms,z*kpc,vz*kms,phi*deg])
        o.integrate(ts,potential.MWPotential2014)
        if ctype == 'Cartesian':
            Coords[:,:,i] = np.column_stack((o.x(ts),o.y(ts),o.z(ts)))
        elif ctype == 'Cylindrical':
            Coords[:,:,i] = np.column_stack((o.R(ts),o.phi(ts),o.z(ts)))
        else:
            print("Type not valid!")
    return Coords


