import numpy as np
from galpy import potential,df
from galpy.orbit import Orbit
from astropy import units as u

kpc = u.kpc
kms = u.km/u.s
deg = u.deg
Gyr = u.Gyr

def orbit_sampling(N_samples,rmin=0.0,rmax=200):
    pot  = potential.NFWPotential(a=19.0*kpc)
    distr_funct = df.isotropicNFWdf(pot=pot,ro=8.0*kpc,vo=220*kms,rmax=rmax*kpc)
    samples = distr_funct.sample(n=N_samples,return_orbit=True,rmin=rmin*kpc)
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
        phi = samples[i].phi()*180/np.pi
        o = Orbit(vxvv=[R*kpc,vR*kms,vT*kms,z*kpc,vz*kms,phi*deg])
        o.integrate(ts,potential.MWPotential2014)
        if ctype == 'Cartesian':
            Coords[:,:,i] = np.column_stack((o.x(ts),o.y(ts),o.z(ts)))
        elif ctype == 'Cylindrical':
            Coords[:,:,i] = np.column_stack((o.R(ts),o.phi(ts),o.z(ts)))
        else:
            print("Type not valid!")
    return Coords


def get_profile(nsamples,bins=500,rmin=0.0,rmax=200):
    '''

    Reproduces the density profile from the sampled miniclusters
    
    Use: r,rho = get_profile(nsamples)
    
    Input: 
        - nsamples: the sample array  

    Output: 
       - r: radial coordinates in kpc
       - rho: density profiles in a.u.
    '''
    N_samples = int(nsamples)
    pot  = potential.NFWPotential(a=19.0*kpc)
    distr_funct = df.isotropicNFWdf(pot=pot,ro=8.0*kpc,vo=220*kms,rmax=rmax*kpc)
    samples = distr_funct.sample(n=N_samples,return_orbit=True,rmin=rmin*kpc)
    R = samples.R()
    dP,rb = np.histogram(np.log10(R),bins=bins,range=[0,np.log10(200)]) 
    rb = 10**rb
    rc = (rb[0:-1]+rb[1:])/2
    dr = rb[1:]-rb[0:-1]
    rho = (1/(4*np.pi*rc**2))*(dP/dr)
    rho = rho/rho[np.argmin(np.abs(rc-19.0))]
    return rc, rho


def disk_encounters(orbits,maxrad=0):
    '''
    Computes the number of disk encounters (crossing z=0) for each orbit, discards encounters with radii larger 
    than maxrad
    Input: 
        - Orbits coordinates from orbit_coordinates function with shape (nframes,3,nsamples)
        - Maximum radius to consider the disk encounter in kpc
    Output:
        - List with number of encounters for each orbit, shape (nsamples)
    '''
    enc_li = []
    for i in range(orbits.shape[-1]):
        z_path = orbits[:,:,i][:,2]
        R_path  = np.sqrt(z_path**2+(orbits[:,:,i][:,1])**2+(orbits[:,:,i][:,0])**2)
        enc_array = np.abs(np.diff(np.sign(z_path)))
        if maxrad > 0:
            mask = np.argwhere(R_path <= maxrad)
            enc_array = enc_array[mask-1] # minus 1
        enc  = int(enc_array.sum()/2)
        enc_li.append(enc)
    return np.array(enc_li)

