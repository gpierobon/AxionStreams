import os, sys, glob
import time
import h5py as h5
import numpy as np
from astropy import units as u
from scipy import stats 
from galpy import potential
from galpy.orbit import Orbit

path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(path)
from AxionStreams import orbit as orb

kpc = u.kpc
kms = u.km/u.s
deg = u.deg
Gyr = u.Gyr

N       = int(sys.argv[1])
T_Gyr   = 13.5
nframes = 2000
dtype   = 'Sun'

if dtype == 'Sun':
    f_pattern = path+"/orbit_data/Sun/orbits_%d_d*.hdf5"%np.log10(N)
else:
    f_pattern = path+"/orbit_data/orbits_%d_d*.hdf5"%np.log10(N)
files = glob.glob(f_pattern)

if files:
    latest = max(glob.glob(f_pattern))
    next_number = int(latest.split("_")[-1][1:3])+1
    fname = "orbits_%d_d%02d.hdf5"%(np.log10(N),next_number)
else:
    fname = "orbits_%d_d00.hdf5"%(np.log10(N))

if dtype == 'Sun':
    print("Creating file: ",path+'/orbit_data/Sun/'+fname)
else:
    print("Creating file: ",path+'/orbit_data/'+fname)


def dump_orbits(fname,N,T_Gyr,nframes):

    samples = orb.orbit_sampling(N,rmin=0.0,rmax=200)
    ts = np.linspace(0.0,T_Gyr*u.Gyr,nframes)

    with h5.File(path+'/orbit_data/'+fname,'w') as f:
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
            fg.attrs.create("Eccentricity", o.e())
            fg.create_dataset('R',data=o.R(ts))
            fg.create_dataset('x',data=o.x(ts))
            fg.create_dataset('y',data=o.y(ts))
            fg.create_dataset('z',data=o.z(ts))
            fg.create_dataset('vx',data=o.vx(ts))
            fg.create_dataset('vy',data=o.vy(ts))
            fg.create_dataset('vz',data=o.vz(ts))

def dump_Sun_orbits(fname,N,T_Gyr,nframes):

    ts = np.linspace(0.0,T_Gyr*u.Gyr,nframes)
    mean = np.array([0,0,0])
    vc = 220
    sig = vc/np.sqrt(2) # km/s
    cov = np.array([[sig**2, 0, 0], [0, sig**2,0],[0,0,sig**2]])
    vels = stats.multivariate_normal.rvs(mean,cov,size=N)
    vin = np.sqrt(vels[:,0]**2+vels[:,1]**2+vels[:,2]**2)
    ave = np.mean(vin)

    # Filter values larger than 530 km/s
    ind = np.where(vin > 530)
    print("Found %d unbound orbits"%(len(ind[0])))

    for j in ind[0]:
        vels[j,0] = np.random.uniform(-1,1)*ave
        vels[j,1] = np.random.uniform(-1,1)*ave
        vels[j,2] = np.random.uniform(-1,1)*ave

    vin = np.sqrt(vels[:,0]**2+vels[:,1]**2+vels[:,2]**2)
    print("Max v_in = %5f"%(np.max(vin)))

    with h5.File(path+'/orbit_data/Sun/'+fname,'w') as f:
        f.create_dataset('TimeSeries',data=ts)
        start = time.time()
        for i in range(N):
            display = np.arange(0,N,np.ceil(N*0.04))
            if i in display:
                print("Walltime %g, %d/%d ..."%(time.time()-start,i,N))
            R = 8
            z = 0
            phi = 0
            vR = vels[i,0]
            vT = vels[i,1]
            vz = vels[i,2]
            
            o = Orbit(vxvv=[R*kpc,vR*kms,vT*kms,z*kpc,vz*kms,phi*deg])
            o.integrate(ts,potential.MWPotential2014)

            fg = f.create_group('Orbit_%.3d'%i)
            #fg.attrs.create("Eccentricity", o.e())
            #fg.create_dataset('R',data=o.R(ts))
            fg.create_dataset('x',data=o.x(ts))
            fg.create_dataset('y',data=o.y(ts))
            fg.create_dataset('z',data=o.z(ts))
            fg.create_dataset('vx',data=o.vx(ts))
            fg.create_dataset('vy',data=o.vy(ts))
            fg.create_dataset('vz',data=o.vz(ts))
            

if __name__ == "__main__":
    if dtype == 'Sun':
        dump_Sun_orbits(fname,N,T_Gyr,nframes)
    else:
        dump_orbits(fname,N,T_Gyr,nframes)
