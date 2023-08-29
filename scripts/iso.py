import os,sys
import time
import scipy 
import h5py as h5
import numpy as np
from galpy import potential,df
from galpy.orbit import Orbit
from astropy import units as u
from tqdm import tqdm

path = os.getcwd()
sys.path.append(path)

# Scripts
from AxionStreams import plot as pl
from AxionStreams import orbit as orb
from AxionStreams import density as dens
import AxionStreams.streams as st

# Ignore RunTimeWarnings, data is filtered after
import warnings
warnings.filterwarnings("ignore")

# Set units
kpc = u.kpc
kms = u.km/u.s
deg = u.deg
Gyr = u.Gyr

INTERACT  = int(sys.argv[1])
N_SAMPLES = int(sys.argv[2])
READ_FILE = int(sys.argv[3])
ISO_MERG  = 1
HDF5      = True


fin     = path+'/orbit_data/Sun/orbits_%d_d%02d.hdf5'%(np.log10(N_SAMPLES),READ_FILE)
fout     = path+'/stream_data/iso/streams_%d_d%02d'%(np.log10(N_SAMPLES),READ_FILE)
ts = st.get_ts(fin)

streams = [st.Stream(i) for i in range(N_SAMPLES)]

if HDF5 == True:
    fileh = fout+'.hdf5'
    if os.path.exists(fileh):
        os.remove(fileh)

fo = open(fout+'.txt','w')

start = time.time()
if INTERACT == 0:
    for i in tqdm(range(N_SAMPLES)):
        streams[i].run(fin,fout,ISO_MERG,ts,hdf5=HDF5)
else:
    for i in range(N_SAMPLES):
        display = np.arange(0,N_SAMPLES,np.ceil(N_SAMPLES*0.04))
        if i in display:
            print("Walltime: %g, %d/%d ..."%(time.time()-start,i,N_SAMPLES))
        streams[i].run(fin,fout,ISO_MERG,ts,hdf5=HDF5)

fo.close()

if INTERACT == 0:
    stream_count = np.array([stream.isstream for stream in streams],dtype=object)
    print("%d/%d Miniclusters are fully disrupted "%(np.sum(stream_count),N_SAMPLES))
