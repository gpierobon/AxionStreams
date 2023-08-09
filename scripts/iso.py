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

fin     = path+'/orbit_data/Sun/orbits_%d_d%02d.hdf5'%(np.log10(N_SAMPLES),READ_FILE)
fout     = path+'/stream_data/iso/streams_%d_d%02d.txt'%(np.log10(N_SAMPLES),READ_FILE)
ts = st.get_ts(fin)

streams = [st.Stream(i) for i in range(N_SAMPLES)]

fo = open(fout,'w')

start = time.time()
if INTERACT == 0:
    for i in tqdm(range(N_SAMPLES)):
        streams[i].run(fin,fout,ISO_MERG,ts)
else:
    for i in range(N_SAMPLES):
        display = np.arange(0,N_SAMPLES,np.ceil(N_SAMPLES*0.04))
        if i in display:
            print("Walltime: %g, %d/%d ..."%(time.time()-start,i,N_SAMPLES))
        streams[i].run(fin,fout,ts)

fo.close()


maxN = np.max([stream.N_encounters for stream in streams])
encount = np.array([stream.N_encounters for stream in streams])
maxenc = np.max([stream.enc_deg for stream in streams])

stream_count = np.array([stream.isstream for stream in streams],dtype=object)

Slocal = np.array([stream.Slocal for stream in streams if stream.Slocal < 1e-5],dtype=object)
maxSLocal = np.max(np.array([stream.Slocal for stream in streams if stream.Slocal < 1e-5],dtype=object))
meanSLocal = np.mean(np.array([stream.Slocal for stream in streams if stream.Slocal < 1e-5],dtype=object))

#print('Max number of encounters: %d'%maxN)
#print('Max encounter deg: %d'%maxenc)
print("%d/%d Miniclusters are fully disrupted "%(np.sum(stream_count),N_SAMPLES))
print("Max Stream_loc %g"%maxSLocal)
print("Mean Stream_loc %g"%meanSLocal)
print("Sum of Stream_loc %g"%(np.sum(Slocal)))