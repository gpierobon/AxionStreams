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
ISO_MERG  = 0

fin     = path+'/orbit_data/Sun/orbits_%d_d%02d.hdf5'%(np.log10(N_SAMPLES),READ_FILE)
fout     = path+'/stream_data/merged/streams_%d_d%02d.txt'%(np.log10(N_SAMPLES),READ_FILE)
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



# Some analysis after build and perturb

#encounters = np.array([stream.N_encounters for stream in streams],dtype=object)
#enc_deg    = np.array([stream.enc_deg for stream in streams],dtype=object)
#ts_enc     = np.array([ts[stream.t_enc_ind] for stream in streams],dtype=object)

#blist      = [stream.blist for stream in streams]



# Mremain

masses = np.array([stream.Mass for stream in streams],dtype=object)
print("Min Mremain %g"%np.min(masses))

# Disruption time

#dis_times = np.array([ts[stream.t_disr_ind] for stream in streams], dtype=object)
#print(dis_times)

# Stream count 

stream_count = np.array([stream.isstream for stream in streams],dtype=object)
print("%d/%d Miniclusters are fully disrupted "%(np.sum(stream_count),N_SAMPLES))

Slocal = np.array([stream.Slocal for stream in streams if stream.Slocal < 1e-5],dtype=object)

#maxMLocal = np.max(np.array([stream.Mlocal for stream in streams if stream.Mlocal < 1e-5],dtype=object))
maxSLocal = np.max(np.array([stream.Slocal for stream in streams if stream.Slocal < 1e-5],dtype=object))
meanSLocal = np.mean(np.array([stream.Slocal for stream in streams if stream.Slocal < 1e-5],dtype=object))
#minMLocal = np.min(np.array([stream.Mlocal for stream in streams if stream.Mlocal < 1e-5],dtype=object))
#minSLocal = np.min(np.array([stream.Slocal for stream in streams if stream.Slocal < 1e-5],dtype=object))
#print("Max M_loc %g"%maxMLocal)
print("Max Stream_loc %g"%maxSLocal)
print("Mean Stream_loc %g"%meanSLocal)
print(np.sum(Slocal))
#print("Min M_loc %g"%minMLocal)
#print("Min Stream_loc %g"%minSLocal)

#print("[ENCOUNTERS] Changed encounter degeneracy %d/%d "%(np.argwhere(enc_deg>1).size,N_SAMPLES))
#print("[ENCOUNTERS] Max encounter degeneracy %d "%(np.max(enc_deg[np.where(enc_deg>1)])))

#Mloss = streams[0].MassLoss
#print(Mloss)
#print(encounters)
#print(enc_deg)
#print(blist[-1]*1e3)
#print(blist[-2]*1e3)


# Some plots
#pl.single_plot_3d(streams[0])
#angles = np.linspace(0,360,50)
#for i,j in zip(range(0,len(ts),20),angles):
#    fig = pl.plot_frame_with_stream(streams,i,lim=60,rot=j,sun=True)
#    fig.savefig('test_plots/frames/frame_%03d'%i,bbox_inches='tight')
