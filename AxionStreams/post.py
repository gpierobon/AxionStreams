import os,glob
import numpy as np

def filter_streamM(arr,verbose=True):
    '''
    '''
    ma = np.logical_or(np.isnan(arr[:,3]),np.logical_or(np.isinf(arr[:,3]),arr[:,3]==0))
    streamM = arr[~ma,3]
    streamm = arr[~ma,0]
    if verbose == True:
        print("Available samples: %d (%.1f percent of raw data)"%(len(streamM),100*len(streamM)/len(arr[:,4])))
    return streamM,streamm

def filter_array(arr,index,verbose=True):
    '''
    '''
    ma = np.logical_or(np.isnan(arr[:,index]),np.logical_or(np.isinf(arr[:,index]),arr[:,index]==0))
    filtered = arr[~ma]
    if verbose == True:
        print("Available samples: %d (%.1f percent of raw data)"%(len(filtered),100*len(filtered)/len(arr[:,index])))
    return filtered

def average_streamN(streamM,void=0.08):
    '''
    '''
    rhoDM = 0.4/37.96*1e-9 # Msun/mpc^3
    rhovoidNorm = void     # 8% of rhoDM is in background already
    rholocal = rhoDM*(1-rhovoidNorm)
    streamN = rholocal/np.mean(streamM)
    return streamN

def add_random_coords(file,localsize=1e-9):
    '''
    Localsize has kpc units, so 1 mpc = 1e-9
    '''
    x,y,z = [],[],[]
    for i in range(len(file[:,0])):
        u1,u2,u3 = np.random.uniform(size=3)
        cos_t = 2*u1 - 1
        sin_t = np.sqrt(1 - cos_t**2)
        phi = 2*np.pi*u2
        r = localsize*u3**(1/3)
        del_x = r*sin_t*np.cos(phi)
        del_y = r*sin_t*np.sin(phi)
        del_z = r*cos_t
        x.append(del_x - 8.2)
        y.append(del_y)
        z.append(del_z)
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    cfile = np.column_stack((file,x,y,z))
    return cfile

def final_samples(size,isomer,localsize=1e-9):
    '''
    '''
    path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
    if isomer == 0:
        f_pattern = path+"/stream_data/merged/streams_%.1d_d*.txt"%size
        ofile = path+'/stream_data/WN_%d.txt'%size
        off = '/stream_data/WN_%d.txt'%size
    else:
        f_pattern = path+"/stream_data/iso/streams_%.1d_d*.txt"%size
        ofile = path+'/stream_data/NL_%d.txt'%size
        off = '/stream_data/NL_%d.txt'%size

    files = glob.glob(f_pattern)
    fsamples = np.empty((0, 14))
    for file in files:
        inputf = np.loadtxt(file) 
        nfile = add_random_coords(inputf,localsize=localsize)
        nfile  = filter_array(nfile,3)
        fsamples = np.vstack((fsamples, nfile))
    header = "# MC/Stream,v_in,vel.disp,Mloc_stream,M_i,M_f,R_i,R_f,vx,vy,vz,x,y,z"
    np.savetxt(ofile, fsamples, header=header,delimiter=',')
    print("Output %s has %d samples "%(off,len(fsamples[:,0])))
    

def draw_Slocal(size,prob=0.7,verbose=False):
    if np.random.random() > prob:
        if verbose == True:
            print('Merged')
        num = np.random.ranint(10)
        fi = np.loadtxt(path+'/stream_data/merged/streams_%d_d%.2d.txt'%(size,num))
    else:
        if verbose == True:
            print('Isolated')
        num = np.random.ranint(10)
        fi = np.loadtxt(path+'/stream_data/iso/streams_%d_d%.2d.txt'%(size,num))

    streamM,_ = post.filter_streamM(fi,verbose=False)
    mass_sample = np.random.choice(streamM)
    return mass_sample


def saturate_DM(size,max_iter=int(np.ceil(1e6))):
    totmass = 0
    totsize = int(10**size)
    malist = []
    for i in range(max_iter):
        ma = draw_Slocal(size)
        malist.append(ma)
        totmass += ma
        if i > totsize-1:
            print("Reached sample size!")
            break
        if totmass > 1e-11:
            print("Saturated the DM total mass with %d streams"%i)
            break
    return np.array(malist),totmass
