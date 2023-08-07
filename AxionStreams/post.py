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


def average_streamN(streamM,void=0.08):
    '''
    '''
    rhoDM = 0.4/37.96*1e-9 # Msun/mpc^3
    rhovoidNorm = void     # 8% of rhoDM is in background already
    rholocal = rhoDM*(1-rhovoidNorm)
    streamN = rholocal/np.mean(streamM)
    return streamN

    
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


def saturate_DM(size,max_iter=int(np.ceil(1e6)))
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
