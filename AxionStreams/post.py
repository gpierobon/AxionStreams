import numpy as np

def filter_streamM(arr):
    '''
    '''
    ma = np.logical_or(np.isnan(arr[:,4]),np.logical_or(np.isinf(arr[:,4]),arr[:,4]==0))
    streamM = arr[~ma,4]
    print("Available samples: %d (%.1f percent of raw data)"%(len(streamM),100*len(streamM)/len(arr[:,4])))
    return streamM


def average_streamN(streamM,void=0.08):
    '''
    '''
    rhoDM = 0.4/37.96*1e-9 # Msun/mpc^3
    rhovoidNorm = void     # 8% of rhoDM is in background already
    rholocal = rhoDM*(1-rhovoidNorm)
    streamN = rholocal/np.mean(streamM)
    return streamN
    

