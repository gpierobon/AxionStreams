import numpy as np

def filter_streamM(arr):
    '''
    '''
    ma = np.logical_or(np.isnan(arr[:,4]),np.logical_or(np.isinf(arr[:,4]),arr[:,4]==0))
    streamM = arr[~ma,4]
    print("Available samples: %d (%.1f percent of raw data)"%(len(streamM),100*len(streamM)/len(arr[:,4])))
    return streamM

