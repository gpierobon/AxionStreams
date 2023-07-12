import scipy
import numpy as np

# Credits to Bradley Kavanagh (https://github.com/bradkav/axion-miniclusters/blob/main/code/tools.py)
def inverse_transform_sampling(function, xmin, xmax, nbins=1000, n_samples=1000, logarithmic=False):
    if (logarithmic):
        bins = np.geomspace(xmin, xmax, num=nbins)
    else:
        bins = np.linspace(xmin, xmax, num=nbins)
    pdf = function(np.delete(bins,-1) + np.diff(bins)/2)
    Norm = np.sum(pdf*np.diff(bins))
    pdf /= Norm
    cumul_values = np.zeros(bins.shape)
    cumul_values[1:] = np.cumsum(pdf*np.diff(bins))
    inv_cdf = scipy.interpolate.interp1d(cumul_values, bins)
    r = np.random.rand(n_samples)
    return inv_cdf(r)

def draw_random_b(bmin,bmax,size):
    P_b = lambda b: 2*b/(bmax**2 - bmin**2)
    blist = inverse_transform_sampling(P_b, bmin, bmax, n_samples=size)
    return blist

def draw_random_Mass(ic):
    '''
     M**pow1*(M<M_break) + M**pow2*(M>M_break)*(M_break**pow1/(M_break**pow2))
    '''
    # Jaxions
    if ic == 0:
        M_min = 1e-15
        M_max = 1e-10
        gamma = 0.8
    # Moore
    elif ic == 1:
        M_min = 1e-16
        M_max = 1e-11
        gamma = 0.9
    # Spax
    elif ic == 2:
        M_min = 1e-18
        M_max = 3e-14
        gamma = 0.68

    norm = 1
    P_M = lambda M: norm*(M)*(-gamma)
    mass = inverse_transform_sampling(P_M,M_min,M_max,n_samples=1,logarithmic=True)
    return mass

def sample_mc_density(ic):
    '''
    Input from N-body simulations
    '''
    # Jaxions
    if  ic == 0:
        rhomean = 1.733e3*1e9
        rhostd = 9.957e4*1e9
    # Moore
    elif ic == 1:
        rhomean = 1.246e3*1e9 # Fix
        rhostd = 1.916e4*1e9  # Fix
    # Spax
    elif ic == 2:
        rhomean = 1.384e3*1e9 # Fix
        rhostd = 3.792e3*1e9  # Fix
    sam = np.abs(scipy.stats.norm.rvs(loc=rhomean, scale=rhostd, size=1))    
    return sam






