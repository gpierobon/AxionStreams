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



