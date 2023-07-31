import scipy
import numpy as np


def draw_random_b(bmin,bmax,nm=1000,n_samples=1):
    bvals = np.geomspace(bmin,bmax,nm)
    PM = 2*bvals/(bmax**2 - bmin**2)
    ic = np.random.choice(np.arange(0,nm),p=PM/np.sum(PM),size=n_samples)
    b_sample = bvals[ic]
    return b_sample


def draw_random_Mass(isomer,ic=1,n_samples=1,nm=1000):
    '''
    '''
    if isomer == 1: # Isolated
        # Jaxions
        if ic == 0:
            M_min = 1e-15
            M_max = 1e-10
            gamma = 0.7
        # Moore
        elif ic == 1:
            M_min = 1e-16
            M_max = 1e-12
            gamma = 0.8
        # Spax
        elif ic == 2:
            M_min = 1e-18
            M_max = 3e-14
            gamma = 0.68

    elif isomer == 0: # Merged
        if ic == 1: 
            M_min = 1e-12
            M_max = 5e-7
            gamma = 0.5

    Mvals = np.geomspace(M_min,M_max,nm)
    PM = Mvals**-gamma
    ic = np.random.choice(np.arange(0,nm),p=PM/np.sum(PM),size=n_samples)
    M_sample = Mvals[ic]
    return M_sample


def draw_random_Mass_old(isomer,ic=1):
    '''
    '''
    if isomer == 0: # Isolated
        # Jaxions
        if ic == 0:
            M_min = 1e-15
            M_max = 1e-10
            gamma = 0.7
        # Moore
        elif ic == 1:
            M_min = 1e-16
            M_max = 1e-12
            gamma = 0.8
        # Spax
        elif ic == 2:
            M_min = 1e-18
            M_max = 3e-14
            gamma = 0.68
    elif isomer == 1: # Merged
        if ic == 1: 
            M_min = 1e-12
            M_max = 1e-6
            gamma = 0.5
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






