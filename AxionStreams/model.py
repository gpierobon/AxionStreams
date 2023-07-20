import os,sys
import numpy as np

path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(path)

from AxionStreams import sample

def get_alpha2(i):
    # Model 0: NFW
    if i == 0:
        al2 = 0.13
    # Model 1: PL 
    elif i == 1:
        al2 = 0.27
    return al2

def get_beta(i):
    # Model 0: NFW
    if i == 0:
        be = 3.47
    # Model 1: PL 
    elif i == 1:
        be = 1.5
    return be

def get_kappa(i):
    # Model 0: NFW
    if i == 0:
        ka = 3.54
    # Model 1: PL 
    elif i == 1:
        ka = 1.73
    return ka

def get_bmax(i):
    # Model 0: NFW
    if i == 0:
        bmax = 0.1*1e-3
    # Model 1: PL
    elif i == 1:
        bmax = 0.075*1e-3
    return bmax

def get_radius(ic,mass):
    rho = sample.sample_mc_density(ic)
    rad = (3*mass/(4*np.pi*rho))**(1./3.)
    return rad,rho
