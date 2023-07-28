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
        bmax = 0.05*1e-3
    return bmax

def get_rho_radius_iso(ic,mass):
    rho = sample.sample_mc_density(ic)
    rad = (3*mass/(4*np.pi*rho))**(1./3.)
    return rad,rho

def get_radius_merged(mass,mtype='Xiao'):
    '''
    '''
    if mtype == 'Dai':
        rad = get_radius_mer_Dai(mass)
    elif mtype == 'Xiao':
        rad = get_radius_mer_Xiao(mass)
    return rad

def get_radius_mer_Xiao(mass,Amp=5e-4,M0=2.47e-10):
    '''
    Relation from Xiao et al. 2019
    '''
    rad = 1e-3*3.7e-3/0.7*(mass/1e-6)**(5/6)*(Amp*M0/1e-11)**(-1/2)
    #print(rad)
    return rad
       

def get_radius_mer_Dai(mass):
    '''
    Formulas from Dai, Miralda-Escude, 2020
    '''
    au_to_pc = 4.848102e-6
    if mass > 9.9e-10:
        c = 4
        #c = np.random.uniform(20,100)
    elif mass > 9.9e-11:
        c = np.random.uniform(5,100)
    else:
        c = np.random.uniform(100,500)
    
    rad = c*1e-3*973/0.7*au_to_pc*(mass/1e-6)**(5/6) 
    print(mass,rad,c)
    return rad

'''
def get_rho_mer(mass):
    rhos = 0.24*0.7**2*(mass/1e-6)**(-3/2) # Msun/pc^3
    rhos = rhos*1e9 # Msun/kpc^3
    rho = rhos/4**(1/3)
    return rho  
'''

def get_radius_mer_K(mass):
    '''
    '''
    m_to_pc  = 3.240756e-17
    rad = 1e12*m_to_pc*(mass/1e-10)**(1/2)*1e-3
    print(mass,rad)
    return rad

def get_rho_mer(mass,rad):
    return 3*mass/(4*np.pi*rad**3)



