import os,sys
import numpy as np

path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(path)

from AxionStreams import sample

def get_perturb_parameters(i):
    # Model 0: NFW
    if i == 0:
        al2 = 0.13
        be = 3.47
        ka = 3.54
        bmax = 0.1*1e-3
    # Model 1: PL 
    elif i == 1:
        al2 = 0.27
        be = 1.5
        ka = 1.73
        bmax = 0.1*1e-3
    return al2,be,ka,bmax

def get_rho_radius_iso(ic,mass):
    # This dilution factor is calibrated by N-body simulations of isolated MCs
    dil = 200
    rho = sample.sample_mc_density(ic)/dil
    rad = (3*mass/(4*np.pi*rho))**(1./3.)
    return rad,rho

def get_radius_merged(mass,mtype='Xiao'):
    '''
    '''
    if mtype == 'Dai':
        rad = get_radius_Dai(mass)
    elif mtype == 'Xiao':
        rad = get_radius_Xiao(mass)
    elif mtype == 'Kavanagh':
        rad = get_radius_Kavanagh(mass)
    return rad

def get_radius_Xiao(mass,Amp=1.3e-3,M0=2.47e-10,z=19):
    '''
    Relation from Xiao et al. 2019
    '''
    rad = 1.4e4/(1+z)/np.sqrt(mass/(Amp*M0))*1e-3*3.7e-3/0.7*(mass/1e-6)**(5/6)*(Amp*M0/1e-11)**(-1/2)
    return rad
       

def get_radius_Dai(mass):
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
    return rad


def get_radius_Kavanagh(mass):
    '''
    Formulas from Kavanagh et al. 2020
    '''
    m_to_pc  = 3.240756e-17
    rad = 1e12*m_to_pc*(mass/1e-10)**(1/2)*1e-3
    return rad

def get_rho_mer(mass,rad):
    return 3*mass/(4*np.pi*rad**3)

def UniformStream(stream,ts,nvals=1000,dL=1e-6):
    '''
    Returns stream length, its corresponding time and mass loss over dL volume
    default: dL = 1e-6
    '''
    s2Gyr = 1/((365*24*3600)*1e9)
    km2kpc = 3.24078e-17
    tenc = ts[stream.t_enc_ind].flatten()
    if tenc[-1] == ts[-1]:
        tenc = tenc[:-1]
    nenc = len(tenc)

    sig = stream.Vdisp.flatten()[:nenc]
    Mloss = stream.MassLoss[:nenc]

    lstr = (sig*km2kpc/s2Gyr)*(ts[-1]-tenc)
    lstr = np.average(lstr, weights=Mloss)
    dMdl = Mloss/lstr

    lvals = np.linspace(-np.amax(lstr)*np.sqrt(2),np.amax(lstr)*np.sqrt(2),nvals)

    nenc = len(tenc.flatten())
    dMtot = np.zeros((nvals))
    for i in range(0,nenc):
        if dMdl[i]*lstr>0:
                dMtot += dMdl[i]
    dMtot = np.sum(Mloss)/np.trapz(dMtot,lvals)*dMtot
    vstr = np.sqrt(stream.vx[-1]**2 + stream.vy[-1]**2 + stream.vz[-1]**2)
    Tvals = lvals/(vstr*km2kpc/s2Gyr)

    return np.max(lvals),np.max(Tvals),dMtot[0]*dL


# Revisit
def GaussianStream(stream,ts,nvals=100,uniform=False):
    s2Gyr = 1/((365*24*3600)*1e9)
    km2kpc = 3.24078e-17
    tenc = ts[stream.t_enc_ind].flatten()
    if tenc[-1] == ts[-1]:
        tenc = tenc[:-1]
    nenc = len(tenc)
    
    sig = stream.Vdisp.flatten()[:nenc]
    Mloss = stream.MassLoss[:nenc]
    
    lstr = (sig*km2kpc/s2Gyr)*(ts[-1]-tenc)
    dMdl = Mloss/lstr

    lvals = np.linspace(-np.amax(lstr)*np.sqrt(2),np.amax(lstr)*np.sqrt(2),nvals)

    nenc = len(tenc.flatten())
    dMtot = np.zeros((nvals))
    for i in range(0,nenc):
        if uniform == True:
            if dMdl[i]*lstr[i]>0:
                dMtot += dMdl[i]
        else:
            if dMdl[i]*lstr[i]>0:
                dM1 = np.exp(-(lvals**2)/(2*(lstr[i]/2)**2))
                if sum(dM1)>0:
                    dM1 = dMdl[i]*dM1/np.trapz(dM1,lvals)
                    dMtot += dM1
    dMtot = sum(Mloss)/np.trapz(dMtot,lvals)*dMtot
    vstr = np.sqrt(stream.vx[-1]**2 + stream.vy[-1]**2 + stream.vz[-1]**2)
    Tvals = lvals/(vstr*km2kpc/s2Gyr)
    if uniform == True:
        return np.max(lvals),np.max(Tvals),dMtot[0]/1e6
    else:
        return lvals,Tvals,dMtot/1e6

