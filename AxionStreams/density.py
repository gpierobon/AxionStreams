import numpy as np

def rhoR(x,y,z):
    '''
    Returns bulge density as a function of cartesian coords.
    Units are kpc for input and SolarMasses/pc^3 for output
    '''
    R = np.sqrt(x**2 + y**2)
    rp = np.sqrt(R**2.0+(z/0.5)**2.0)
    return 95.6*1000/((1.0+(rp/0.075))**1.8)*np.exp(-(rp/2.1)**2.0)

def rhoz1(x,y,z):
    '''
    Returns thin disk density as a function of cartesian coords.
    Units are kpc for input and SolarMasses/pc^3 for output
    '''
    R = np.sqrt(x**2 + y**2+z**2.0)
    return 816.6/(2*0.3)*np.exp(-np.abs(z)/0.3 - R/2.6)

def rhoz2(x,y,z):
    '''
    Returns thick disk density as a function of cartesian coords.
    Units are kpc for input and SolarMasses/pc^3 for output
    '''
    R = np.sqrt(x**2 + y**2)
    return 209.5/(2*0.9)*np.exp(-np.abs(z)/0.9 - R/3.6)

def stellar_density(x,y,z):
    '''
    Returns the full stellar density as a function of cartesian coords.
    Units are kpc for input and SolarMasses/pc^3 for output
    '''
    return rhoR(x,y,z)+rhoz1(x,y,z)+rhoz2(x,y,z)

def get_ts_encounters(ts,rho,samplesize):
    deg = 1
    new_samplesize = samplesize
    if samplesize > 0.3*len(ts):
        new_samplesize = int(0.3*len(ts))
        deg = int(round(samplesize/new_samplesize)) 
    sampled_ts = np.random.choice(ts, p=rho/np.sum(rho), size=new_samplesize, replace=False)
    sampled_ts = np.sort(sampled_ts)
    ts_indices = np.where(np.isin(ts, sampled_ts))[0]
    return new_samplesize,ts_indices,deg



def sample_from_stellar_density(ts,Orbits,j,samplesize):
    '''
    Sample density values from the calculated stellar density
    Returns sampled timestamps and encounter degeneracy 

    FIX: for now it onl takes a single orbit,
         there is an issue with time units (has to be dimensionless)
    '''
    deg = 1
    new_samplesize = samplesize
    
    if samplesize > 0.2*len(ts):
        new_samplesize = 0.2*len(ts)
        deg = samplesize/new_samplesize 

    rho = stellar_density(Orbits[:,0,j],Orbits[:,1,j],Orbits[:,2,j])
    sorted_indices = np.argsort(ts)
    sampled_ts = np.random.choice(ts, p=rho/np.sum(rho), size=new_samplesize, replace=False)
    sampled_indices = np.searchsorted(ts[sorted_indices], sampled_ts)
    sorted_sampled_indices = sorted_indices[sampled_indices]
    sampled_rho = rho[sorted_sampled_indices]
    return sampled_ts,deg
