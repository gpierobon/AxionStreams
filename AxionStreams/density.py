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
