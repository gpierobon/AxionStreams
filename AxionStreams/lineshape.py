import numpy as np

# Constants

v_pec = np.array([11.1,12.2,7.3])
vv_earthrev = 29.79
eccentricity = 0.016722
e1 = np.array([0.9941,0.1088,0.0042])
e2 = np.array([-0.0504,0.4946,-0.8677])
w_p = 2*np.pi/365 # orbital freq.
t1 = 79
ve = 29.79 # Earth's revolution
vrot = 0.47 # Earth's rotation
AstronomicalUnit = 1.49597892e11 # Astronomical Unit

def LabVelocitySimple(day,v_LSR=233.0):
    '''
    day measured from Jan1
    '''
    vsun = np.array([0.0,v_LSR,0.0])+v_pec
    v_lab = vsun + EarthVelocity(day)
    return v_lab

def EarthVelocity(day):
    '''
    Second order in eccentricity
    day measured from Jan1
    '''
    lambda_p = 102.93*np.pi/180.0
    th = w_p*(day-t1)
    v_E = np.cos(th)*(e1-2*eccentricity*np.sin(lambda_p)*e2) \
          +np.sin(th)*(e2+2*eccentricity*np.sin(lambda_p)*e1) \
          -eccentricity*(np.cos(2*th)*(np.cos(lambda_p)*e1-np.sin(lambda_p)*e2) \
          +np.sin(2*th)*(np.sin(lambda_p)*e1+np.cos(lambda_p)*e2))
    return vv_earthrev*v_E

def SpeedDist_Isotropic(v,day,v_LSR=233.0,sig=164.75,v_esc=528.0,\
                        v_shift=np.array([0.0,0.0,0.0])):
    '''
    ## Isotropic speed distribution
    v = speeds (in km/s)
    v_LSR = Local standard of rest
    sig = 1d dispersion
    v_esc = escape speed
    v_shift = any shift to v_lab needed (used for the stream velocity)
    '''
    # Analytic form for the speed distribution:
    v_lab = LabVelocitySimple(day,v_LSR=v_LSR)-v_shift
    v_e = np.sqrt(sum(v_lab**2.0))
    v0 = sig*np.sqrt(2.0)
    fv1 = (1.0/(np.sqrt(2*np.pi)))*(v/(v_e*sig))\
        *(np.exp(-(v**2.0+v_e**2.0-2.0*v*v_e)/(2*sig**2.0))\
        -np.exp(-(v**2.0+v_e**2.0+2.0*v*v_e)/(2*sig**2.0)))\
        *((v)<(v_esc+v_e))
    fv1 = fv1/np.trapz(fv1,v)
    return fv1


def get_A(streams):
    rhoDM = 0.01
    return streams[:,2]*1e-3/(np.pi*streams[:,-4]*(streams[:,6])**2)/rhoDM


def get_bkg(n=3000):
    c_km = 3e8/1000
    m_a = 3e-6
    m_a_s = m_a/6.58e-16
    vmin = 0.1
    vmax = 700

    omega_min = (m_a_s)*(1+(vmin/c_km)**2.0/2.0)
    omega_max = (m_a_s)*(1+(vmax/c_km)**2.0/2.0)
    axionBW = (omega_max-omega_min)
    omega = np.linspace(omega_min,omega_max+0.21*axionBW,n)
    omega_min = omega[0]
    omega_max = omega[-1]
    domega = omega[1]-omega[0]
    v = c_km*np.sqrt(2*(omega-m_a_s)/omega)
    v[omega<=m_a_s] = 0.0
    dv = (1.0/m_a_s)*(c_km/v)*c_km
    dv[omega<=m_a_s] = 0.0
    f0 = SpeedDist_Isotropic(v,day=0)
    line0 = f0/v
    return 1e6*(omega/m_a_s-1), 1e6*line0

def get_stream_lines(streams,stream_count,Nst,ndays=1000,n=3000,a_noise=1e-2):
    c_km = 3e8/1000
    m_a = 3e-6
    m_a_s = m_a/6.58e-16
    vmin = 0.1
    vmax = 700

    omega_min = (m_a_s)*(1+(vmin/c_km)**2.0/2.0)
    omega_max = (m_a_s)*(1+(vmax/c_km)**2.0/2.0)
    axionBW = (omega_max-omega_min)
    omega = np.linspace(omega_min,omega_max+0.21*axionBW,n)
    omega_min = omega[0]
    omega_max = omega[-1]
    domega = omega[1]-omega[0]
    v = c_km*np.sqrt(2*(omega-m_a_s)/omega)
    v[omega<=m_a_s] = 0.0
    dv = (1.0/m_a_s)*(c_km/v)*c_km
    dv[omega<=m_a_s] = 0.0
    f0 = SpeedDist_Isotropic(v,day=0)
    line0 = f0/v

    alpha_noise = np.random.normal(size=n)*a_noise

    v_str = [np.array([streams[i,7],streams[i,8],streams[i,9]]) for i in range(Nst)]
    sigma = [streams[i,1] for i in range(Nst)]
    f_str  = [np.zeros(shape=(ndays,n)) for i in range(Nst)]
    days = np.linspace(0,2*365-2*365/ndays,ndays)
    fstr = np.array(f_str)[:,0]

    Amp = get_A(streams)
    fs = np.nansum(np.array(fstr)*Amp[:,np.newaxis],axis=0)

    return 1e6*(omega/m_a_s-1),1e6*(1+alpha_noise)*(line0+fs/v)
