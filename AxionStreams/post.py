import os,glob
import numpy as np

def filter_streamM(arr,verbose=True):
    '''
    '''
    ma = np.logical_or(np.isnan(arr[:,3]),np.logical_or(np.isinf(arr[:,3]),arr[:,3]==0))
    streamM = arr[~ma,3]
    streamm = arr[~ma,0]
    if verbose == True:
        print("Available samples: %d (%.1f percent of raw data)"%(len(streamM),100*len(streamM)/len(arr[:,4])))
    return streamM,streamm

def filter_array(arr,index,verbose=True):
    '''
    '''
    ma = np.logical_or(np.isnan(arr[:,index]),np.logical_or(np.isinf(arr[:,index]),arr[:,index]==0))
    filtered = arr[~ma]
    if verbose == True:
        print("Available samples: %d (%.1f percent of raw data)"%(len(filtered),100*len(filtered)/len(arr[:,index])))
    return filtered

def average_streamN(streamM,void=0.08):
    '''
    '''
    rhoDM = 0.4/37.96*1e-9 # Msun/mpc^3
    rhovoidNorm = void     # 8% of rhoDM is in background already
    rholocal = rhoDM*(1-rhovoidNorm)
    streamN = rholocal/np.mean(streamM)
    return streamN

def add_random_coords(file,localsize=1e-6):
    '''
    Localsize has kpc units, so 1 mpc = 1e-6
    '''
    x,y,z = [],[],[]
    for i in range(len(file[:,0])):
        u1,u2,u3 = np.random.uniform(size=3)
        cos_t = 2*u1 - 1
        sin_t = np.sqrt(1 - cos_t**2)
        phi = 2*np.pi*u2
        r = localsize*u3**(1/3)
        del_x = r*sin_t*np.cos(phi)
        del_y = r*sin_t*np.sin(phi)
        del_z = r*cos_t
        x.append(del_x - 8.2)
        y.append(del_y)
        z.append(del_z)
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    cfile = np.column_stack((file,x,y,z))
    return cfile

def final_samples(size,isomer,fout='',localsize=1e-6):
    '''
    '''
    path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
    if isomer == 0:
        f_pattern = path+"/stream_data/merged/streams_%.1d_d*.txt"%size
        ofile = path+'/stream_data/WN%s.txt'%fout
        off = '/stream_data/WN%s.txt'%fout
    else:
        f_pattern = path+"/stream_data/iso/streams_%.1d_d*.txt"%size
        ofile = path+'/stream_data/NL%s.txt'%fout
        off = '/stream_data/NL%s.txt'%fout

    files = glob.glob(f_pattern)
    fsamples = np.empty((0, 15))
    for file in files:
        inputf = np.loadtxt(file) 
        nfile = add_random_coords(inputf,localsize=localsize)
        nfile  = filter_array(nfile,2)
        fsamples = np.vstack((fsamples, nfile))
    header = "# MC/Stream,vel.disp,Mloc_stream,M_i,M_f,R_i,R_f,vx,vy,vz,lmax,lstr,x,y,z"
    np.savetxt(ofile, fsamples, header=header,delimiter=',',fmt='%d %.2e %.2e %.2e %.2e %.2e %.2e %d %d %d %.2e %.2e %.2e %.2e %.2e')
    print("Output %s has %d samples "%(off,len(fsamples[:,0])))
    
def draw_sample(wn,nl,prob=0.7,verbose=False):
    path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
    if np.random.random() > prob:
        if verbose == True:
            print('Merged')
        ind = np.random.randint(0, wn.shape[0])
        sample = wn[ind]
        status = 0
    else:
        if verbose == True:
            print('Isolated') 
        ind = np.random.randint(0, nl.shape[0])
        sample = nl[ind]
        status = 1
    return sample,status

def get_initial_DM_mass(wn,nl,void=0.075,verbose=False):
    merg_count = 0
    totmass = 0
    for i in range(1000):
        sample,status = draw_sample(wn,nl,verbose=verbose)
        if status == 0:
            merg_count += 1
        mass = sample[3]
        totmass += mass
        if totmass > 4.1e-11*(1-void):
            if verbose == True:
                print("Saturated the DM total mass with %d minclusters (%d are merged)"%(i,merg_count))
            break
    return i,merg_count

def get_DM_mass(wn,nl,void=0.075,inputsize=2000000,verbose=False):
    '''
    sample has: 
        MC/Stream,vel.disp,Mloc_stream,M_i,M_f,R_i,R_f,vx,vy,vz,lmax,lstream,x,y,z
        
    Returns the samples
    '''
    totmass = 0
    merg_count = 0
    stream_count = 0
    mc_count = 0
    merg_zero = 0
    size = int(0.7*nl.shape[0])
    slist = []
    loop_size = np.minimum(size,inputsize)

    if verbose == True:
        print("Loop size: %d"%loop_size)
    for i in range(loop_size):
        sample,status = draw_sample(wn,nl)
        lmax = sample[10]
        lstr = sample[11]
        if status == 0:
            merg_count += 1
        slist.append(sample)

        rand_pos = np.random.uniform(0,lmax) 
        if rand_pos > lstr:
                mass = 0
                if status == 0:
                    merg_zero +=1 
        else:
            if np.abs(rand_pos) < 1e-6:
                if sample[4] > 0.0:
                    if verbose == True:
                        print('Found of a minicluster of mass %.2e'%sample[4])
                    mc_count += 1
                mass = sample[2]+sample[4]
            else:
                mass = sample[2]
                stream_count += 1 
        totmass += mass
        if totmass > 4.1e-11*(1-void):
            if verbose == True:
                print("Saturated the DM total mass with %d streams, of which %d are minclusters (total entiers : %d)"%(stream_count,mc_count,i))
            break
    return np.array(slist),merg_count,merg_zero,mc_count,stream_count,i,totmass

def get_DM_mass_2(wn,nl,void=0.075,inputsize=1000000,verbose=False):
    '''
    '''
    totmass = 0
    iso_count = 0
    size = int(0.7*nl.shape[0])
    slist = []
    cflag = 0
    loop_size = np.minimum(size,inputsize)
    if verbose == True:
        print("Loop size: %d"%loop_size)
    for i in range(loop_size):
        sample,status = draw_sample(wn,nl)
        if status == 1:
            iso_count += 1
        slist.append(sample)
        # Column number 3 is Mloc_stream, column 5 is remaining MC mass
        if cflag == 0:
            mass = sample[3] + sample[5]
        else:
            mass = sample[3]
        totmass += mass
                    
        if totmass > 0.7*4.1e-11*(1-void):
            cflag = 1
        
        if totmass > 4.1e-11*(1-void):
            if verbose == True:
                print("Saturated the DM total mass with %d streams"%i)
            break
    return np.array(slist),iso_count














