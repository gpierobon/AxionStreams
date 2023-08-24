import os,sys
import numpy as np

path1 = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
path2 = os.getcwd()
sys.path.append(path1)
sys.path.append(path2)

class Perturb:
    def __init__(self,dens_prof):
        self.profile = dens_prof

    def load_perturb_data(self):
        # FIX
        try:
            dE, dM, fej, fub = np.loadtxt(path1+"/perturb_data/Perturbations_%s.txt"%self.profile, unpack=True)
        except:
            dE, dM, fej, fub = np.loadtxt(path2+"/perturb_data/Perturbations_%s.txt"%self.profile, unpack=True)
        self.dE = dE
        self.dM = dM
        self.fej = fej
        self.fub = fub

    def dM_interp(self,x):
        return np.interp(x, self.dE, self.dM, left=0.0, right=1.0)

    def fej_interp(self,x):
        return np.interp(x, self.dE, self.fej, left=0.0, right=1.0)

    def fub_interp(self,x):
        return np.interp(x, self.dE, self.fub, left=0.0, right=1.0)
    

