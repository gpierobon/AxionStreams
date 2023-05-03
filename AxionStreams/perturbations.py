import os,sys
import numpy as np

#path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
path = os.getcwd()
sys.path.append(path)

class Perturb:
    def __init__(self,dens_prof):
        self.profile = dens_prof

    def load_perturb_data(self):
        dE, dM, fej, fub = np.loadtxt(path+"/perturb_data/Perturbations_%s.txt"%self.profile, unpack=True)
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
    

