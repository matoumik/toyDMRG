from numpy import array, zeros
import numpy as np

class DMRG:
    def __init__(self, DMRG_info):
        self.info = DMRG_info
        self.iter = 0
        self.sweep = 0
        self.nleft = 0
        self.nright = self.nsites
        
    def H_contr_L(self, ncontr):
        for i in range(ncontr):
            if i == 0:
                pass

    
    
class DMRG_info:
    def __init__(self):
        self.M_start = 20
        self.nsites = 16


class MPS:
    def __init__(self, DMRG_info):
        self.nsites = DMRG_info.nsites
        self.sitedim = 2
        self.nbonds = self.nsites - 1
        self.init_M(DMRG_info.M_start)
        self.DBSS= False #TODO      
        

    def init_M(self, Mstart):
        M_temp = list()
        k = self.nbonds
        for i in range(self.nbonds):
            if (self.sitedim**(i+1)<Mstart):
                M_temp.append(self.sitedim**(i+1))
            elif self.sitedim**(k-i)<Mstart:
                M_temp.append(self.sitedim**(k-i))
            else:
                M_temp.append(Mstart)
        self.M = array(M_temp)
        print(self.M)
        
    def init_tensors(self):
        tensors_temp = list()
        for i in range(self.nsites):
            if i == 0:
                tensors_temp.append(zeros((self.sitedim,self.M[0])))
            elif i == self.nsites:
                tensors_temp.append(zeros((self.M[self.nbonds],self.sitedim)))
            else:
                tensors_temp.append(zeros((self.M[i-1],self.sitedim,self.M[0])))
        self.tensors = tensors_temp
        
    def init_guess_rand(self, ampl):
        for i,T in enumerate(self.tensors):
            self.tensors[i] = ampl*(np.random.random(T.shape))
            
                
        
    