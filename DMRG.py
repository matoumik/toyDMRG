from numpy import array, zeros
import numpy as np
from H_Heis import H_heis
from scipy.sparse.linalg import eigsh

class DMRG:
    def __init__(self, DMRG_info):
        self.info = DMRG_info
        self.iter = 0
        self.sweep = 0
        self.nleft = 0
        self.nright = self.info.nsites
        self.MPS = MPS(self.info)
        self.H = H_heis(self.info.nsites)
        
    def H_contr_L(self, ncontr):
        for i in range(ncontr):
            if i == 0:
                A=self.MPS.tensors[0]
                #print(A.shape)
                A=np.tensordot(A,self.H.tensors[0],(0,0))
                #print(A.shape)
                A=np.tensordot(A,self.MPS.tensors[0],(1,0))
                #print(A.shape)
            else:
                #print(i,":")
                A=np.tensordot(A,self.MPS.tensors[i],(0,0))
                #print(A.shape)
                A=np.tensordot(A,self.H.tensors[i],((0,2),(0,1)))
                #print(A.shape)
                A=np.tensordot(A,self.MPS.tensors[i],((0,2),(0,1)))
                #print(A.shape)
                
        return A
                
    def H_contr_R(self, ncontr):
        a =0
        for i in range(self.info.nsites-1,self.info.nsites-ncontr-1,-1):
            #print(i,":")
            if i-self.info.nsites == -1:
                A=self.MPS.tensors[i]
                #print(A.shape)
                A=np.tensordot(A,self.H.tensors[i],(1,1))
                #print(A.shape)
                A=np.tensordot(A,self.MPS.tensors[i],(2,1))
                #print(A.shape)
            else:
                A=np.tensordot(A,self.MPS.tensors[i],(0,2))
                #print(A.shape)
                A=np.tensordot(A,self.H.tensors[i],((3,0),(1,3)))
                #print(A.shape)
                A=np.tensordot(A,self.MPS.tensors[i],((0,3),(2,1)))
                #print(A.shape)
            a+=1
        #print("contracted ",a," sites")
        return A
    
    def H_onesite(self, sitenum):
        if sitenum == 0:
            A = np.tensordot(self.H.tensors[0],self.H_contr_R(self.info.nsites-1),(2,1))
        elif sitenum == self.info.nsites - 1:
            A = np.tensordot(self.H_contr_L(self.info.nsites -1 ),self.H.tensors[self.info.nsites - 1], (1,0))
        else:
            A = np.tensordot(self.H_contr_L(sitenum),self.H.tensors[sitenum], (1,0))
            #print(A.shape)
            A = np.tensordot(A, self.H_contr_R(self.info.nsites - sitenum - 1), (4,1))
        
        #print(A.shape)
            
        return A
            
    def H_onesite_mat(self, sitenum):
        H = self.H_onesite(sitenum)
        if sitenum == 0:
            H = H.transpose(0,2,1,3)
            H = H.reshape(H.shape[0]*H.shape[1],-1)
        elif sitenum == self.info.nsites - 1:
            H = H.transpose(0,2,1,3)
            H = H.reshape(H.shape[0]*H.shape[1],-1)
        else:
            H = H.transpose(0,2,4,1,3,5)
            #print(H.shape)
            H = H.reshape(H.shape[0]*H.shape[1]*H.shape[2],-1)
        
        #print(H.shape)
        return H
    
    def solve_onesite(self, sitenum):
        Psi = self.MPS.tensors[sitenum]
        indexes = Psi.shape
        H = self.H_onesite_mat(sitenum)
        Psi.reshape(-1)
        eigs, eigvals = eigsh(H,k=3, which = 'SA', v0=Psi, tol = 1e-6)
        print(eigs)
        #print(eigvals.shape)
        Psi = eigvals[:,0].reshape(indexes)
        self.MPS.tensors[sitenum] = Psi
        
    def sweep_half(self):
        for i in range(self.info.nsites//2, self.info.nsites):
            self.solve_onesite(i)
    
    def sweep_forw(self):
        for i in range(self.info.nsites):
            self.solve_onesite(i)
            
    def sweep_backw(self):
        for i in range(self.info.nsites-1, 0, -1):
            self.solve_onesite(i)
            
    def do_sweep(self):
        self.sweep_forw()
        self.sweep_backw()


class DMRG_info:
    def __init__(self):
        self.M_start = 100
        self.nsites = 10


class MPS:
    def __init__(self, DMRG_info):
        self.nsites = DMRG_info.nsites
        self.sitedim = 2
        self.nbonds = self.nsites - 1
        self.init_M(DMRG_info.M_start)
        self.init_tensors()
        self.DBSS= False #TODO      
        

    def init_M(self, Mstart):
        M_temp = list()
        k = self.nbonds
        for i in range(self.nbonds):
            if (self.sitedim**(i+1)<Mstart and self.sitedim**(i+1)<self.sitedim**(k-i)):
                M_temp.append(self.sitedim**(i+1))
            elif self.sitedim**(k-i)<Mstart:
                M_temp.append(self.sitedim**(k-i))
            else:
                M_temp.append(Mstart)
        self.M = array(M_temp)
        
    def init_tensors(self):
        tensors_temp = list()
        for i in range(self.nsites):
            if i == 0:
                tensors_temp.append(zeros((self.sitedim,self.M[0])))
            elif i == self.nsites-1:
                tensors_temp.append(zeros((self.M[-1],self.sitedim)))
            else:
                tensors_temp.append(zeros((self.M[i-1],self.sitedim,self.M[i])))
        self.tensors = tensors_temp
        
    def init_guess_rand(self, ampl):
        for i,T in enumerate(self.tensors):
            self.tensors[i] = ampl*(np.random.random(T.shape))
            
                
        
    