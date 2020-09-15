from H import H
from numpy import array

Sz = array([[0.5,0],[0,-0.5]])
Sp = array([[0,1],[0,0]])
Sm = array([[0,0],[1,0]])
O = array([[0,0],[0,0]])
I = array([[1,0],[0,1]])

j=1.0





print(Sp.view())

class H_heis(H):
    def __init__(self,nsites):
        tensors_temp = list()
        for i in range(nsites):
            if i == 0:
                A = array([O,0.5*j*Sm, 0.5*j*Sp, j*Sz, I])
                tensors_temp.append(A.transpose(1,2,0))
            elif i == nsites-1:
                A = array([I,Sp,Sm,Sz,O])
                tensors_temp.append(A)
                
            else: 
                A=array([[I,O,O,O,O],[Sp,O,O,O,O],[Sm,O,O,O,O],[Sz,O,O,O,O],[O,0.5*j*Sm,0.5*j*Sp, j*Sz, I]])
                tensors_temp.append(A.transpose(0,2,3,1))
            self.tensors = tensors_temp
            
H = H_heis(4)

MPO = H.tensors