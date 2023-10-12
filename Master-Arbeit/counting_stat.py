# /bin/env/python
from cProfile import label
from calendar import c
from ctypes import c_char_p
from email import message_from_file
import imp
from tkinter import messagebox
from turtle import color, title
from typing import List
from IPython.display import display
from re import A, U
from sys import displayhook
import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit, Aer, transpile, assemble
from qiskit.visualization import plot_histogram
from math import gcd
from numpy.random import randint
import pandas as pd
from fractions import Fraction
import qutip
print("conda activate")
from qutip import mesolve as mesolve
from qutip import basis as basis
print("import succesful")
from qutip import tensor as tensor
from qutip import dag as dag
from qutip import dims as dims

#from qutip import inv as inv
from qutip import steadystate as steadystate
from qutip import *
from qutip import ptrace 
import csv
from Loup_for_different_coupling import Diverse_Loups as Diverse_Loups
import multiprocessing as mp
import csv
from numpy.linalg import matrix_rank
#Konstante Grössen
##############################
# ##########################################################################
qutip.settings.has_mkl = False
omega_1=0
omega_2=30
omega_3=150

omega_f= omega_2 - omega_1
omega_h=omega_3-omega_1  # frequency of the atomic transition that is coupled to the hot bath
omega_c=omega_3-omega_2

omega_d=30

h=1
nph=9  # Maximale Photonen im cavity 
 
Th=100.    # temperature of the hot bath
Tc=20.     # temperature of the cold bath
Tenv=0.0000000000000000000000000001



nh=10
nc=0.02

nf=0.02    #Beschreibt den cavity/Photonen. 


f =0.1
global kappa

gamma_h=1
gamma_c=1
kappa=0.2
#kappa=0.028
kb=1
global g
g=14*kappa
#g=f=0


b_fock=qutip.states.fock(nph,0) #m)/fock(N,#m)
b_atom=basis(3)
b_comp=tensor( b_atom, b_fock)

#rho2=tensor(basis(6,4),basis(6,4).dag())


# hier ist ein wenig gebastel mit den transitionsoperatoren

va=qutip.Qobj(qutip.qutrit_basis()[2])
vb=qutip.Qobj(qutip.qutrit_basis()[1])
vg=qutip.Qobj(qutip.qutrit_basis()[0])

Trans_13=tensor(vg*va.dag(),qutip.identity(nph))
Trans_23=tensor(vb*va.dag(),qutip.identity(nph))
Trans_12=tensor(vg*vb.dag(),qutip.identity(nph))

proj_1=tensor(vg*vg.dag(),qutip.identity(nph))
proj_2=tensor(vb*vb.dag(),qutip.identity(nph))
proj_3=tensor(va*va.dag(),qutip.identity(nph))

a=qutip.tensor(qutip.identity(3),qutip.destroy(nph))




#print(qutip.to_super(np.sqrt(kappa_6)*A5))
#L=qutip.to_super(np.sqrt(kappa_6)*A5)


def Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g):
    H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a

    H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

    V=f*a.dag()+f*a #das got glaub nid

    H=H_free+H_int 

    Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)

    return Hdilde
Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g)
global DichteMatrix
def DichteMatrix(nh, nc, nf, Hami):
    gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
    gamma_2=(nh)*gamma_h
    gamma_3=(nc+1)*gamma_c
    gamma_4=(nc)*gamma_c
    kappa_5=(nf+1)*2*kappa ####goes to zero
    kappa_6=(nf)*2*kappa

    A1=Trans_13
    A2=Trans_13.dag()
    A3=Trans_23
    A4=Trans_23.dag()
    A5=a
    A6=a.dag()
########################################################################################################
      
    c_op_list=[]

    c_op_list.append(np.sqrt(gamma_1)*A1)
    c_op_list.append(np.sqrt(gamma_2)*A2)
    c_op_list.append(np.sqrt(gamma_3)*A3)
    c_op_list.append(np.sqrt(gamma_4)*A4)
    c_op_list.append(np.sqrt(kappa_5)*A5)
    c_op_list.append(np.sqrt(kappa_6)*A6)
    rho = steadystate(Hami, c_op_list)
    return rho


global colaps
def colaps(nh, nc, nf):

    gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
    gamma_2=(nh)*gamma_h
    gamma_3=(nc+1)*gamma_c
    gamma_4=(nc)*gamma_c
    kappa_5=(nf+1)*2*kappa ####goes to zero
    kappa_6=(nf)*2*kappa

    A1=Trans_13
    A2=Trans_13.dag()
    A3=Trans_23
    A4=Trans_23.dag()
    A5=a
    A6=a.dag()
########################################################################################################
      
    c_op_list=[]

    c_op_list.append(np.sqrt(gamma_1)*A1)
    c_op_list.append(np.sqrt(gamma_2)*A2)
    c_op_list.append(np.sqrt(gamma_3)*A3)
    c_op_list.append(np.sqrt(gamma_4)*A4)
    c_op_list.append(np.sqrt(kappa_5)*A5)
    c_op_list.append(np.sqrt(kappa_6)*A6)

    return c_op_list





vk_list=[-1,1]
#print("J==========",qutip.to_super(c_op_list[0]))
def J_sup(nh, nc, nf,vk_list):
    J=[]
    c_op_list=colaps(nh, nc, nf)
    J=vk_list[0]*(qutip.to_super(c_op_list[0]))+vk_list[1]*(qutip.to_super(c_op_list[1]))
    #J = [qutip.spre(c_op_list[0]) + qutip.spost(c_op_list[0].dag())]
    #J=np.array(J)
    return J

#print("J===",J_sup(nh, nc, nf,vk_list))
def K_trace(nh, nc, nf,vk_list):
    K=[]
    rhoV=qutip.operator_to_vector(DichteMatrix(nh, nc, nf, Hdilde))
    c_op_list=colaps(nh, nc, nf)
    K=(vk_list[0]**2)*(np.trace(qutip.to_super(c_op_list[0])*rhoV)+vk_list[1]**2)*(np.trace(qutip.to_super(c_op_list[1])*rhoV)) 
    
    return K



print (J_sup(nh, nc, nf,vk_list),"ktrace==",K_trace(nh, nc, nf,vk_list))




from scipy.sparse.linalg import dsolve
"""""

def solveZ2( nh,nc,nf, vk_list,Hdilde): #Z = the Drazin inverse times  alpha
    
    rho=DichteMatrix(nh, nc, nf, Hdilde)
    L=qutip.liouvillian(Hdilde,rho)
    rhoV=qutip.operator_to_vector(rho)
    rhoV
    IdV=(qutip.operator_to_vector(tensor(qutip.identity(3),qutip.identity(nph))))#bereits transponiert
    
    #O=(qutip.operator_to_vector(tensor(qutip.Qobj(np.zeros(s1)),qutip.Qobj(np.zeros(s2)))))
    O=np.ones(nph^2+1)
    I=np.array(IdV)
   
    
    alpha=J_sup(nh, nc, nf,vk_list)*rhoV
    alpha=np.real(alpha)
    #print(((alpha-rhoV*IdV.trans()*alpha)))
    Rechts=alpha-rhoV*IdV.dag()*alpha
    
    #RechtsForcea=np.array(Rechts.full())
    RechtsForce=np.vstack((Rechts.data,[[0]]))
    #LN=np.matrix(L.full())
    IdVD=IdV.data
    LForce=np.vstack((L.data,IdVD.trans()))
    
    print("Force=====",LForce)
    
    
    
    
    #Z = scipy.linalg.solve_triangular(L,Rechts) #uf der rechte site muess e vektor  stoh
    
    #Zs=np.linalg.lstsq(LForce,RechtsForce,rcond=None)
    Z2=dsolve.spsolve(L.data,Rechts.data,use_umfpack=False)
    Z2=Z2-rhoV.data*(IdVD.trans()+Z2)
    Z=qutip.Qobj(np.array(Z2)[0],dims=[[[3, nph], [3, nph]], [1]])
    return Z

print("mis z isch =========== lieg do ",solveZ2( nh,nc,nf, vk_list,Hdilde))


"""""
def solveZ( nh,nc,nf, vk_list,Hdilde): #Z = the Drazin inverse times  alpha
    
    rho=DichteMatrix(nh, nc, nf, Hdilde)
    L=qutip.liouvillian(Hdilde,rho)
    rhoV=qutip.operator_to_vector(rho)
    global IdV
    IdV=(qutip.operator_to_vector(tensor(qutip.identity(3),qutip.identity(nph))))#bereits transponiert


    alpha=J_sup(nh, nc, nf,vk_list)*rhoV
        
    #print(((alpha-rhoV*IdV.trans()*alpha)))
    Rechts=alpha-rhoV*IdV.trans()*alpha
    Rechts=np.real(Rechts)
    #RechtsForcea=np.array(Rechts.full())
    RechtsForce=np.vstack((Rechts,[[0]]))
    #LN=np.matrix(L.full())
    LForce=np.r_[L.full(),IdV.trans()]
    
    #print("Force=====",LForce)
    
   
    
    
    Zs=np.linalg.lstsq(LForce,RechtsForce,rcond=None)
    
    Z=qutip.Qobj(Zs[0],dims=[[[3, nph], [3, nph]], [1]])
    return Z

#print("Z====",solveZ(nh,nc,nf, vk_list,Hdilde))


"""with open('LLLL.csv', 'w', newline='') as yfile:
        writer = csv.writer(file)
        writer.writerows((LForc))"""
    


def averrageJ(nh,nc,nf, vk_list,Hdilde):


    rho=DichteMatrix(nh, nc, nf, Hdilde)
    L=qutip.liouvillian(Hdilde,rho)
    rhoV=qutip.operator_to_vector(rho)
    
    IdV=(qutip.operator_to_vector(tensor(qutip.identity(3),qutip.identity(nph))))#bereits transponiert
    alpha=J_sup(nh, nc, nf,vk_list)*rhoV
    alpha=np.real(alpha)
    averrageJ=IdV.trans()*alpha
    return averrageJ

def Dcalc(nh,nc,nf, vk_list,Hdilde):
    IdV=(qutip.operator_to_vector(tensor(qutip.identity(3),qutip.identity(nph))))

    D= K_trace(nh, nc, nf,vk_list)-2*(IdV.dag()*J_sup(nh,nc,nf, vk_list)*solveZ(nh,nc,nf, vk_list,Hdilde))
    return(np.real(D))

#print(np.imag( D(nh,nc,nf, vk_list,Hdilde))+np.imag(D(nh,nc,nf,vk_list,Hdilde)))
print(( np.real(Dcalc(nh,nc,nf, vk_list,Hdilde))))

#Ja = np.array([np.real((IdV.trans() * L1i * rhovec)[0,0]) for L1i in L1])
import scipy as scipy
rho=DichteMatrix(nh, nc, nf, Hdilde)
matrix=qutip.liouvillian(Hdilde,rho)

def compute_drazin_inverse(matrix):
    # Step 1: Compute the eigenvalues and eigenvectors
    eigenvalues, eigenvectors = matrix.eigenstates()

    # Step 2: Identify the blocks corresponding to eigenvalue 0
    zero_eigenvalue_indices = [i for i, eigenvalue in enumerate(eigenvalues) if eigenvalue == 0]
    zero_blocks = []

    # Step 3: Compute the matrix exponentials for each block and form the Drazin inverse
    drazin_inverse = qutip.Qobj(np.zeros_like(matrix))

    for index in zero_eigenvalue_indices:
        block = eigenvectors[index][0]
        block_exp = scipy.linalg.expm(block)
        block_inverse =np.linalg.inv(block)
        #block_inverse =block.qutip.inv()
        block= qutip.Qobj(block)
        drazin_block = block_inverse - qutip.qeye(1)
        drazin_inverse += block * drazin_block * block_exp

    return drazin_inverse




D_list=[]
nh_list=[]
anzahl=30
nh=0
step=0.5
for i in range(anzahl):
    D=Dcalc(nh,nc,nf, vk_list,Hdilde)
    D_list.append((np.trace(D)))
    print("bla==rho",np.trace(D))


    nh_list.append(nh)
    nh+=step


fig3, ax = plt.subplots()

ax.set_xlabel(r' $n_h$', fontsize=19)
ax.set_ylabel('D')
plt.title(r' D vs $n_h$ ')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(D_list)[:anzahl],label=r' $\langle{J_h}\rangle$',color='red')

#plt.plot(np.asarray(nh_list3)[:100],np.asarray(Entropy2)[:100,3],'--',label=r' $\frac{J_{cav}}{T_{cav}}$',color='orange')
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')
plt.show()#
print(K_trace(100, nc, nf, Hdilde))













from sympy import Matrix
#j=matrix.jordan_form()
def Drazin(matrix):
    
    A=np.matrix(matrix)
    
    Ak=A
    k=0
    Aold=A
    Ad=0
    Ad1=0
   
    while matrix_rank(Ak) != matrix_rank(Ak*A):
            
        
        k+=1
        Ak=Ak*A
        
    

    
    if matrix_rank(Ak) == matrix_rank(Ak*A):
        
        Ad=np.linalg.lstsq(Ak*A,Ak,rcond=None)
        #Ad1=np.linalg.solve(Ak*A,Ak)
        
                
    else :
        print("error")
        
    
    Ad=qutip.Qobj(Ad[0],dims = [[[3, nph], [3, nph]], [[3, nph], [3, nph]]])

    return Ad






#matrix = qutip.Qobj([[1, 2, 3], [0, 0, 4], [4, 0, 1]])

matrix = [[1, 2, 3], [0, 0, 4], [4, 0, 1]]
L=qutip.liouvillian(Hdilde,rho)#andere liouvillian neh
L.data
L=L.full()
print(L)
IdV=(qutip.operator_to_vector(tensor(qutip.identity(3),qutip.identity(nph))))
drazin_inv = Drazin(L)
#drazin_inv = Drazin(matrix)
print("Drazin inverse:\n", drazin_inv)
Lstrich=J_sup(nh,nc,nf,vk_list)
D=drazin_inv

with open('LLLL.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows((Lstrich*D*Lstrich))

d=IdV.trans()*Lstrich*D*Lstrich*qutip.to_super(rho)*qutip.operator_to_vector(rho)

print(d)
#L=qutip.liouvillian(Hdilde,rho)
#print(steadystate(L))


















def tilted_liouvillian(H, L, chi, v):
    
    # Only works for one jump operator
    H_vec = spre(H) - spost(H)
    L_vec = np.exp(1j*chi*v)*to_super(L) - 0.5*(spre(L.dag() * L) + spost(L.dag()*L))
    
    return -1j*H_vec + L_vec


# Compute vectorised density operator
rhovec = operator_to_vector(rho)

# Create chi space
chi = np.linspace(-np.pi, np.pi, 100)
dchi = chi[1]-chi[0]

t = [1, 10, 20, 30]

# Compute tilted Liouvillians
c_op_list=colaps(nh, nc, nf)
liouvs = [tilted_liouvillian(Hdilde, c_op_list[0], i, 1) for i in chi]
rhochi = [[vector_to_operator((liou * ti).expm()*rhovec) for liou in liouvs] for ti in t]

# Take trace of all rhochis
pchis = np.array([[rho.tr() for rho in rhoix] for rhoix in rhochi])


