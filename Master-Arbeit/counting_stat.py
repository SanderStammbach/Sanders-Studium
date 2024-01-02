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
from Loup_for_different_coupling import Diverse_Loups
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
nph=10 # Maximale Photonen im cavity for fsc it works fine with 19
 
Th=100.    # temperature of the hot bath
Tc=20.     # temperature of the cold bath
Tenv=0.0000000000000000000000000001


nH=5

nh_fix=nh=0.0013

nc_fix=nc=0.027

nf=nf_fix=0.07 #Beschreibt den cavity/Photonen. 

f_fix=f =0.0


gamma_h=0.1
gamma_c=2
#kappa=0.2
#kappa=0.028
kb=1

g_fix=g=2.8
#g=14*kappa
#g=f=0
#kappa=0.07
kappa_fix=kappa=0.002

epsilon=0.15













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
H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a
#H_free=1*proj_1+1*proj_2+1*proj_3+1*(a.dag()*a)


#print(qutip.to_super(np.sqrt(kappa_6)*A5))
#L=qutip.to_super(np.sqrt(kappa_6)*A5)


def Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d):
    H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a

    H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

    V=f*a.dag()+f*a #das got glaub nid

    H=H_free+H_int 

    Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)
    
    #Hdilde = f*(Trans_12+Trans_12.dag())

    return Hdilde

def DichteMatrix(nh, nc, nf, Hami,kappa):
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


def colaps(nh, nc, nf,kappa):
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


B=1
A=0

""""
vk_list=[-1,1]
#print("J==========",qutip.to_super(c_op_list[0]))
def J_sup(nh, nc, nf,vk_list,collaps_list):
    J=[]
    c_op_list=collaps_list
    J=vk_list[0]*(qutip.spre(collaps_list[A])+qutip.spost(collaps_list[A].dag()))+vk_list[1]*(qutip.spre(collaps_list[B])+qutip.spost(collaps_list[B].dag()))
    #J = [qutip.spre(c_op_list[0]) + qutip.spost(c_op_list[0].dag())]
    #J=np.array(J)
    return J
"""
vk_list=[-1,1]
#print("J==========",qutip.to_super(c_op_list[0]))
def J_sup(nh, nc, nf,vk_list,collaps_list):
    J=[]
    c_op_list=collaps_list
    J=vk_list[0]*(qutip.to_super(collaps_list[A]))+vk_list[1]*(qutip.to_super(collaps_list[B]))
    #J = [qutip.spre(c_op_list[0]) + qutip.spost(c_op_list[0].dag())]
    #J=np.array(J)
    return J


"""
def Entropy(nh, nc, nf,vk_list,collaps_list,omega):
    def T(n,omega):
            
            T=h*omega/(kb*(np.log((1/n)+1)))
            return T
        Entropy=J=vk_list[0]*(qutip.to_super(collaps_list[A]))+vk_list[1]*(qutip.to_super(collaps_list[B]))

    return Entropy
    
"""
#print("J===",J_sup(nh, nc, nf,vk_list))
def K_trace(nh, nc, nf,vk_list,Hdilde): 
    K=[]
    rhoV=qutip.operator_to_vector(DichteMatrix(nh, nc, nf, Hdilde,kappa))
    c_op_list=colaps(nh, nc, nf,kappa)
    K=np.trace((vk_list[0]**2)*(qutip.to_super(c_op_list[A])*rhoV)+(vk_list[1]**2)*(qutip.to_super(c_op_list[B])*rhoV))
    
    return np.real(K)

def K_Trace2(nh, nc, nf,vk_list,Hdilde,collaps_list): 
    K=[]
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa)
    c_op_list=colaps(nh, nc, nf,kappa)
    K=np.trace((vk_list[0]**2)*((c_op_list[A].dag()*c_op_list[A])*rho)+(vk_list[1]**2)*((c_op_list[B].dag()*c_op_list[B])*rho))
    
    return np.real(K)

#print("K===========",K_trace(nh, nc, nf,vk_list))








def EnergieCalculator_mit_faktor(c_op_list ,g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h):
   
    V=f*a.dag()+f*a
    H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())
    omega_1=0
       
    H_free1=H_free 
    Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)   
    rho = steadystate(Hdilde, c_op_list) ######## Are you sure its with only photons H_free?
    rho_f=rho.ptrace(1)
        
    def D(c_op_list,rho):
        D=[]
        for i in range(6):
                
            D.append(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho+rho*c_op_list[i].dag()*c_op_list[i]))
            
        return D

    Liste_von_Q=[] # ExpectValue for Thermal Energy

    Liste_von_Q.append(np.trace(H_free1*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1]))/omega_h)
    Liste_von_Q.append(np.trace(H_free1*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3]))/omega_h)
    Liste_von_Q.append(np.trace(H_free1*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5]))/omega_h)
    #Liste_von_Q.append(g)  g in der liste anfügen

    float_list= list(np.float_(Liste_von_Q))
           
    Liste_von_Q=float_list

    return(Liste_von_Q)





def Entropy_ohne_omega(c_op_list,nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3):
        

        H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a
    
        
        V=f*a.dag()+f*a
        H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())
        
        
        H_free1=H_free 
        Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)      
        rho = steadystate(Hdilde, c_op_list) ######## Are you sure its with only photons H_free?
        rho_f=rho.ptrace(1)

        def D(c_op_list,rho):
            D=[]
            for i in range(6):
                
                D.append(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho+rho*c_op_list[i].dag()*c_op_list[i]))
            
            return D
        
        def T(n,omega):
            
            T=h*omega/(kb*(np.log((1/n)+1)))
            return T
        
        omega_c=120
        omega_f=30
        omega_h=150

        
        Liste_von_Q=[] # ExpectValue for Thermal Energy
     
        Liste_von_Q.append((np.trace(-H_free1*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1])))/(T(nh,omega_h)))
        Liste_von_Q.append((np.trace(-H_free1*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3])))/(T(nc,omega_c)))
        Liste_von_Q.append((np.trace(-H_free1*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5])))/(T(nf,omega_f)))
        Liste_von_Q.append((np.trace(-H_free1*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1])))/(T(nh,omega_h))+(np.trace(-H_free1*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3])))/(T(nc,omega_c))+(np.trace(-H_free1*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5])))/(T(nf,omega_f)))
        Liste_von_Q.append((np.trace(-H_free1*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1])))/(T(nh,omega_h))+(np.trace(-H_free1*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3])))/(T(nc,omega_c)))
        print("T=",T(nc,omega_c))
        #Liste_von_Q.append(g)  g in der liste anfügen

        float_list= list(np.float_(Liste_von_Q))
            
        Liste_von_Q=float_list

        return(Liste_von_Q)

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
    c_op_list=colaps(nh, nc, nf,kappa)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa)
    L=qutip.liouvillian(Hdilde,c_op_list)
    rhoV=qutip.operator_to_vector(rho)
    global IdV
    IdV=(qutip.operator_to_vector(tensor(qutip.identity(3),qutip.identity(nph))))#bereits transponiert


    alpha=J_sup(nh, nc, nf,vk_list,c_op_list)*rhoV
        
    #print(((alpha-rhoV*IdV.trans()*alpha)))
    Rechts=alpha-rhoV*IdV.trans()*alpha
    #Rechts=np.real(Rechts)
    
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
    beta=J_sup(nh, nc, nf,vk_list,collaps_list)*rhoV
    beta=np.real(beta)
    averrageJ=IdV.trans()*beta
    return averrageJ


"""
#print(np.imag( D(nh,nc,nf, vk_list,Hdilde))+np.imag(D(nh,nc,nf,vk_list,Hdilde)))
#print(( np.real(Dcalc(nh,nc,nf, vk_list,Hdilde))))

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
"""

def Dcalc(nh,nc,nf, vk_list,Hdilde):
    collaps_list=colaps(nh, nc, nf,kappa)
    IdV=(qutip.operator_to_vector(tensor(qutip.identity(3),qutip.identity(nph))))
    print(K_trace(nh,nc,nf,vk_list,Hdilde))
    D= K_Trace2(nh, nc, nf,vk_list,Hdilde,collaps_list)-2*np.real(IdV.dag()*J_sup(nh,nc,nf, vk_list,collaps_list)*solveZ(nh,nc,nf, vk_list,Hdilde))
    return(D)


def Potts(nh,nc,gammac,gammah,epsilon):
    Q=(((1 + nc)*nh + nc*(1 + nh))*np.log(((1 + nc)*nh)/(nc*(1 + nh))))/(-nc + nh) - (2*epsilon**2*gammac*gammah*(-nc + nh)*(gammac*nc + gammah*nh)*(-2/(gammac*nc + gammah*nh) + (gammac*gammah*(gammac*nc + gammah*nh)*(nc + nh + 3*nc*nh) + (gammac + gammah + 2*(gammac*nc + gammah*nh))*(4*epsilon**2 + (gammac*nc + gammah*nh)**2/4.))/ ((gammac*gammah*(gammac*nc + gammah*nh)**2*(nc + nh + 3*nc*nh))/4. + 2*epsilon**2*(gammac*nc + gammah*nh)*(gammac + gammah + (3*(gammac*nc + gammah*nh))/2.)))*np.log(((1 + nc)*nh)/(nc*(1 + nh))))/ ((gammac*gammah*(gammac*nc + gammah*nh)**2*(nc + nh + 3*nc*nh))/4. + 2*epsilon**2*(gammac*nc + gammah*nh)*(gammac + gammah + (3*(gammac*nc + gammah*nh))/2.))
    return Q


HdildeTest = Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,0.3,g,omega_d)
rhoTest=DichteMatrix(nh,nc,nf,HdildeTest,kappa)
if rhoTest.ptrace(1)[[nph-1], [nph-1]]>10**(-4):
    print("Error",rhoTest.ptrace(1)[[nph-1], [nph-1]])
    
        
anzahl=100
step =0.0018
step=0.2/anzahl
nc_Line=np.linspace(0, 0.6, 100, endpoint=True)
Q=Potts(nh,nc_Line,gamma_c,gamma_h,epsilon)

nc=0.0001

D_list=[]
nh_list=[]
Jh_list1=[]
Entropy_list1=[]
Q_list1=[]
Energy_list1=[]
Q_List_PRE=[]
for i in range(anzahl):
    collapse_list=colaps(nh, nc, nf,kappa)
    Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa)
    D=Dcalc(nh,nc,nf, vk_list,Hdilde)
    D_list.append(D[0])
    print("bla==rho",np.trace(D))
    Lstrich=J_sup(nh, nc, nf,vk_list,collapse_list)
    Jh_list1.append(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])#noesomega h  dezue
    print(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])
    Entropy_list1.append(Entropy_ohne_omega(collapse_list,nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3)[3])
   
    Energy_list1.append(EnergieCalculator_mit_faktor(collapse_list,g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h))
    Q_list1.append((Entropy_list1[i])*(D_list[i]/(Jh_list1[i]**2)))
    Q_List_PRE.append(Potts(nH,nc,gamma_c,gamma_h,epsilon))
    print(Q_List_PRE[i])
    #Q_list1.append((Entropy_list1[i])*(D_list[i]/(Jh_list1[i]**2)))
    nh_list.append(nc)
    nc+=step
    print(i-anzahl, D_list[i])







    

g=g_fix
nh=nh_fix
f=0.02
omega_f=30
omega_d=29
Delta_list=[]
step=2/anzahl

D_list_D=[]
D_list_D=[]
Jh_list_D=[]
Entropy_list_D=[]
Q_list_D=[]
Energy_list_D=[]
for i in range(anzahl):
    
    collapse_list=colaps(nh, nc, nf,kappa)
    Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa)
    D=Dcalc(nh,nc,nf, vk_list,Hdilde)
    D_list_D.append(D[0])
    print("bla==rho",np.trace(D))
    Lstrich=J_sup(nh, nc, nf,vk_list,collapse_list)
    Jh_list_D.append(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])#noesomega h  dezue
    print(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])
    Energy_list_D.append(EnergieCalculator_mit_faktor(collapse_list,g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h))
    Entropy_list_D.append(Entropy_ohne_omega(collapse_list,nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3)[3])
    Q_list_D.append((Entropy_list_D[i])*(D_list_D[i]/(Energy_list_D[i][0]**2)))

    Delta_list.append(omega_d-omega_f)
    omega_d+=step
    print(i-anzahl)
    
  
g=g_fix   
f=f_fix 
nh=nh_fix
omega_d=30   
D_list_nf=[]
nf_list=[]
Jh_list_nf=[]
Entropy_list_nf=[]
Q_list_nf=[]
Energy_list_nf=[]
step=0.5/anzahl
nf=0.00001
for i in range(anzahl):
    g=g_fix   
    f=f_fix 
    nh=nh_fix
    nc=nc_fix
    omega_d=30  
    collapse_list=colaps(nh, nc, nf,kappa)
    Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa)
    D=Dcalc(nh,nc,nf, vk_list,Hdilde)
    D_list_nf.append(D[0])
    print("bla==rho",np.trace(D))
    Lstrich=J_sup(nh, nc, nf,vk_list,collapse_list)
    Jh_list_nf.append(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])#noesomega h  dezue
    print(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])
    Entropy_list_nf.append(Entropy_ohne_omega(collapse_list,nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3)[3])
   
   
    Energy_list_nf.append(EnergieCalculator_mit_faktor(collapse_list,g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h))
    Q_list_nf.append((Entropy_list_nf[i])*(D_list_nf[i]/(Jh_list_nf[i]**2)))
    nf_list.append(nf)
    nf+=step
    print(i-anzahl, D_list[i],"parameter=",nh,nc,nf,g,kappa,f)


label1=r' $Var \langle J_{h} \rangle $'
label12=r' $\langle J_{h} \rangle $'
label3=r' $\mathcal{Q}'
label4=r' $\dot{\sigma}$'
fig, (ax1,ax2,ax4) = plt.subplots(3)

ax1.set_xlabel(r' $n_c$', fontsize=19)
ax1.set_ylabel('D')
#plt.title(r' D vs $n_h$ ')
#ax1.plot(np.asarray(nc_Line)[:anzahl],np.asarray(Q)[:anzahl],'*',label=r' $\mathcal{Q}$',color='blue',alpha=0.5)
ax1.plot(np.asarray(nh_list)[:anzahl],np.asarray(D_list)[:anzahl],label=label1,color='black')
ax1.plot(np.asarray(nh_list)[:anzahl],np.asarray(Jh_list1)[:anzahl],label=label12,color='red')
ax1.plot(np.asarray(nh_list)[:anzahl],np.asarray(Q_list1)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
ax1.plot(np.asarray(nh_list)[:anzahl],np.asarray(Q_List_PRE)[:anzahl],'-',label=r' $\mathcal{Q} PRE $',color='blue',alpha=0.5)
ax1.plot(np.asarray(nh_list)[:anzahl],np.asarray(Entropy_list1)[:anzahl],'-',label=label4,color='orange')
legend = ax1.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')
ax1.axhline(y=2)
ax1.grid()



ax2.set_xlabel(r' $n_f$', fontsize=19)
ax2.set_ylabel(r' $\mathcal{Q}$')
#plt.title(r' $\mathcal{Q} vs f$ ')
ax2.plot(np.asarray(nf_list)[:anzahl],np.asarray(D_list_nf)[:anzahl],label=label1,color='black')
ax2.plot(np.asarray(nf_list)[:anzahl],np.asarray(Jh_list_nf)[:anzahl],label=label12,color='red')
ax2.plot(np.asarray(nf_list)[:anzahl],np.asarray(Q_list_nf)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
ax2.plot(np.asarray(nf_list)[:anzahl],np.asarray(Entropy_list_nf)[:anzahl],'-',label=r' $\dot{\sigma}$',color='orange')
legend = ax2.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')
ax2.axhline(y=2)
ax2.grid()






ax4.set_xlabel(r' $\Delta$', fontsize=19)
ax4.set_ylabel(r' $\mathcal{Q}$')
#plt.title(r' $\mathcal{Q} vs f$ ')
ax4.plot(np.asarray(Delta_list)[:anzahl],np.asarray(D_list_D)[:anzahl],label=label1,color='black')
ax4.plot(np.asarray(Delta_list)[:anzahl],np.asarray(Jh_list_D)[:anzahl],label=label12,color='red')
ax4.plot(np.asarray(Delta_list)[:anzahl],np.asarray(Q_list_D)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
ax4.plot(np.asarray(Delta_list)[:anzahl],np.asarray(Entropy_list_D)[:anzahl],'-',label=label4,color='orange')

ax4.axhline(y=2)
ax4.grid()
legend = plt.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')




plt.show()#



fig, ax= plt.subplots()

plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(Energy_list_nf)[:anzahl,0],'-',color='red',label=r' $\frac{J_{h}}{\hbar\gamma_h \omega_{h}}$')

plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(Jh_list_nf)[:anzahl],'*',label=r' $\langle \langle 1|\mathcal{J}|\rho \rangle \rangle $',color='red')

plt.show()











f=f_fix
nf=nf_fix
nc=nc_fix
nh=nh_fix
omega_d=30

kappa_list=[]
D_list_f=[]
Jh_list_kappa=[]
Entropy_list_kappa=[]
Q_list_f=[]


step=4/anzahl
Energy_list_kappa=[]
D_list_kappa=[]
Q_list_kappa=[]
kappa=0.000000001
for i in range(anzahl):
    collapse_list=colaps(nh, nc, nf,kappa)
    Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa)
    D=Dcalc(nh,nc,nf, vk_list,Hdilde)
    D_list_kappa.append(D[0])
    print("bla==rho",np.trace(D))
    Lstrich=J_sup(nh, nc, nf,vk_list,collapse_list)
    Jh_list_kappa.append(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])#noesomega h  dezue
    print(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])
    Energy_list_kappa.append(EnergieCalculator_mit_faktor(collapse_list ,g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h))
    Entropy_list_kappa.append(Entropy_ohne_omega(collapse_list,nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3)[3])
    Q_list_kappa.append((Entropy_list_kappa[i])*(D_list_kappa[i]/(Energy_list_kappa[i][0]**2)))

    kappa_list.append(kappa)
    kappa+=step
    print(i-anzahl,print(nh,nc,nf,f,g,kappa))

kappa=kappa_fix 



f=f_fix
nh=nh_fix
nc=nc_fix
g_list=[]
D_list_g=[]
Jh_list_g=[]
Entropy_list_g=[]
Q_list_g=[]
g=0
step=1/anzahl
Energy_list_g=[]
kappa=kappa_fix
for i in range(anzahl):
    collapse_list=colaps(nh, nc, nf,kappa)
    Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa)
    D=Dcalc(nh,nc,nf, vk_list,Hdilde)
    D_list_g.append(D[0])
    print("bla==rho",np.trace(D))
    Lstrich=J_sup(nh, nc, nf,vk_list,collapse_list)
    Jh_list_g.append(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])#noesomega h  dezue
    print(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])
    Energy_list_g.append(EnergieCalculator_mit_faktor(collapse_list,g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h))
    Entropy_list_g.append(Entropy_ohne_omega(collapse_list,nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3)[3])
    Q_list_g.append((Entropy_list_g[i])*(D_list_g[i]/(Energy_list_g[i][0]**2)))
    print("QList g=",Q_list_g[i])
    g_list.append(g)
    g+=step
    print(i-anzahl)    
    
    
    


nc=nc_fix
nh=nh_fix
f_list=[]
D_list_f=[]
Jh_list_f=[]
Entropy_list_f=[]
Q_list_f=[]
f=0.001
step=0.3/anzahl

Energy_list_f=[]

for i in range(anzahl):
    nf=0.5
    Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa)
    D=Dcalc(nh,nc,nf, vk_list,Hdilde)
    D_list_f.append(D[0])
    collapse_list=colaps(nh, nc, nf,kappa)
    print("bla==rho",np.trace(D))
    Lstrich=J_sup(nh, nc, nf,vk_list,collapse_list)
    Jh_list_f.append(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])#noesomega h  dezue
    print(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])
    Energy_list_f.append(EnergieCalculator_mit_faktor(collapse_list,g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h))
    Entropy_list_f.append(Entropy_ohne_omega(collapse_list,nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3)[3])
    Q_list_f.append((Entropy_list_f[i])*(D_list_f[i]/(Jh_list_f[i]**2)))

    f_list.append(f)
    f+=step
    print(i-anzahl)



fig2, (ax1,ax2,ax3) = plt.subplots(3)

ax1.set_xlabel(r' $\kappa$', fontsize=19)
ax1.set_ylabel(r' $\mathcal{Q}$')
#plt.title(r' D vs $n_h$ ')
ax1.plot(np.asarray(kappa_list)[:anzahl],np.asarray(D_list_kappa)[:anzahl],label=label1,color='black')
ax1.plot(np.asarray(kappa_list)[:anzahl],np.asarray(Jh_list1)[:anzahl],label=label12,color='red')
ax1.plot(np.asarray(kappa_list)[:anzahl],np.asarray(Q_list_kappa)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
ax1.plot(np.asarray(kappa_list)[:anzahl],np.asarray(Entropy_list_kappa)[:anzahl],'-',label=r' $\dot{\sigma}$',color='orange')
legend = ax1.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')
ax1.axhline(y=2)
ax1.grid()



ax2.set_xlabel(r' $f/\gamma_h$', fontsize=19)
ax2.set_ylabel(r' $\mathcal{Q}$')
#plt.title(r' $\mathcal{Q} vs f$ ')
ax2.plot(np.asarray(f_list)[:anzahl],np.asarray(D_list_f)[:anzahl],label=label1,color='black')
ax2.plot(np.asarray(f_list)[:anzahl],np.asarray(Jh_list_f)[:anzahl],label=label12,color='red')
ax2.plot(np.asarray(f_list)[:anzahl],np.asarray(Q_list_f)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
ax2.plot(np.asarray(f_list)[:anzahl],np.asarray(Entropy_list_f)[:anzahl],'-',label=label4,color='orange')
legend = ax2.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')





legend = plt.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')
ax2.axhline(y=2)
ax2.grid()

ax3.set_xlabel(r' $g$', fontsize=19)
ax3.set_ylabel(r' $\mathcal{Q}$')
#plt.title(r' $\mathcal{Q} vs f$ ')
ax3.plot(np.asarray(g_list)[:anzahl],np.asarray(D_list_g)[:anzahl],label=label1,color='black')
ax3.plot(np.asarray(g_list)[:anzahl],np.asarray(Jh_list_g)[:anzahl],label=label12,color='red')
ax3.plot(np.asarray(g_list)[:anzahl],np.asarray(Q_list_g)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
ax3.plot(np.asarray(g_list)[:anzahl],np.asarray(Entropy_list_g)[:anzahl],'-',label=label4,color='orange')

ax3.axhline(y=2)
ax3.grid()
legend = ax3.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')


plt.show()#


nf=nf_fix







































































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
        print(k,Ak)
        
    if matrix_rank(Ak) == matrix_rank(Ak*A):
        
        Ad=np.linalg.lstsq(Ak*A,Ak,rcond=None)
        #Ad1=np.linalg.solve(Ak*A,Ak)
        
                
    else :
        print("error")
        
    
    Ad=qutip.Qobj(Ad[0],dims = [[[3, nph], [3, nph]], [[3, nph], [3, nph]]])

    return Ad


from scipy.linalg import null_space

import numpy as np
from sympy import Matrix
from scipy.sparse.linalg import dsolve
from numpy.linalg import matrix_rank
import scipy  as scipy






def Drazin2(matrix):
    
    A=np.matrix(matrix)
    #A=matrix
    Ak=A
    k=0
    Aold=A
    Ad=0
    Ad1=0
   
    while matrix_rank(Ak) != matrix_rank(Ak*A):
            
        
        k+=1
        Ak=Ak*A
        print(k,Ak)
        
    if matrix_rank(Ak) == matrix_rank(Ak*A):
        
        NSu=scipy.linalg.null_space(Ak)
        NSv=np.transpose( scipy.linalg.null_space(np.transpose(Ak)))
        dimZ1, dimZ2= NSv.shape
        
        Zeros=np.zeros((dimZ1,dimZ1))

        D=np.block([[A,NSu],[NSv,Zeros]])
        D_inv=np.linalg.inv(D)

        
        num_rows, num_cols = A.shape
        Drazin_inv=D_inv[0:num_rows,0:num_cols]
                
    else :
        print("error")
        
    Ad=qutip.Qobj(Drazin_inv,dims = [[[3, nph], [3, nph]], [[3, nph], [3, nph]]])
    

    return Ad




""""

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






"""

#do bruchbar



IdV=(qutip.operator_to_vector(tensor(qutip.identity(3),qutip.identity(nph))))
f=0.02
omega_f=30
anzahl=30
step =0.1
nh=0.0001
nh_list=[]
d_list=[]
Jh_list2=[]
Entropy_list=[]
Q_list=[]
Energy_list1=[]
for i in range(anzahl):
    Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa)
    c_op_list=colaps(nh, nc, nf,kappa)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa)
    L=qutip.liouvillian(Hdilde,c_op_list)
    print(L)
    Lstrich=J_sup(nh,nc,nf,vk_list,c_op_list)
    d_list.append(K_Trace2(nh, nc, nf,vk_list,Hdilde,c_op_list)-2*np.real((IdV.trans()*Lstrich*Drazin2(L)*Lstrich*qutip.operator_to_vector(rho)))[0])
    Jh_list2.append(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])#noesomega h  dezue
    print(np.imag((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])
    Energy_list1.append(EnergieCalculator_mit_faktor(c_op_list,g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h))
    Entropy_list.append(Entropy_ohne_omega(c_op_list,nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3)[3])
    Q_list.append((Entropy_list[i])*(d_list[i]/(Energy_list1[i][0]**2)))
    
    nh_list.append(nh)
    nh=nh+step
    







fig3, ax = plt.subplots()

ax.set_xlabel(r' $n_h$', fontsize=19)
ax.set_ylabel('D')
plt.title(r' D vs $n_h$ ')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(d_list)[:anzahl],label=r' $Var \langle J \rangle $',color='black')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(Jh_list2)[:anzahl],label=r' $\langle \langle 1|\mathcal{J}|\rho \rangle \rangle $',color='red')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(Q_list)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
#plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(Entropy_list)[:anzahl],'-',label=r' $H_{free} \mathcal{L}_h $',color='purple')
#plt.plot(np.asarray(nh_list3)[:100],np.asarray(Entropy2)[:100,3],'--',label=r' $\frac{J_{cav}}{T_{cav}}$',color='orange')
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')

plt.axhline(y=2)
plt.grid()

fig3.set_figheight(9)
fig3.set_figwidth(13)
plt.show()#



""""

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


"""




