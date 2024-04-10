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
import cupy as cp
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams.update({
    'text.usetex': True,
})


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

nph=17 # Maximale Photonen im cavity for fsc it works fine with 19
 
Th=100.    # temperature of the hot bath
Tc=20.     # temperature of the cold bath
Tenv=0.0000000000000000000000000001

tz=13
nH=5
epsilon=0.15
nh_fix=nh=0.0013


nc_fix=nc=0.028

nf_fix=nf=0.05

#nf=nf_fix=0.0401 #Beschreibt den cavity/Photonen. 


anzahl=50
#anzahl=5
gamma_h=0.1
gamma_c=2

kb=1

g_fix=g=2.8

#kappa=0.07
kappa_fix=kappa=0.002
f_fix=f =0

#g_fix=g=14*kappa
#parameters without cavity entropy
"""
g_fix=g=0.15

#kappa=0.07
kappa_fix=kappa=2
f_fix=f =2
nh_fix=nh=5
"""



"""
Matrix = [[1, 2, 3], [0, 0, 4], [4, 0, 1]]
d=cp.asarray(Matrix)


e=cp.asnumpy(d)
print(e)

"""



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

S=qutip.squeezing(a,a,-1j*1)
#print(qutip.to_super(np.sqrt(kappa_6)*A5))
#L=qutip.to_super(np.sqrt(kappa_6)*A5)


def Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d):
    H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a

    H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

    V=f*a.dag()+f*a #das got glaub nid

    H=H_free+H_int 

    Hdilde=H_int+ S*V*S.dag() +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)
    
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




#print("J===",J_sup(nh, nc, nf,vk_list))
def K_trace(nh, nc, nf,vk_list,Hdilde,collaps_list): 
    K=[]
    rhoV=qutip.operator_to_vector(DichteMatrix(nh, nc, nf, Hdilde,kappa))
    c_op_list=collaps_list
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
    
   
    """""
    LForce_Matrix=np.matrix(LForce)
    RechtsForce_Matrix=np.matrix(RechtsForce)
    LForce_GPU=cp.asarray(LForce_Matrix)
    RechtsForce_GPU=cp.asarray(RechtsForce_Matrix)
    
    Zs_GPU=cp.linalg.lstsq(LForce_GPU,RechtsForce_GPU,rcond=None)
   
    Zs=cp.ndarray.get(Zs_GPU[0],stream=None,order='A',out=None)
   
    Z=qutip.Qobj(Zs,dims=[[[3, nph], [3, nph]], [1]])"""
    
    Zs=np.linalg.lstsq(LForce,RechtsForce)
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










def Dcalc(nh,nc,nf, vk_list,Hdilde):
    collaps_list=colaps(nh, nc, nf,kappa)
    IdV=(qutip.operator_to_vector(tensor(qutip.identity(3),qutip.identity(nph))))
    print(K_trace(nh,nc,nf,vk_list,Hdilde,collaps_list))
    D= K_Trace2(nh, nc, nf,vk_list,Hdilde,collaps_list)-2*np.real(IdV.dag()*J_sup(nh,nc,nf, vk_list,collaps_list)*solveZ(nh,nc,nf, vk_list,Hdilde))
    return(D)


def Potts(nh,nc,gammac,gammah,epsilon):
    Q=(((1 + nc)*nh + nc*(1 + nh))*np.log(((1 + nc)*nh)/(nc*(1 + nh))))/(-nc + nh) - (2*epsilon**2*gammac*gammah*(-nc + nh)*(gammac*nc + gammah*nh)*(-2/(gammac*nc + gammah*nh) + (gammac*gammah*(gammac*nc + gammah*nh)*(nc + nh + 3*nc*nh) + (gammac + gammah + 2*(gammac*nc + gammah*nh))*(4*epsilon**2 + (gammac*nc + gammah*nh)**2/4.))/ ((gammac*gammah*(gammac*nc + gammah*nh)**2*(nc + nh + 3*nc*nh))/4. + 2*epsilon**2*(gammac*nc + gammah*nh)*(gammac + gammah + (3*(gammac*nc + gammah*nh))/2.)))*np.log(((1 + nc)*nh)/(nc*(1 + nh))))/ ((gammac*gammah*(gammac*nc + gammah*nh)**2*(nc + nh + 3*nc*nh))/4. + 2*epsilon**2*(gammac*nc + gammah*nh)*(gammac + gammah + (3*(gammac*nc + gammah*nh))/2.))
    return Q


HdildeTest = Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,0,g,omega_d)
rhoTest=DichteMatrix(nh,nc,nf,HdildeTest,kappa)
if rhoTest.ptrace(1)[[nph-1], [nph-1]]>10**(-4):
    print("Error",rhoTest.ptrace(1)[[nph-1], [nph-1]])
    
        


step =0.08/anzahl
#step=20/anzahl
#step=0.2/anzahl
nc_Line=np.linspace(0, 0.6, 100, endpoint=True)
Q=Potts(nh,nc_Line,gamma_c,gamma_h,epsilon)
epsilon_line=np.linspace(0, 0.2, 100, endpoint=True)
Qf=Potts(nh,nc_Line,gamma_c,gamma_h,epsilon_line*g/kappa)



nc=0.0004
D_list=[]
nh_list=[]
Jh_list1=[]
Entropy_list1=[]
Q_list1=[]
Energy_list1=[]
Q_List_PRE=[]
#step=0.3
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
    Q_list1.append(np.float128(Entropy_list1[i])*(np.float128(D_list[i])/(np.float128(Jh_list1[i])*np.float128(Jh_list1[i]))))
    Q_List_PRE.append(Potts(nH,nc,gamma_c,gamma_h,epsilon))
    print(Q_List_PRE[i])
    #Q_list1.append((Entropy_list1[i])*(D_list[i]/(Jh_list1[i]**2)))
    nh_list.append(nc)
    nc+=step
    print(i-anzahl, D_list[i])








print(D_list)
print(Entropy_list1)


fig1, ax = plt.subplots()
ax.set_xlabel(r' $\frac{nc}{\gamma}$', fontsize=21)
#ax.set_ylabel(r' $\langle J \rangle $', fontsize=21)
plt.title('')


plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(Q_list1)[:anzahl],'-',color='blue',linewidth=3)
ax.axhline(y=2)

plt.show()
    

    
nc=nc_fix
g=g_fix   
f=f_fix 
kappa=kappa_fix
nh=nh_fix
omega_d=30   
D_list_nf=[]
nf_list=[]
Jh_list_nf=[]
Entropy_list_nf=[]
Q_list_nf=[]
Energy_list_nf=[]
step=0.2/anzahl
nf=0.0001
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
    print(i-anzahl,"parameter=",nh,nc,nf,g,kappa,f)


nc=nc_fix
nh=nh_fix
nf=nf_fix
f_list=[]
D_list_f=[]
Jh_list_f=[]
Entropy_list_f=[]
Q_list_f=[]
nc=nc_fix
g=g_fix
nh=nh_fix
f=0.0
step=0.3/anzahl
Potts_f=[]
Energy_list_f=[]

for i in range(anzahl):
    nf=0.09
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
    epsi=(f*g)/kappa
    Potts_f.append(Potts(nh,nc,gamma_c,gamma_h,epsi))
    f_list.append(f)
    f+=step
    print(i-anzahl)
    
    
    
    
nc=nc_fix
g=g_fix
nh=nh_fix
f=0.001
omega_f=30
omega_d=29.2
Delta_list=[]
step=1.6/anzahl

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
    Q_list_D.append((Entropy_list_D[i])*(D_list_D[i]/(Jh_list_D[i][0]**2)))

    Delta_list.append(omega_d-omega_f)
    omega_d+=step
    print(i-anzahl)


label1=r' $Var  \left( \frac{ J_{h} }{(~ \gamma_h \omega_{h})^2} \right)  $'
label12=r' $\frac{\langle J_{h} \rangle }{~\gamma_h \omega_{h}} $'
label3=r' $\mathcal{Q}'
label4=r' $\dot{\sigma}$'
fig, axs = plt.subplots(2,3)
axs[0,0].tick_params( labelsize=tz)
axs[0,1].tick_params( labelsize=tz)
axs[1,1].tick_params( labelsize=tz)
axs[1,0].tick_params( labelsize=tz)

axs[0,2].tick_params( labelsize=tz)
axs[1,2].tick_params( labelsize=tz)

axs[0,0].set_title(r' a)', fontsize=23)
axs[0,1].set_title(r' b)', fontsize=23)
axs[0,2].set_title(r' c)', fontsize=23)

axs[1,0].set_xlabel(r' $n_c$', fontsize=23)
#axs[0,0].set_ylabel(r'$\mathcal{Q}$', fontsize=23)
#axs[1,0].set_ylabel(r'$\mathcal{Q}$', fontsize=23)
#plt.title(r' D vs $n_h$ ')
#ax1.plot(np.asarray(nc_Line)[:anzahl],np.asarray(Q)[:anzahl],'*',label=r' $\mathcal{Q}$',color='blue',alpha=0.5)
axs[1,0].plot(np.asarray(nh_list)[:anzahl],np.asarray(D_list)[:anzahl],label=label1,color='black')
axs[1,0].plot(np.asarray(nh_list)[:anzahl],np.asarray(Jh_list1)[:anzahl],label=label12,color='red')

axs[0,0].plot(np.asarray(nh_list)[:anzahl],np.asarray(Q_list1)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')

#axs[0,0].plot(np.asarray(nh_list)[:anzahl],np.asarray(Q_List_PRE)[:anzahl],'-',label=r' $\mathcal{Q}  $',color='aqua')

axs[1,0].plot(np.asarray(nh_list)[:anzahl],np.asarray(Entropy_list1)[:anzahl],'-',label=label4,color='teal')
legend = axs[0,0].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend = axs[1,0].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend.get_frame().set_facecolor('white')

axs[0,0].axhline(y=2)

axs[0,0].grid()
axs[1,0].grid()



axs[1,1].set_xlabel(r' $n_f$', fontsize=23)
axs[0,1].set_ylabel(r' ')
#plt.title(r' $\mathcal{Q} vs f$ ')
axs[1,1].plot(np.asarray(nf_list)[:anzahl],np.asarray(D_list_nf)[:anzahl],label=label1,color='black')
axs[1,1].plot(np.asarray(nf_list)[:anzahl],np.asarray(Jh_list_nf)[:anzahl],label=label12,color='red')
axs[0,1].plot(np.asarray(nf_list)[:anzahl],np.asarray(Q_list_nf)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
axs[1,1].plot(np.asarray(nf_list)[:anzahl],np.asarray(Entropy_list_nf)[:anzahl],'-',label=r' $\dot{\sigma}$',color='teal')
legend = axs[0,1].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend = axs[1,1].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend.get_frame().set_facecolor('white')
axs[0,1].axhline(y=2)
axs[1,1].grid()
axs[0,1].grid()






axs[1,2].set_xlabel(r' $f/\gamma_h$', fontsize=23)
axs[0,2].set_ylabel(r' ')
#plt.title(r' $\mathcal{Q} vs f$ ')
axs[1,2].plot(np.asarray(f_list)[:anzahl],np.asarray(D_list_f)[:anzahl],label=label1,color='black')
axs[1,2].plot(np.asarray(f_list)[:anzahl],np.asarray(Jh_list_f)[:anzahl],label=label12,color='red')
axs[0,2].plot(np.asarray(f_list)[:anzahl],np.asarray(Q_list_f)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
#axs[0,2].plot(np.asarray(f_list),np.asarray(Potts_f),'-',label=r' $\mathcal{Q}$',color='aqua')
axs[1,2].plot(np.asarray(f_list)[:anzahl],np.asarray(Entropy_list_f)[:anzahl],'-',label=label4,color='teal')

axs[0,2].axhline(y=2)
axs[0,2].grid()
axs[1,2].grid()

legend = axs[0,2].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend = axs[1,2].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend.get_frame().set_facecolor('white')


fig.set_figheight(9)
fig.set_figwidth(13)
plt.savefig('plot.png', dpi=1500)

#plt.figure(dpi=300)
plt.show()#



fig, ax= plt.subplots()

plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(Energy_list1)[:anzahl,0],'-',color='red',label=r' $\frac{J_{h}}{~\gamma_h \omega_{h}}$')

plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(Jh_list1)[:anzahl],'*',label=r' $\langle \langle 1|\mathcal{J}|\rho \rangle \rangle $',color='red')

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

g=g_fix
step=1/anzahl
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
    print("parameter",nh,nc,f,g,kappa)
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


fig, ax= plt.subplots()

plt.plot(np.asarray(kappa_list)[:anzahl],np.asarray(Energy_list_kappa)[:anzahl,0],'-',color='red',label=r' $\frac{J_{h}}{~\gamma_h \omega_{h}}$')

plt.plot(np.asarray(kappa_list)[:anzahl],np.asarray(Jh_list_kappa)[:anzahl],'*',label=r' $\langle \langle 1|\mathcal{J}|\rho \rangle \rangle $',color='red')

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
    
    
    




fig2, axs2 = plt.subplots(2,3)

axs2[1,0].set_xlabel(r' $\kappa$', fontsize=23)
#axs2[0,0].set_ylabel(r' $\mathcal{Q}$')
#plt.title(r' D vs $n_h$ ')
axs2[0,0].tick_params( labelsize=tz)
axs2[0,1].tick_params( labelsize=tz)
axs2[1,1].tick_params( labelsize=tz)
axs2[1,0].tick_params( labelsize=tz)

axs2[0,2].tick_params( labelsize=tz)
axs2[1,2].tick_params( labelsize=tz)

axs2[1,0].plot(np.asarray(kappa_list)[:anzahl],np.asarray(D_list_kappa)[:anzahl],label=label1,color='black')
axs2[1,0].plot(np.asarray(kappa_list)[:anzahl],np.asarray(Jh_list_kappa)[:anzahl],label=label12,color='red')
axs2[0,0].plot(np.asarray(kappa_list)[:anzahl],np.asarray(Q_list_kappa)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
axs2[1,0].plot(np.asarray(kappa_list)[:anzahl],np.asarray(Entropy_list_kappa)[:anzahl],'-',label=r' $\dot{\sigma}$',color='orange')
legend = axs2[0,1].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend = axs2[1,1].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend.get_frame().set_facecolor('white')
axs2[0,0].axhline(y=2)
axs2[1,0].grid()
axs2[0,0].grid()



axs2[1,1].set_xlabel(r' $\frac{\Delta_1}{\gamma_h}$', fontsize=23)
#axs2[1,1].set_ylabel(r' $\mathcal{Q}$')
#plt.title(r' $\mathcal{Q} vs f$ ')
axs2[1,1].plot(np.asarray(Delta_list)[:anzahl],np.asarray(D_list_D)[:anzahl],label=label1,color='black')
axs2[1,1].plot(np.asarray(Delta_list)[:anzahl],np.asarray(Jh_list_D)[:anzahl],label=label12,color='red')
axs2[0,1].plot(np.asarray(Delta_list)[:anzahl],np.asarray(Q_list_D)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
axs2[1,1].plot(np.asarray(Delta_list)[:anzahl],np.asarray(Entropy_list_D)[:anzahl],'-',label=label4,color='teal')
legend = axs2[0,1].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend = axs2[1,1].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend.get_frame().set_facecolor('white')
#ax2.axhline(y=2)
axs2[0,1].grid()
axs2[0,1].axhline(y=2)
axs2[1,1].grid()


legend = plt.legend(loc='upper right', shadow=True, fontsize='xx-large')
legend.get_frame().set_facecolor('white')


axs2[1,2].set_xlabel(r' $g$', fontsize=23)
#axs2[1,2].set_ylabel(r' $\mathcal{Q}$')
#plt.title(r' $\mathcal{Q} vs f$ ')
axs2[1,2].plot(np.asarray(g_list)[:anzahl],np.asarray(D_list_g)[:anzahl],label=label1,color='black')
axs2[1,2].plot(np.asarray(g_list)[:anzahl],np.asarray(Jh_list_g)[:anzahl],label=label12,color='red')
axs2[0,2].plot(np.asarray(g_list)[:anzahl],np.asarray(Q_list_g)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
axs2[1,2].plot(np.asarray(g_list)[:anzahl],np.asarray(Entropy_list_g)[:anzahl],'-',label=label4,color='orange')

axs2[0,2].axhline(y=2)
axs2[1,2].grid()
axs2[0,2].grid()
legend = axs2[0,2].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend = axs2[1,2].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend.get_frame().set_facecolor('white')


fig2.set_figheight(9)
fig2.set_figwidth(13)
plt.savefig('plot.png', dpi=1500)
plt.show()


nf=nf_fix



















nc=nc_fix
nh=nh_fix
nf=nf_fix
f_list1=[]
D_list_f1=[]
Jh_list_f1=[]
Entropy_list_f1=[]
Q_list_f1=[]
nc=nc_fix
g=g_fix
nh=nh_fix
f=0.0
step=3.5/anzahl
Potts_f2=[]
Energy_list_f1=[]

for i in range(anzahl):
    g_fix=g=0.15


    kappa_fix=kappa=2
    
    nh_fix=nh=5
    nf=0.09
    Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa)
    D=Dcalc(nh,nc,nf, vk_list,Hdilde)
    D_list_f1.append(D[0])
    collapse_list=colaps(nh, nc, nf,kappa)
    print("bla==rho",np.trace(D))
    Lstrich=J_sup(nh, nc, nf,vk_list,collapse_list)
    Jh_list_f1.append(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])#noesomega h  dezue
    print(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])
    Energy_list_f1.append(EnergieCalculator_mit_faktor(collapse_list,g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h))
    Entropy_list_f1.append(Entropy_ohne_omega(collapse_list,nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3)[4])
    Q_list_f1.append((Entropy_list_f1[i])*(D_list_f1[i]/(Jh_list_f1[i]**2)))

    
    f_list1.append(f)
    f+=step
    print(i-anzahl)









f_list2=[]
D_list_f2=[]
Jh_list_f2=[]
Entropy_list_f2=[]
Q_list_f2=[]


nh=nh_fix
f=0.0
step=3.5/anzahl
Potts_f=[]
Energy_list_f2=[]






for i in range(anzahl):
    g_fix=g=0.015


    kappa=0.2
   
    nh=5
    nf=0.09
    Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa)
    D=Dcalc(nh,nc,nf, vk_list,Hdilde)
    D_list_f2.append(D[0])
    collapse_list=colaps(nh, nc, nf,kappa)
    print("bla==rho",np.trace(D))
    Lstrich=J_sup(nh, nc, nf,vk_list,collapse_list)
    Jh_list_f2.append(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])#noesomega h  dezue
    print(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])
    Energy_list_f2.append(EnergieCalculator_mit_faktor(collapse_list,g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h))
    Entropy_list_f2.append(Entropy_ohne_omega(collapse_list,nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3)[3])
    Q_list_f2.append((Entropy_list_f2[i])*(D_list_f2[i]/(Jh_list_f2[i]**2)))
    epsi=(f*g)/kappa
    Potts_f2.append(Potts(nh,nc,gamma_c,gamma_h,epsi))
    f_list2.append(f)
    f+=step
    print(i-anzahl)














fig4, axs = plt.subplots(2,2)

axs[1,0].set_xlabel(r' $f/\gamma_h$', fontsize=23)
#axs2[0,0].set_ylabel(r' $\mathcal{Q}$')
#plt.title(r' D vs $n_h$ ')
axs[0,0].set_title(r' a)', fontsize=23)
axs[0,1].set_title(r' b)', fontsize=23)
axs[0,0].tick_params( labelsize=tz)
axs[0,1].tick_params( labelsize=tz)
axs[1,1].tick_params( labelsize=tz)
axs[1,0].tick_params( labelsize=tz)

axs[1,0].set_xlabel(r' $f/\gamma_h$', fontsize=23)
axs[0,0].set_ylabel(r'', fontsize=23)
#plt.title(r' $\mathcal{Q} vs f$ ')
axs[1,0].plot(np.asarray(f_list1)[:anzahl],np.asarray(D_list_f1)[:anzahl],label=label1,color='black')
axs[1,0].plot(np.asarray(f_list1)[:anzahl],np.asarray(Jh_list_f1)[:anzahl],label=label12,color='red')
axs[0,0].plot(np.asarray(f_list1)[:anzahl],np.asarray(Q_list_f1)[:anzahl],'-',label=r' $\mathcal{Q\prime}$',color='blue')
axs[0,0].plot(np.asarray(f_list1),np.asarray(Potts_f2),'-',label=r' $\mathcal{Q\prime \prime}$',color='aqua')
axs[1,0].plot(np.asarray(f_list1)[:anzahl],np.asarray(Entropy_list_f1)[:anzahl],'-',label=label4,color='teal')

axs[0,0].axhline(y=2)
axs[0,0].grid()
axs[1,0].grid()

legend = axs[0,0].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend = axs[1,0].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend.get_frame().set_facecolor('white')







axs[1,1].set_xlabel(r' $f/\gamma_h$', fontsize=23)
axs[0,1].set_ylabel(r' ', fontsize=23)
#plt.title(r' $\mathcal{Q} vs f$ ')
axs[1,1].plot(np.asarray(f_list2)[:anzahl],np.asarray(D_list_f2)[:anzahl],label=label1,color='black')
axs[1,1].plot(np.asarray(f_list2)[:anzahl],np.asarray(Jh_list_f2)[:anzahl],label=label12,color='red')
axs[0,1].plot(np.asarray(f_list2)[:anzahl],np.asarray(Q_list_f2)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
axs[0,1].plot(np.asarray(f_list2),np.asarray(Potts_f2),'-',label=r' $\mathcal{Q\prime \prime}$',color='aqua')
axs[1,1].plot(np.asarray(f_list2)[:anzahl],np.asarray(Entropy_list_f2)[:anzahl],'-',label=label4,color='teal')

axs[0,1].axhline(y=2)
axs[0,1].grid()
axs[1,1].grid()

legend = axs[0,1].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend = axs[1,1].legend(loc='upper right', shadow=True, fontsize='xx-large')
legend.get_frame().set_facecolor('white')
fig4.set_figheight(9)
fig4.set_figwidth(13)
plt.show()









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
nh_fix=nh=0.0013
f=0
kappa=0.002
nc_fix=nc=0.028
g=2.8
nf_fix=nf=0.05
omega_f=30
anzahl=40
step =0.08/anzahl

nc=0.0004
omega_d=30
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
    Q_list.append((Entropy_list[i])*(d_list[i]/(Jh_list2[i][0]**2)))
    
    nh_list.append(nc)
    nc=nc+step
    







fig3, ax = plt.subplots()

ax.set_xlabel(r' $n_c$', fontsize=23)
ax.set_ylabel('D')
plt.title(r' D vs $n_c$ ')
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




