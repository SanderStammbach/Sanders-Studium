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

qutip.settings.has_mkl = False
omega_1=0
omega_2=30
omega_3=150

omega_f= omega_2 - omega_1
omega_h=omega_3-omega_1  # frequency of the atomic transition that is coupled to the hot bath
omega_c=omega_3-omega_2

omega_d=30

h=1
nph=15 # Maximale Photonen im cavity for fsc it works fine with 19
 
Th=100.    # temperature of the hot bath
Tc=20.     # temperature of the cold bath
Tenv=0.0000000000000000000000000001



nh=0.0013
nc=0.027

nf=0.5  #Beschreibt den cavity/Photonen. 

f =0.0
global kappa

#gamma_h=1
#gamma_c=20
#kappa=0.2
#kappa=0.028
kb=1
global g
g=2.8
#g=14*kappa
#g=f=0
#kappa=0.07
kappa=0.002
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



#print(qutip.to_super(np.sqrt(kappa_6)*A5))
#L=qutip.to_super(np.sqrt(kappa_6)*A5)


def Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d):
    H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a

    H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

    V=f*a.dag()+f*a #das got glaub nid

    H=H_free+H_int 

    Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)

    return Hdilde

def DichteMatrix(nh, nc, nf, Hami,kappa,gamma_h,gamma_c):
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


def colaps(nh, nc, nf,kappa,gamma_h,gamma_c):
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

vk_list=[-1,1]
#print("J==========",qutip.to_super(c_op_list[0]))
def J_sup(nh, nc, nf,vk_list,gamma_h,gamma_c): #J = the Drazin inverse times  alpha
    J=[]
    c_op_list=colaps(nh, nc, nf,kappa,gamma_h,gamma_c)
    J=vk_list[0]*(qutip.to_super(c_op_list[A]))+vk_list[1]*(qutip.to_super(c_op_list[B]))
    #J = [qutip.spre(c_op_list[0]) + qutip.spost(c_op_list[0].dag())]
    #J=np.array(J)
    return J


def K_Trace2(nh, nc, nf,vk_list,Hdilde,gamma_h,gamma_c): #K = the Drazin inverse times  alpha
    K=[]
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa,gamma_h,gamma_c)
    c_op_list=colaps(nh, nc, nf,kappa,gamma_h,gamma_c)
    K=np.trace((vk_list[0]**2)*((c_op_list[A].dag()*c_op_list[A])*rho)+(vk_list[1]**2)*((c_op_list[B].dag()*c_op_list[B])*rho))
    
    return np.real(K)



def solveZ( nh,nc,nf, vk_list,Hdilde,gamma_h,gamma_c): #Z = the Drazin inverse times  alpha
    c_op_list=colaps(nh, nc, nf,kappa,gamma_h,gamma_c)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa,gamma_h,gamma_c)
    L=qutip.liouvillian(Hdilde,c_op_list)
    rhoV=qutip.operator_to_vector(rho)
    global IdV
    IdV=(qutip.operator_to_vector(tensor(qutip.identity(3),qutip.identity(nph))))#bereits transponiert


    alpha=J_sup(nh, nc, nf,vk_list,gamma_h,gamma_c)*rhoV
        
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






def EnergieCalculator_mit_faktor(nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3):
   
   
    H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a
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





def Entropy_ohne_omega(nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3):
        
        H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a
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




def Dcalc(nh,nc,nf, vk_list,Hdilde,gamma_h,gamma_c):
    IdV=(qutip.operator_to_vector(tensor(qutip.identity(3),qutip.identity(nph))))
    D= K_Trace2(nh, nc, nf,vk_list,Hdilde,gamma_h,gamma_c)-2*np.real(IdV.dag()*J_sup(nh,nc,nf, vk_list,gamma_h,gamma_c)*solveZ(nh,nc,nf, vk_list,Hdilde,gamma_h,gamma_c))
    return(D)


anzahl=50
step=0.5
D_list=[]
nh_list=[]
Jh_list1=[]
x_line=[]
Q_list1=[]
Entropy_list1=[]
gamma_h=0.00001
Energy_list1=[]
for i in range(anzahl):
    gamma_c=51-gamma_h
    Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa,gamma_h,gamma_c)
    D=Dcalc(nh,nc,nf, vk_list,Hdilde,gamma_h,gamma_c)
    D_list.append(D[0])
    
    Lstrich=J_sup(nh, nc, nf,vk_list,gamma_h,gamma_c)
    Jh_list1.append(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])#noesomega h  dezue
    print(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])
    Entropy_list1.append(Entropy_ohne_omega(nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3)[3])
   
    Energy_list1.append(EnergieCalculator_mit_faktor(nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3))
    Q_list1.append((Entropy_list1[i])*(D_list[i]/(Energy_list1[i][0]**2)))
    x_line.append(gamma_c/gamma_h)
    gamma_h+=step
    print(i-anzahl, D_list[i])


D_list2=[]

Jh_list2=[]
x_line2=[]
Q_list2=[]
Entropy_list2=[]
gamma_h=0.00001
Energy_list2=[]
for i in range(anzahl):
    gamma_c=31/gamma_h
    Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
    rho=DichteMatrix(nh, nc, nf, Hdilde,kappa,gamma_h,gamma_c)
    D=Dcalc(nh,nc,nf, vk_list,Hdilde,gamma_h,gamma_c)
    D_list2.append(D[0])
    
    Lstrich=J_sup(nh, nc, nf,vk_list,gamma_h,gamma_c)
    Jh_list2.append(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])#noesomega h  dezue
    print(np.real((IdV.trans()*Lstrich*qutip.operator_to_vector(rho)))[0])
    Entropy_list2.append(Entropy_ohne_omega(nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3)[3])
   
   
    Energy_list2.append(EnergieCalculator_mit_faktor(nh,Trans_12,a, kb,h,g,proj_3,proj_1,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f,omega_3))
    Q_list2.append((Entropy_list2[i])*(D_list[i]/(Energy_list2[i][0]**2)))
    x_line2.append(gamma_c/gamma_h)
    gamma_c+=step
    print(i-anzahl, D_list[i])


label1=r' $Var \langle J_{h} \rangle $'
label12=r' $\langle J_{h} \rangle $'
label3=r' $\mathcal{Q}'
label4=r' $\dot{\sigma}-\frac{J_{cav}}{T_{cav}}$'
fig, (ax1,ax2) = plt.subplots(2)

ax1.set_xlabel(r' $\gamma_c/\gamma_h$', fontsize=19)
ax1.set_ylabel('D')
#plt.title(r' D vs $n_h$ ')
ax1.plot(np.asarray(x_line)[:anzahl],np.asarray(D_list2)[:anzahl],label=label1,color='black')
ax1.plot(np.asarray(x_line)[:anzahl],np.asarray(Jh_list2)[:anzahl],label=label12,color='red')
ax1.plot(np.asarray(x_line)[:anzahl],np.asarray(Q_list2)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
ax1.plot(np.asarray(x_line)[:anzahl],np.asarray(Entropy_list2)[:anzahl],'-',label=label4,color='orange')
legend = ax1.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')
ax1.axhline(y=2)
ax1.grid()

ax2.set_xlabel(r' $\gamma_c/\gamma_h$', fontsize=19)
ax2.set_ylabel('D')
#plt.title(r' D vs $n_h$ ')
ax2.plot(np.asarray(x_line)[:anzahl],np.asarray(D_list)[:anzahl],label=label1,color='black')
ax2.plot(np.asarray(x_line)[:anzahl],np.asarray(Jh_list1)[:anzahl],label=label12,color='red')
ax2.plot(np.asarray(x_line)[:anzahl],np.asarray(Q_list1)[:anzahl],'-',label=r' $\mathcal{Q}$',color='blue')
ax2.plot(np.asarray(x_line)[:anzahl],np.asarray(Entropy_list1)[:anzahl],'-',label=label4,color='orange')
legend = ax2.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')
ax2.axhline(y=2)
ax2.grid()

plt.show()