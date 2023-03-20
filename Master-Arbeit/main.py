# /bin/env/python
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
from qutip import steadystate as steadystate

import math
from qutip.qobj import *
from qutip.states import *
from qutip.operators import *
from qutip import *
from qutip import ptrace
from qutip import variance as variance
from Loup_for_different_coupling import Diverse_Loups as Diverse_Loups
import multiprocessing as mp
import csv
from numpy import log as ln
from IPython.display import display, Latex
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
#Konstante Grössen
########################################################################################################
omega_1=0
omega_2=30
omega_3=150

omega_f= omega_2 - omega_1
omega_h=omega_3-omega_1  # frequency of the atomic transition that is coupled to the hot bath
omega_c=omega_3-omega_2
print(ln(3))
omega_d=omega_1/(ln(2)-1)

h=1
nph=30    # Maximale Photonen im cavity 
 
Th=100.    # temperature of the hot bath
Tc=20.     # temperature of the cold bath
Tenv=0.0000000000000000000000000001



nh=2.6
nc=0.02

nf=0.02    #Beschreibt den cavity/Photonen. 

f =100

gamma_h=1
gamma_c=1
kappa=0.028
kb=1
g=14*kappa
 

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

H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

V=f*a.dag()+f*a #das got glaub nid

H=H_free+H_int




Hdilde=H_free+H_int -omega_d*(a.dag()+a*proj_2) + f*(a+a.dag()) 



########################################################################################################
A1=Trans_13
A2=Trans_13.dag()
A3=Trans_23
A4=Trans_23.dag()
A5=a
A6=a.dag()

"""
gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
gamma_2=(nh)*gamma_h
gamma_3=(nc+1)*gamma_c
gamma_4=(nc)*gamma_c
kappa_5=(nf+1)*kappa ####goes to zero
kappa_6=(nf)*kappa

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
"""
def DichteMatrix(nh, nc, nf, Hami):
    gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
    gamma_2=(nh)*gamma_h
    gamma_3=(nc+1)*gamma_c
    gamma_4=(nc)*gamma_c
    kappa_5=(nf+1)*kappa ####goes to zero
    kappa_6=(nf)*kappa

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

rho = DichteMatrix(nh,nc,nf,Hdilde)

#print(c_op_list)


#print(rho)
#qutip.plot_wigner_fock_distribution(rho)
#plt.show()



rho_f=rho.ptrace(1)  ### State in the cavity
print(rho)
qutip.plot_wigner_fock_distribution(rho_f,colorbar='colorbar')
plt.show()


#Hamilton
print(variance(a.dag()*a, rho))

print(type(rho))

"""Mk=qutip.superop_reps.to_kraus(Trans_12, tol=1e-09)
Mkt=qutip.superop_reps.to_kraus(Trans_12.dag(), tol=1e-09)
print(Mk)
print(Mkt)
rhok=Mk*rho*Mkt
print(rhok)
"""

L_List=[A1,A3,A5]
Lt_list=[A2,A4,A6]
Mk_List=[]
Mkt_List=[]
for i in range (2):
    Mk_List.append(qutip.superop_reps.to_kraus(L_List[i], tol=1e-09))
    Mkt_List.append(qutip.superop_reps.to_kraus(Lt_list[i], tol=1e-09))


print(np.vdot(Mk_List[1],Mkt_List[1])*rho)

def CalculateK(Mk_List,Mkt_List,rho):
    K=0.0
    for i in range (2):
        K=K+np.trace(np.vdot(Mk_List[i],Mkt_List[i])*rho) #vk² muess no  dezue
        
    return K

def CalculateJ(Mk_List,Mkt_List,rho):
    J=0.0
    for i in range (2):
        J=J+np.trace(np.vdot(Mk_List[i],Mkt_List[i])*rho) #2*vk muess no dezue
        
    return 2*J

print("mis k  sött " ,CalculateK(Mk_List,Mkt_List,rho))

print(np.trace(rho))

#D =j+k und das isch d zitableitig vo de varianz
#print(L_List)
