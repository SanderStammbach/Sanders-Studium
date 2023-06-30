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
import matplotlib.pyplot as plt
from IPython.display import display, Latex
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

omega_1=0
omega_2=30
omega_3=150

omega_f= omega_2 - omega_1
omega_h=omega_3-omega_1  # frequency of the atomic transition that is coupled to the hot bath
omega_c=omega_3-omega_2
print(ln(3))
omega_d=30

h=1
nph=30    # Maximale Photonen im cavity 
 
Th=100.    # temperature of the hot bath
Tc=20.     # temperature of the cold bath
Tenv=0.0000000000000000000000000001



nh=5
nc=0.02

nf=0.02    #Beschreibt den cavity/Photonen. 

f =0.03

gamma_h=1
gamma_c=1
kappa=0.028
kb=1
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

H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a

H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

V=f*a.dag()+f*a #das got glaub nid

H=H_free+H_int
Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)  

A1=Trans_13
A2=Trans_13.dag()
A3=Trans_23
A4=Trans_23.dag()
A5=a
A6=a.dag()


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

Delta1=Delta2=0
gamma_h = gamma_c = 1
g = 0.5
nc = ncav = 0.0
kappa = 0.2


Delta1=0
Delta2=0
anzahl=100
nh=0.5
nc=ncav=0
n_list=[]
nh_list=[]
f=0
Photonnumber_list=[]
nh2 = np.linspace(0, 70, 100)
for i in range(anzahl):
    n_list.append(np.abs(Diverse_Loups.EquationOfMotion2(Delta1 , Delta2 , f , nh, ncav , nc, gamma_c, gamma_h, g , kappa)))
    nh_list.append(nh)
    
    
    ist_temp=[]
    list_temp=Diverse_Loups.Photonnumber2(nh,a,proj_1,proj_2,proj_3,Trans_12,nc,ncav,gamma_h,gamma_c,kappa,A1,A2,A3,A4,A5,A6,omega_d,omega_f,omega_1,omega_2,f,g)
    #g_list.append(i/100)  #Erstellt eine Liste mit Wären von g 
    Photonnumber_list.append(list_temp)
    print(n_list[i],Photonnumber_list[i])
    nh=nh+0.7
    
    
    

fig4, ax = plt.subplots()
ax.set_xlabel(r' $n_h$', fontsize=21)
ax.set_ylabel(r' $\langle n \rangle$', fontsize=21)
plt.title(r' Photonnumber vs $n_h$',fontsize=21)
#plt.plot(np.asarray(nh_list2)[:100],np.asarray(Photonnumber_list)[:100],color='red',label='f=1')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(Photonnumber_list)[:anzahl],color='orange',label='numerical')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(n_list)[:anzahl],'--',color='black',label=r'analytical')
plt.plot(nh2, Diverse_Loups.N_Analytic2(gamma_h,kappa,g,nh2,ncav,nc),color='red',label='analytical')
legend = ax.legend(loc='center right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('C0')
plt.show()    
print(n_list,nh)
