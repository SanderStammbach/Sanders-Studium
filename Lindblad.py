# /bin/env/python
from ctypes import c_char_p
from email import message_from_file
import imp
from tkinter import messagebox
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
from qutip import *
#Konstante Grössen
##############################################################################
omega_1=0
omega_2=30
omega_3=150

omega_f= omega_2 - omega_1
omega_h=omega_3-omega_1  # frequency of the atomic transition that is coupled to the hot bath
omega_c=omega_3-omega_2

h=1
nph=30     # Maximale Photonen im cavity 
 
Th=100.    # temperature of the hot bath
Tc=20.     # temperature of the cold bath
Tenv=0. 
g=5

b_fock=qutip.states.fock(nph,0) #m)/fock(N,#m)
b_atom=basis(3)
b_comp=tensor( b_atom, b_fock)



#psi1=basis(b_atom,1)
#psi2=basis(b_atom,2)
#psi3=basis(b_atom,3)

# hier ist ein wenig gebastel mit den transitionsoperatoren

va=qutip.Qobj(qutip.qutrit_basis()[2])
vb=qutip.Qobj(qutip.qutrit_basis()[1])
vg=qutip.Qobj(qutip.qutrit_basis()[0])

Trans_13=tensor(vg*va.dag(),qutip.identity(nph))
Trans_23=tensor(vg*vb.dag(),qutip.identity(nph))
Trans_12=tensor(va*vb.dag(),qutip.identity(nph))

proj_1=tensor(vg*vg.dag(),qutip.identity(nph))
proj_2=tensor(vb*vb.dag(),qutip.identity(nph))
proj_3=tensor(va*va.dag(),qutip.identity(nph))

a=qutip.tensor(qutip.destroy(3),qutip.identity(nph))

################################################################
#implementierung von dem Hamilton
H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a

H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

H=H_free+H_int

print(H_int,H_free,H)
#H=Hfree+Hint
################################################################
kappa=0.1
gamma_h=1
gamma_c=1
gamma_f=1
kb=1


def n(omega,T):
    n=1/(np.exp(h*omega/(kb*T))-1)
    return n

gamma_1=(n(gamma_h,Th)+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
gamma_2=(n(gamma_h,Th))*gamma_h
gamma_3=(n(gamma_c,Th)+1)*gamma_c
gamma_4=(n(gamma_c,Th))*gamma_c
gamma_5=(n(gamma_f,Th)+1)*gamma_f
gamma_6=(n(gamma_f,Th))*gamma_f

print(gamma_1)
############################################################################################
#A1=qutip.tensor(Trans_13,qutip.identity(nph))
#A2=qutip.tensor(Trans_13.dag(),qutip.identity(nph))
#A3=qutip.tensor(Trans_23,qutip.identity(nph))
#A4=qutip.tensor(Trans_23.dag(),qutip.identity(nph))
#A5=qutip.tensor(qutip.identity(nph),a) 
#A6=qutip.tensor(qutip.identity(nph),a.dag())

A1=Trans_13
A2=Trans_13.dag()
A3=Trans_23
A4=Trans_23.dag()
A5=a
A6=a.dag()
##################################################################################################
c_op_list=[]

c_op_list.append(np.sqrt(gamma_1)*A1)
c_op_list.append(np.sqrt(gamma_2)*A2)
c_op_list.append(np.sqrt(gamma_3)*A3)
c_op_list.append(np.sqrt(gamma_4)*A4)
c_op_list.append(np.sqrt(gamma_5)*A5)
c_op_list.append(np.sqrt(gamma_6)*A6)


print(c_op_list)

rho = steadystate(H, c_op_list)
print(rho)
qutip.plot_wigner_fock_distribution(rho)
plt.show()

#########################################################
#random testing

H = 2*np.pi * 0.1 * qutip.sigmax()
psi0 = basis(2, 0)
times = np.linspace(0.0, 10.0, 100)

result = qutip.sesolve(H, psi0, times, [qutip.sigmaz(), qutip.sigmay()])
fig, ax = plt.subplots()
ax.plot(result.times, result.expect[0])
ax.plot(result.times, result.expect[1])
ax.set_xlabel('Time')
ax.set_ylabel('Expectation values')
ax.legend(("Sigma-Z", "Sigma-Y"))
#plt.show()
ket = basis(5,2)
print(ket*ket.dag())
#result=mesolve(H, rho0, tlist)
