# /bin/env/python
from calendar import c
from ctypes import c_char_p
from email import message_from_file
import imp
from tkinter import messagebox
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
from qutip import *
from qutip import ptrace 
from Loup_for_different_coupling import Energie as Energie
import csv
#Konstante Grössen
########################################################################################################
omega_1=0
omega_2=30
omega_3=150

omega_f= omega_2 - omega_1
omega_h=omega_3-omega_1  # frequency of the atomic transition that is coupled to the hot bath
omega_c=omega_3-omega_2

h=1
nph=60    # Maximale Photonen im cavity 
 
Th=1000.    # temperature of the hot bath
Tc=20.     # temperature of the cold bath
Tenv=1. 
g=1


nh=10
nc=0
nf=3


kappa=0.1
gamma_h=1
gamma_c=1
gamma_f=1
kb=1

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

a=qutip.tensor(qutip.identity(3),qutip.destroy(nph))


################################################################
#implementierung von dem Hamilton
H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a

H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

H=H_free+H_int

print(H-H.dag(),H_int-H_int.dag(),H_free-H_free.dag(),"sollte null geben!!!!!!!!!!!!!!!!!!!!!!!!")

#print(H_int,H_free,H)
#H=Hfree+Hint
#########################################################################################################



def n(omega,T):
    n=1/(np.exp(h*omega/(kb*T))-1)
    return n

#gamma_1=(n(omega_h,Th)+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
#gamma_2=(n(omega_h,Th))*gamma_h
#gamma_3=(n(omega_c,Tc)+1)*gamma_c
#gamma_4=(n(omega_c,Tc))*gamma_c
#gamma_5=(n(omega_f,Tenv)+1)*gamma_f ####goes to zero
#gamma_6=(n(omega_f,Tenv))*gamma_f

gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
gamma_2=(nh)*gamma_h
gamma_3=(nc+1)*gamma_c
gamma_4=(nc)*gamma_c
gamma_5=(nf+1)*gamma_f ####goes to zero
gamma_6=(nf)*gamma_f
print(gamma_1)


######################################################################################################
#Vorfaktoren rechenr
def T(omega,n):
    T=h*omega/(kb*(np.log((1/n))+1))
    return T

print("Die Temperatur des warmen Bades ist: ",T(omega_h,Th))
######################################################################################################
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
########################################################################################################
c_op_list=[]

c_op_list.append(np.sqrt(gamma_1)*A1)
c_op_list.append(np.sqrt(gamma_2)*A2)
c_op_list.append(np.sqrt(gamma_3)*A3)
c_op_list.append(np.sqrt(gamma_4)*A4)
c_op_list.append(np.sqrt(gamma_5)*A5)
c_op_list.append(np.sqrt(gamma_6)*A6)


#print(c_op_list)

rho = steadystate(H, c_op_list)
#print(rho)
#qutip.plot_wigner_fock_distribution(rho)
#plt.show()



rho_f=rho.ptrace(1)  ### State in the cavity
print(rho_f)
qutip.plot_wigner_fock_distribution(rho_f)
plt.show()
##########################################################################################################
#Kapitel Temperatur
#c_op_list.remove(np.sqrt(gamma_6)*A6)
#c_op_list.remove(np.sqrt(gamma_5)*A5)
#P=[]
#for i in range(nph):
#   P[i]=np.trace(b_fock[i]*b_fock[i].dag()*rho_f)
#print("erstes element", c_op_list[1])



##########################################################################################################
#Berechnen der Wärme als Tr(H*rho)
def D(c_op_list,rho):
    D=[]
    for i in range(6):
        D.append(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho-rho*c_op_list[i].dag()*c_op_list[i]))
    return D


    
Liste_von_Q=[]

Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1])))
Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3])))
Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5])))
print(Liste_von_Q)


float_list= list(np.float_(Liste_von_Q))
print(float_list)    
Liste_von_Q=float_list


g_list=[]

Energie_VS_g=[]
for i in range(100):
    list_temp=[]
    list_temp=Energie.EnergieCalculator(i/10, H_free,Trans_12,Trans_13, Trans_23,a,nh,nf,nc,h,kb,kappa,gamma_h,gamma_c,gamma_f,c_op_list)
    g_list.append(i/10)  #Erstellt eine Liste mit Wären von g 
    Energie_VS_g.append(list_temp)

#Liste von Stings in floats konvertieren
float_list2=list(np.float_(Energie_VS_g))
print(Energie_VS_g)  

#Speicherern der Liste in csv datei
with open('Speicherort.csv','w') as temp_file:
    for item in Energie_VS_g:
        temp_file.write("%s\n" % item)

################################################################################################################################################
#Plotten Wärme VS g 

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.set_xlabel('Kopplungskonstante g')
ax.set_ylabel('Wärme_Energie')
plt.plot(g_list,Energie_VS_g,"ob")
plt.show()

#with open("Speicherort.csv", "wb") as f:
#    writer = csv.writer(f)
#    writer.writerows(Energie_VS_g)

######################################################################################################################################################################
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
#print(ket*ket.dag())






#result=mesolve(H, rho0, tlist)
#print(D(c_op_list,rho)[3])


