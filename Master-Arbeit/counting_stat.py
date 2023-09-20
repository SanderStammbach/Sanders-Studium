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
from qutip import *
from qutip import ptrace 
import csv
from Loup_for_different_coupling import Diverse_Loups as Diverse_Loups
import multiprocessing as mp
import csv
#Konstante Grössen
########################################################################################################
omega_1=0
omega_2=30
omega_3=150

omega_f= omega_2 - omega_1
omega_h=omega_3-omega_1  # frequency of the atomic transition that is coupled to the hot bath
omega_c=omega_3-omega_2

omega_d=30

h=1
nph=30    # Maximale Photonen im cavity 
 
Th=100.    # temperature of the hot bath
Tc=20.     # temperature of the cold bath
Tenv=0.0000000000000000000000000001



nh=5
nc=0.02

nf=0.02    #Beschreibt den cavity/Photonen. 

f =0.8
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
""""
H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a

H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

V=f*a.dag()+f*a #das got glaub nid

H=H_free+H_int


#print(a, a.dag(), a*a*a.dag()-a*a.dag()*a)

#Hdilde=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 

Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)   

########################################################################################################
A1=Trans_13
A2=Trans_13.dag()
A3=Trans_23
A4=Trans_23.dag()
A5=a
A6=a.dag()

global gamma_1, gamma_2, gamma_3, gamma_4, kappa_5, kappa_6
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

#Fullcounting
A7=Trans_12

########################################################################################################
global c_op_list
c_op_list=[]

c_op_list.append(np.sqrt(gamma_1)*A1)
c_op_list.append(np.sqrt(gamma_2)*A2)
c_op_list.append(np.sqrt(gamma_3)*A3)
c_op_list.append(np.sqrt(gamma_4)*A4)
c_op_list.append(np.sqrt(kappa_5)*A5)
c_op_list.append(np.sqrt(kappa_6)*A6)
"""
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





vk_list=[np.sqrt(1/2*gamma_h)]
#print("J==========",qutip.to_super(c_op_list[0]))
def J_sup(nh, nc, nf,vk_list):
    J=[]
    c_op_list=colaps(nh, nc, nf)
    J=vk_list[0]*(qutip.to_super(c_op_list[0]))
    

    return J


def K_trace(nh, nc, nf,vk_list):
    K=[]
    rhoV=qutip.operator_to_vector(DichteMatrix(nh, nc, nf, Hdilde))
    c_op_list=colaps(nh, nc, nf)
    K=(vk_list[0]**2)*(np.trace(qutip.to_super(c_op_list[0])*rhoV)) 
    
    return K



print (J_sup(nh, nc, nf,vk_list),"ktrace==",K_trace(nh, nc, nf,vk_list))



rho=DichteMatrix(nh, nc, nf, Hdilde)
L=qutip.liouvillian(Hdilde,rho)
with open('studentsq.csv', 'w', newline='') as file:
     writer = csv.writer(file)
     writer.writerows(L)




def solveZ():

    Z = np.linalg.solve(L,rhoV) #uf der rechte site muess e vektor  stoh

    return Z