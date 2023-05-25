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
import matplotlib.pyplot as plt
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
#g=0

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


print(a, a.dag(), a*a*a.dag()-a*a.dag()*a)

#Hdilde=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 

Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)   
print("dfjk",omega_d)
########################################################################################################
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

print("D lösig isch:",np.trace((a*(c_op_list[4]*rho*c_op_list[4].dag()-1/2*(c_op_list[4].dag()*c_op_list[4]*rho-rho*c_op_list[4].dag()*c_op_list[4]))))-np.trace(a*a.dag()*a),"        und       ",a*(c_op_list[4]*rho*c_op_list[4].dag()-1/2*(c_op_list[4].dag()*c_op_list[4]*rho-rho*c_op_list[4].dag()*c_op_list[4])))

rho_f=rho.ptrace(1)  ### State in the cavity
print(rho)
qutip.plot_wigner_fock_distribution(rho_f,colorbar='colorbar')
plt.show()


#Hamilton
#print(variance(a.dag()*a, rho))

#print(type(rho))

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
#do weiss i noni gnau welli liste öb mk mit superop to kraus oder ni
def CalculateK(Mk_List,Mkt_List,rho):
    K=0.0
    k=1
    for i in range (2):
        if i==0:
            vk=1/(np.sqrt(2*k))
        else:
            vk=1/np.sqrt(k)

        K=K+np.trace(np.vdot(Mk_List[i],Mkt_List[i])*rho) #vk² muess no  dezue
        
    return K

def CalculateJ(Mk_List,Mkt_List,rho):
    J=0.0
    k=1
    for i in range (2):
        if i==0:
            vk=1/(np.sqrt(2*k))
        else:
            vk=1/np.sqrt(k)
            
        J=J+ 2*vk*np.trace(np.vdot(Mk_List[i],Mkt_List[i])*rho) #2*vk muess no dezue
        
    return 2*J

print("mis k  sött " ,CalculateK(Mk_List,Mkt_List,rho),"print j",CalculateJ(Mk_List,Mkt_List,rho))

print(np.trace(rho))

#D =j+k und das isch d zitableitig vo de varianz
#print(L_List)

########################################################################################################################################################################
g_list=[]

Energie_VS_g=[]
for i in range(200):
    #list_temp=[]
    #list_temp=Diverse_Loups.EnergieCalculator(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,omega_2,V,omega_1,omega_d,proj_2,omega_f,c_op_list)
    g_list.append(i/100)  #Erstellt eine Liste mit Wären von g 
    #Energie_VS_g.append(list_temp)

"""
fig, ax = plt.subplots()
ax.set_xlabel(r' $\frac{g}{\gamma_h}$', fontsize=23)
ax.set_ylabel(r' Heat current', fontsize=15)
plt.title('current/energy flux vs coupling constant')
plt.plot(np.asarray(g_list)[:200],np.asarray(Energie_VS_g)[:200,0],label=r' $\frac{J_h}{\gamma_h \omega_h}$')
plt.plot(np.asarray(g_list)[:200],np.asarray(Energie_VS_g)[:200,1],label=r' $\frac{J_c}{\gamma_c \omega_c}$')
plt.plot(np.asarray(g_list)[:200],np.asarray(Energie_VS_g)[:200,2],label=r' $\frac{J_{cav}}{\gamma_{cav} \omega_{cav}}$')
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('C0')
plt.show()
"""

P_list=Diverse_Loups.P(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,f,omega_2)
print("liste von P",P_list)
PowerPlot, ax = plt.subplots() 
ax.set_xlabel(r' g', fontsize=23)
ax.set_ylabel(r' Power', fontsize=15)
plt.title("Power vs  coupling-constant")
plt.plot(np.asarray(g_list)[:200],np.asarray(P_list)[:200],label=r' Kurve')


plt.show()

"""

g=0
g_li=[]
P_li=[]
for i in range(200):
    g=g+1/80
    H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

    H=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 
        
    V=f*a.dag()+f*a
    Hdilde=H_int+V +(30-(30+omega_d))*(a.dag()*a)+(omega_f-omega_d)*proj_2  
    rho = steadystate(Hdilde, c_op_list) ######## Are you sure its with only photons H_free?
        
    g_li.append(g)
    P_li.append(Diverse_Loups.P4(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,omega_2,g,f))

    print(P_li)





PowerPlot2, ax = plt.subplots() 
ax.set_xlabel(r' g', fontsize=23)
ax.set_ylabel(r' Power', fontsize=15)
plt.plot(np.asarray(g_li)[:200],np.asarray(P_li)[:200],label=r' Kurve')
plt.show()



"""




#Entropy. production
def T(omega,n):
    T=h*omega/(kb*(np.log((1/n)+1)))
    return T


"""
nh_list=[]
Trace_list=[]
nh=0.1 #set nh again to zero
for j in range(100):
    Trace_list_temp=Diverse_Loups.Funktion(nh,proj_1,proj_2,proj_3,H,nc,nf,gamma_h,gamma_c,kappa,A1,A2,A3,A4,A5,A6)
    Trace_list.append(Trace_list_temp)

    nh_list.append(nh)
    nh=nh+0.3
"""
nh2=0.1
nh_list2=[]
Entropy=[]
for i in range(100):
    list_temp=[]
    list_temp=Diverse_Loups.Entropy(nh2,Trans_12,a, kb,h,g,H,H_free,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_c,omega_h,omega_f)
    #g_list.append(i/100)  #Erstellt eine Liste mit Wären von g 
    Entropy.append(list_temp)
    nh2=nh2+0.3
    nh_list2.append(nh2)

#Liste von Stings in floats konvertieren
#float_list2=list(np.float_(Energie_VS_g))
print(Entropy) 

#result=mesolve(H, rho0, tlist)
#print(D(c_op_list,rho)[3])


print("Die Temperatur des warmen Bades ist: ",T(omega_h,nh))
print("Die Temperatur des kalten Bades ist: ",T(omega_c,nc))


fig3, ax = plt.subplots()

ax.set_xlabel(r' $n_h$', fontsize=19)
ax.set_ylabel('Entropy production rate')
plt.title(r' Entropy Production  rate vs $n_h$ ')
plt.plot(np.asarray(nh_list2)[:100],np.asarray(Entropy)[:100,0],label=r' $\frac{J_h}{T_h}+\frac{J_{cav}}{T_{cav}}+\frac{J_c}{T_c}$',color='red')
plt.plot(np.asarray(nh_list2)[:100],np.asarray(Entropy)[:100,1],label=r' $\frac{J_h}{T_h}$',color='green')
plt.plot(np.asarray(nh_list2)[:100],np.asarray(Entropy)[:100,2],label=r' $\frac{J_c}{T_c}$',color='pink')
plt.plot(np.asarray(nh_list2)[:100],np.asarray(Entropy)[:100,3],label=r' $\frac{J_{cav}}{T_{cav}}$',color='orange')
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('C0')
#Linien in plt
"""plt.axvline(x=2.6)
plt.axvline(x=2.6)
plt.axvline(x=5.5)
plt.axvline(x=0.17)
plt.axvline(x=20)
plt.axvline(x=1.7)"""

plt.show()

################################################################

g=14*kappa
f=0
f1=f
f_list=[]
for i in range(200):
    f1=f1+1/80
    f_list.append(f1)


anzahl=200 #anzahl iterationen im loop

#f against the  power
P_list=Diverse_Loups.P3(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,omega_2,g,f,anzahl)
print("liste von P",P_list)


PowerPlot, ax = plt.subplots() 
ax.set_xlabel(r'$f$', fontsize=23)
ax.set_ylabel(r' Power', fontsize=15)
plt.title('power vs driven field-strenght')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(P_list)[:anzahl],label=r' Kurve')
plt.show()


#f against power


Energie_VS_f=Diverse_Loups.current(H_free, Trans_12, a, h,c_op_list,omega_d,omega_f ,proj_2,g,f,anzahl)

PundJ=[]
Energie_VS_f2=np.array(Energie_VS_f)
P_list2=np.array(P_list)
for i in range(anzahl):
    PundJ.append(Energie_VS_f2[i,0]+Energie_VS_f2[i,1]+Energie_VS_f2[i,2]+P_list2[i])



fig, ax = plt.subplots()
ax.set_xlabel(r' $\frac{f}{\gamma_h}$', fontsize=23)
ax.set_ylabel(r' Heat current or power ', fontsize=15)
plt.title('current/power vs driven field')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(P_list)[:anzahl],'--',label=r'$ \frac{P}{\hbar \gamma_h \omega_{h}}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energie_VS_f)[:anzahl,0],label=r' $\frac{J_h}{\hbar \gamma_h \omega_h}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energie_VS_f)[:anzahl,1],label=r' $\frac{J_c}{\hbar\gamma_h \omega_h}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energie_VS_f)[:anzahl,2],label=r' $\frac{J_{cav}}{\hbar\gamma_h \omega_{h}}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(PundJ)[:anzahl],'*',label=r'$\hbar \frac{P+J_c+J_h+J_{cav}}{\hbar \gamma_h \omega_{h}}$')
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')
plt.show()


print(f_list)





###############################################################################################


f=0.3
g_list=[]
g=0
Energie_VS_g=[]
for i in range(200):
    g=g+1/120
    list_temp=[]
    list_temp=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2)
    g_list.append(g)  #Erstellt eine Liste mit Wären von g 
    Energie_VS_g.append(list_temp)

g=0
P_list=Diverse_Loups.P4(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,omega_2,g,f,anzahl)
print("liste von P",P_list)
PowerPlot, ax = plt.subplots() 
ax.set_xlabel(r' g', fontsize=23)
ax.set_ylabel(r' Power', fontsize=15)
plt.title("Power vs  coupling-constant")
plt.plot(np.asarray(g_list)[:200],np.asarray(P_list)[:200],label=r' Kurve')

fig, ax = plt.subplots()
ax.set_xlabel(r' $\frac{g}{\gamma_h}$', fontsize=23)
ax.set_ylabel(r' Heat current', fontsize=15)
plt.title('current/power vs g')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(P_list)[:anzahl],'--',label=r'$ \frac{P}{\hbar \gamma_h \omega_{h}}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energie_VS_f)[:anzahl,0],label=r' $\frac{J_h}{\hbar \gamma_h \omega_h}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energie_VS_f)[:anzahl,1],label=r' $\frac{J_c}{\hbar\gamma_h \omega_h}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energie_VS_f)[:anzahl,2],label=r' $\frac{J_{cav}}{\hbar\gamma_h \omega_{h}}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(PundJ)[:anzahl],'*',label=r'$\hbar \frac{P+J_c+J_h+J_{cav}}{\hbar \gamma_h \omega_{h}}$')
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')
plt.show()




plt.show()