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
from Loup_for_different_coupling import Diverse_Loups as Diverse_Loups
import multiprocessing as mp
import csv

import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import numpy as np
#Konstante Grössen
########################################################################################################
plt.rcParams.update({
    'text.usetex': True,
})
fz=25
FZ=29
tz=15
global omega_1, omega_2,  omega_d,omega_f
Schrrift=30
omega_1=0
omega_2=30
omega_3=150

omega_f= omega_2 - omega_1
omega_h=omega_3-omega_1  # frequency of the atomic transition that is coupled to the hot bath
omega_c=omega_3-omega_2

omega_d=30

h=1
nph=30  # Maximale Photonen im cavity 
 
Th=100.    # temperature of the hot bath
Tc=20.     # temperature of the cold bath
Tenv=0.0000000000000000000000000001



nh=5
nc=0.0002

nf=0.0002#Beschreibt den cavity/Photonen. 
f =0

global kappa

gamma_h=1
gamma_c=1
kappa=0.07
#kappa=0.2
#kappa=0.028
kb=1
global g
g=2.8
#g=14*kappa
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

"""
aa=operator_to_vector=tensor(Trans_13)
IdV=(qutip.operator_to_vector(tensor(qutip.identity(3),qutip.identity(nph))))

Id=np.matrix(IdV)
c= tensor(aa.trans(),aa)
print("c====",np.transpose(Id)*c*Id)
print(1/0)"""

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


def Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d):
    H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a

    H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

    V=f*a.dag()+f*a #das got glaub nid
    S=qutip.squeezing(a,a,-1j*1)
    
    H=H_free+H_int 

    Hdilde=H_int+V.dag() +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)

    return Hdilde


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

g=2.8
kappa=0.07
f=0.01
gammac=gammah=1
nh=5
#Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)

rho = DichteMatrix(nh,nc,nf,Hdilde)

nh_Line=np.linspace(0, 200, 70)
f_Line=np.linspace(0, 1, 70)
JH=Diverse_Loups.High_f_Approx(kappa,0.1,gamma_c,gamma_h,nc,nh_Line,g,nf)[0]
JC=omega_c/omega_h*Diverse_Loups.High_f_Approx(kappa,0.1,gamma_c,gamma_h,nc,nh_Line,g,nf)[1]
P=omega_f/omega_h*Diverse_Loups.High_f_Approx(kappa,0.1,gamma_c,gamma_h,nc,nh_Line,g,nf)[5]
JCAV=omega_f/omega_h*Diverse_Loups.High_f_Approx(kappa,0.1,gamma_c,gamma_h,nc,nh_Line,g,nf)[2]
TOT=JH+JC+JCAV+P
JH_f=Diverse_Loups.High_f_Approx(kappa,f_Line,gamma_c,gamma_h,nc,5,g,nf)[0]
JC_f=omega_c/omega_h*Diverse_Loups.High_f_Approx(kappa,f_Line,gamma_c,gamma_h,nc,5,g,nf)[1]
#JCAV=Diverse_Loups.High_f_Approx(kappa,f,gamma_c,gamma_h,nc,nh_Line,g,nf)[2]
#JH=(4*f**2*g**2*gammac*gammah*(-nc + nh_Line))/(gammac*gammah*kappa**2*(gammac*nc + gammah*nh_Line)*(nc + nh_Line + 3*nc*nh_Line) + 4*f**2*g**2*(gammac*(2 + 3*nc) + gammah*(2 + 3*nh_Line)))
fig3,ax1 = plt.subplots()
plt.title('')
ax1.set_xlabel(r' $n_h$', fontsize=FZ)
ax1.set_ylabel('current')
plt.rc ('xtick', labelsize = 15) #Schriftgröße der x-Tick-Labels
plt.rc ('ytick', labelsize = 15) #Schriftgröße der y-Tick-Labels


print(JH)

ax1.plot(nh_Line,JC,color='blue',label=r'$\frac{\langle J_{h}\rangle }{~\gamma_h \omega_{h}} $ ')
ax1.plot(nh_Line,JH,color='red',label=r'$\frac{\langle J_{c}\rangle }{~\gamma_h \omega_{h}} $ ')
ax1.plot(nh_Line,JCAV,color='orange',label=r'$\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}} $ ')
ax1.plot(nh_Line,P,color='green',label=r'$\frac{\langle P\rangle }{~\gamma_h \omega_{h}} $ ')
ax1.plot(nh_Line,TOT,color='black',label=r'$\frac{\langle Tot\rangle }{~\gamma_h \omega_{h}} $ ')


axins= ax1.inset_axes([0.5,0.5,0.5,0.5])
axins.plot(nh_Line,JH)

mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="0.5")

plt.draw()
plt.show()

#ax1.plot(f_Line,JC_f,color='red',label=r'$\langle J_h \rangle $ high $f$ approx')
#ax1.plot(f_Line,JH_f,color='red',label=r'$\langle J_c \rangle $ high $f$ approx')
legend = ax1.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
fig3.set_figheight(9)
fig3.set_figwidth(13)
#plt.plot(nh_Line,JCAV,color='pink',label=r'$\langle a\rangle$ "from paper"')
plt.show()


"""""
def get_demo_image():
    from matplotlib.cbook import get_sample_data
    import numpy as np
    f = get_sample_data("axes_grid/bivariate_normal.npy", asfileobj=False)
    z = np.load(f)
    # z is a numpy array of 15x15
    return z, (-3,4,-4,3)

fig, ax = plt.subplots(figsize=[5,4])

# prepare the demo image
Z, extent = get_demo_image()
Z2 = np.zeros([150, 150], dtype="d")
ny, nx = Z.shape
Z2[30:30+ny, 30:30+nx] = Z

# extent = [-3, 4, -4, 3]
ax.imshow(Z2, extent=extent, interpolation="nearest",
         origin="lower")

axins = zoomed_inset_axes(ax, 6, loc=1) # zoom = 6
axins.imshow(Z2, extent=extent, interpolation="nearest",
            origin="lower")

# sub region of the original image
x1, x2, y1, y2 = -1.5, -0.9, -2.5, -1.9
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

plt.xticks(visible=False)
plt.yticks(visible=False)

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

plt.draw()
plt.show()
"""

f=0.2
Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,0.3,0,omega_d)
rho1 = DichteMatrix(0.01,0.01,0.01,Hdilde)
rho_c=rho1.ptrace(1)
qutip.plot_wigner_fock_distribution(rho_c,colorbar='colorbar')



xvec = np.linspace(-5,5,200)

f=0
Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
rho = DichteMatrix(nh,nc,nf,Hdilde)
rho_f=rho.ptrace(1)
W_coherent = qutip.wigner(rho_f, xvec, xvec)

f=0.1
Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
rho = DichteMatrix(nh,nc,nf,Hdilde)
rho_f=rho.ptrace(1)
W_thermal = qutip.wigner(rho_f, xvec, xvec)


f=0.3
Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
rho = DichteMatrix(nh,nc,nf,Hdilde)
rho_f=rho.ptrace(1)
W_fock = qutip.wigner(rho_f, np.linspace(-10,10,200), np.linspace(-10,10,200))
#print("wignerplot an stelle 9.9 ",W_coherent[2][0.0])
f=0.5
Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
rho = DichteMatrix(nh,nc,nf,Hdilde)
rho_f=rho.ptrace(1)
WignerPlot = qutip.wigner(rho_f, np.linspace(-12,12,200), np.linspace(-12,12,200))

# plot the results

fig, axes = plt.subplots(2, 2,figsize=(12,10))
plt.rc ('axes', titlesize = 17)
cont0 = axes[0,0].contourf(xvec, xvec, W_coherent, 100)
fig.colorbar( plt.cm.ScalarMappable(),ax=axes[0,0])
lbl0 = axes[0,0].set_title(r' $f=0\gamma$')

cont1 = axes[0,1].contourf(xvec, xvec, W_thermal, 100)
fig.colorbar( plt.cm.ScalarMappable(),ax=axes[0,1])

lbl1 = axes[0,1].set_title(r' $f=0.1\gamma$')

cont0 = axes[1,0].contourf(np.linspace(-10,10,200), np.linspace(-10,10,200), W_fock, 100,colobar=True)
fig.colorbar( plt.cm.ScalarMappable(),ax=axes[1,0])

lbl2 = axes[1,0].set_title(r' $f=0.3\gamma$')

cont3 = axes[1,1].contourf(np.linspace(-12,12,200), np.linspace(-12,12,200), WignerPlot, 100,colobar=True)

lbl2 = axes[1,1].set_title(r' $f=0.5\gamma$')


fig.colorbar( plt.cm.ScalarMappable())
plt.show()

#############################################################################################################################################################################

xvec = np.linspace(-5,5,200)
f=0
nh=0.002
Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
rho = DichteMatrix(nh,nc,nf,Hdilde)
rho_f=rho.ptrace(1)
W_coherent = qutip.wigner(rho_f, xvec, xvec)

nh=4
Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
rho = DichteMatrix(nh,nc,nf,Hdilde)
rho_f=rho.ptrace(1)
W_thermal = qutip.wigner(rho_f, np.linspace(-4,4,200), np.linspace(-4,4,200))


nh=11
Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
rho = DichteMatrix(nh,nc,nf,Hdilde)
rho_f=rho.ptrace(1)
W_fock = qutip.wigner(rho_f, xvec, xvec)
#print("wignerplot an stelle 9.9 ",W_coherent[2][0.0])
nh=195
Hdilde=Hamilton(omega_1,proj_1,omega_2,proj_2,omega_3,proj_3,h,omega_f,a,f,g,omega_d)
rho = DichteMatrix(nh,nc,nf,Hdilde)
rho_f=rho.ptrace(1)
WignerPlot = qutip.wigner(rho_f, xvec, xvec)

# plot the results

fig, axes = plt.subplots(2, 2,figsize=(12,10))

cont0 = axes[0,0].contourf(xvec, xvec, W_coherent, 100)
fig.colorbar( plt.cm.ScalarMappable(),ax=axes[0,0])
lbl0 = axes[0,0].set_title(r'$n_h=0.002$')

cont1 = axes[0,1].contourf(xvec, xvec, W_thermal, 100)
fig.colorbar( plt.cm.ScalarMappable(),ax=axes[0,1])

lbl1 = axes[0,1].set_title(r'$n_h=4$')

cont0 = axes[1,0].contourf(xvec, xvec, W_fock, 100,colobar=True)
fig.colorbar( plt.cm.ScalarMappable(),ax=axes[1,0])

lbl2 = axes[1,0].set_title(r'$n_h=11$')

cont3 = axes[1,1].contourf(xvec, xvec, WignerPlot, 100,colobar=True)

lbl2 = axes[1,1].set_title(r'$n_h=195$')


fig.colorbar( plt.cm.ScalarMappable())
plt.show()

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




########################################################################################################################################################################
g_list=[]

Energie_VS_g=[]
for i in range(200):
    #list_temp=[]
    #list_temp=Diverse_Loups.EnergieCalculator(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,omega_2,V,omega_1,omega_d,proj_2,omega_f,c_op_list)
    g_list.append(i/100)  #Erstellt eine Liste mit Wären von g 
    #Energie_VS_g.append(list_temp)

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
anzahl=200
step=200/anzahl
nh2=0.001
nh_list2=[]
Entropy=[]
Entropy_tot=[]
for i in range(anzahl):
    list_temp=[]
    list_temp=Diverse_Loups.Entropy(nh2,Trans_12,a, kb,h,g,H_free,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,0,nph)
    #g_list.append(i/100)  #Erstellt eine Liste mit Wären von g 
    Entropy.append(list_temp)
    
    nh_list2.append(nh2)
    nh2=nh2+step
    Entropy_tot.append(Diverse_Loups.Entropy(nh2,Trans_12,a, kb,h,g,H_free,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,0,nph)[3])


nh2=0.001
nh_list3=[]
Entropy2=[]
Entropy_tot2=[]
for i in range(anzahl):
    list_temp2=[]
    list_temp2=Diverse_Loups.Entropy(nh2,Trans_12,a, kb,h,g,H_free,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,0.3,nph)
    #g_list.append(i/100)  #Erstellt eine Liste mit Wären von g 
    Entropy2.append(list_temp2)
    nh_list3.append(nh2)
    nh2=nh2+step
    
    Entropy_tot2.append(Diverse_Loups.Entropy(nh2,Trans_12,a, kb,h,g,H_free,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,0.3,nph)[3])


#Liste von Stings in floats konvertieren
#float_list2=list(np.float_(Energie_VS_g))
print(Entropy) 

#result=mesolve(H, rho0, tlist)
#print(D(c_op_list,rho)[3])


print("Die Temperatur des warmen Bades ist: ",T(omega_h,nh))
print("Die Temperatur des kalten Bades ist: ",T(omega_c,nc))


fig3, ax = plt.subplots()

ax.set_xlabel(r' $n_h$', fontsize=FZ)
ax.set_ylabel(r' $ \dot{\sigma} $' ,fontsize=FZ)

#plt.title(r' Entropy Production  rate vs $n_h$ ')
plt.plot(np.asarray(nh_list2)[:anzahl],np.asarray(Entropy)[:anzahl,0],label=r' $\frac{\langle J_h \rangle }{T_h}$',color='red')
plt.plot(np.asarray(nh_list2)[:anzahl],np.asarray(Entropy)[:anzahl,1],'-',label=r' $\frac{\langle J_c\rangle }{T_c}$',color='blue')
plt.plot(np.asarray(nh_list2)[:anzahl],np.asarray(Entropy)[:anzahl,2],label=r' $\frac{\langle J_{cav}\rangle }{T_{cav}}$',color='orange')
plt.plot(np.asarray(nh_list2)[:anzahl],np.asarray(Entropy_tot)[:anzahl],label=r' $\dot{\sigma}$',color='teal')
plt.plot(np.asarray(nh_list2)[:anzahl],np.asarray(Entropy_tot2)[:anzahl],color='teal',alpha=0.4)
plt.plot(np.asarray(nh_list3)[:anzahl],np.asarray(Entropy2)[:anzahl,0],'-',alpha=0.4,color='red')
plt.plot(np.asarray(nh_list3)[:anzahl],np.asarray(Entropy2)[:anzahl,1],'-',alpha=0.4,color='blue')
plt.plot(np.asarray(nh_list3)[:anzahl],np.asarray(Entropy2)[:anzahl,2],'-',alpha=0.4,color='orange')
#plt.plot(np.asarray(nh_list3)[:100],np.asarray(Entropy2)[:100,3],'--',label=r' $\frac{J_{cav}}{T_{cav}}$',color='orange')

legend = ax.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
fig3.set_figheight(9)
fig3.set_figwidth(13)
#Linien in plt
"""plt.axvline(x=2.6)
plt.axvline(x=2.6)
plt.axvline(x=5.5)
plt.axvline(x=0.17)
plt.axvline(x=20)
plt.axvline(x=1.7)"""



################################################################
"""""
g=14*kappa
f=0
f1=f
f_list=[]
nh=5
step=0.1
anzahl=200
for i in range(anzahl):
    f1=f1+step
    f_list.append(f1)


anzahl=200 #anzahl iterationen im loop

#f against the  power
P_list=Diverse_Loups.P3(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,omega_2,g,f,anzahl,step)
print("liste von P",P_list)


PowerPlot, ax = plt.subplots() 
ax.set_xlabel(r'$f$', fontsize=23)
ax.set_ylabel(r' Power', fontsize=15)
plt.title('power vs driven field-strenght')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(P_list)[:anzahl],label=r' Kurve')
#plt.show()



Energie_VS_f=Diverse_Loups.current(H_free, Trans_12, a, h,c_op_list,omega_d,omega_f ,proj_2,g,f,anzahl)

PundJ=[]
Energie_VS_f2=np.array(Energie_VS_f)
P_list2=np.array(P_list)
for i in range(anzahl):
    PundJ.append(Energie_VS_f2[i,0]+Energie_VS_f2[i,1]+Energie_VS_f2[i,2])



fig, ax = plt.subplots()
ax.set_xlabel(r' $\frac{f}{\gamma_h}$', fontsize=23)
ax.set_ylabel(r' Heat current or power ', fontsize=15)
plt.title('current/power vs driven field')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(P_list)[:anzahl],'--',label=r'$ \frac{P}{~ \gamma_h \omega_{h}}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energie_VS_f)[:anzahl,0],color='red',label=r' $\frac{J_h}{~ \gamma_h \omega_h}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energie_VS_f)[:anzahl,1],color='green',label=r' $\frac{J_c}{~\gamma_h \omega_h}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energie_VS_f)[:anzahl,2],color='yellow',label=r' $\frac{J_{cav}}{~\gamma_h \omega_{h}}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(PundJ)[:anzahl],'.',label=r'$~ \frac{P+J_c+J_h+J_{cav}}{~ \gamma_h \omega_{h}}$')
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')

#plt.show()


print(f_list)

"""""



###############################################################################################
"""""
nh5=5.5

g_list=[]
g=0
Energie_VS_g=[]
nc=nf=0.01

for i in range(200):
    g=g+1/12
    list_temp=[]
    list_temp=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh5,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2)
    g_list.append(g)  #Erstellt eine Liste mit Wären von g 
    Energie_VS_g.append(list_temp)

g=0
P_list=Diverse_Loups.P4(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,omega_2,g,f,anzahl)

f=0
g_list=[]
g=0
Energie_VS_g2=[]
for i in range(200):
    g=g+1/12
    list_temp=[]
    list_temp=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh5,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2)
    g_list.append(g)  #Erstellt eine Liste mit Wären von g 
    Energie_VS_g2.append(list_temp)

g=0
P_list2=Diverse_Loups.P4(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,omega_2,g,f,anzahl)



print("liste von P",P_list)
PowerPlot, ax = plt.subplots() 
ax.set_xlabel(r' g', fontsize=23)
ax.set_ylabel(r' Power', fontsize=15)
plt.title("Power vs  coupling-constant")
plt.plot(np.asarray(g_list)[:200],np.asarray(P_list)[:200],label=r' Kurve')

PundJ=[]
Energie_VS_g2=np.array(Energie_VS_g2)
P_list2=np.array(P_list)
for i in range(anzahl):
    PundJ.append(Energie_VS_g2[i,0]+Energie_VS_g2[i,1]+Energie_VS_g2[i,2]+P_list2[i])
fig, ax = plt.subplots()

ax.set_xlabel(r' $\frac{g}{\gamma_h}$', fontsize=23)
ax.set_ylabel(r' Heat current', fontsize=15)
plt.title('current/power vs g')
plt.plot(np.asarray(g_list)[:anzahl],np.asarray(P_list)[:anzahl],'.',label=r'$ \frac{P}{~ \gamma_h \omega_{h}}$')
plt.plot(np.asarray(g_list)[:anzahl],np.asarray(Energie_VS_g)[:anzahl,0],label=r' $\frac{J_h}{~ \gamma_h \omega_h}$')
plt.plot(np.asarray(g_list)[:anzahl],np.asarray(Energie_VS_g)[:anzahl,1],label=r' $\frac{J_c}{~\gamma_h \omega_h}$')
plt.plot(np.asarray(g_list)[:anzahl],np.asarray(Energie_VS_g)[:anzahl,2],label=r' $\frac{J_{cav}}{~\gamma_h \omega_{h}}$')
plt.plot(np.asarray(g_list)[:anzahl],np.asarray(P_list2)[:anzahl],'--',label=r'$ \frac{P}{~ \gamma_h \omega_{h}}$')
plt.plot(np.asarray(g_list)[:anzahl],np.asarray(Energie_VS_g2)[:anzahl,0],'--',label=r' $\frac{J_h}{~ \gamma_h \omega_h}$')
plt.plot(np.asarray(g_list)[:anzahl],np.asarray(Energie_VS_g2)[:anzahl,1],'--',label=r' $\frac{J_c}{~\gamma_h \omega_h}$')
plt.plot(np.asarray(g_list)[:anzahl],np.asarray(Energie_VS_g2)[:anzahl,2],'--',label=r' $\frac{J_{cav}}{~\gamma_h \omega_{h}}$')

plt.plot(np.asarray(g_list)[:anzahl],np.asarray(PundJ)[:anzahl],'*',label=r'$~ \frac{P+J_c+J_h+J_{cav}}{~ \gamma_h \omega_{h}}$')
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')



"""

""""

step=10
anzahl=100
##############################
#nh und power 
f=0
nh=0
g=14*kappa

omega_f=30
omega_h=100
Energie_VS_nh=[]
for i in range(anzahl):
    list_temp=[]
    list_temp=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2)
    nh=nh+step
    Energie_VS_nh.append(list_temp)
nh=0
P_list3,nh_list3=Diverse_Loups.P5(g,H_free, Trans_12, Trans_13, Trans_23, a,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2, anzahl,step)

print(P_list3,nh_list3)

Delta1=Delta2=0
f=1
nh=0.1
g=14*kappa
n_list=[]
nh_list=[]
omega_f=30
omega_h=100
Energie_VS_nh4=[]
for i in range(anzahl):
   
    list_temp=[]
    list_temp=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2)
     
    Energie_VS_nh4.append(list_temp)
    nh_list.append(nh)
    nh=nh+step
nh=0
P_list4,nh_list3=Diverse_Loups.P5(g,H_free, Trans_12, Trans_13, Trans_23, a,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2, anzahl,step)

print(P_list4,nh_list3)
PundJ=[]
Energie_VS_nh4=np.array(Energie_VS_nh4)
P_list2=np.array(P_list)
for i in range(anzahl):
    PundJ.append(Energie_VS_nh4[i,0]+Energie_VS_nh4[i,1]+Energie_VS_nh4[i,2]+P_list4[i])



    
    
fig, ax = plt.subplots()
ax.set_xlabel(r' $n_h$', fontsize=23)
ax.set_ylabel(r' Heat current', fontsize=15)
plt.title('current/power vs n_h')
plt.plot(np.asarray(nh_list3)[:anzahl],np.asarray(P_list3)[:anzahl],label=r'$ \frac{P}{~ \gamma_h \omega_{h}}$')
plt.plot(np.asarray(nh_list3)[:anzahl],np.asarray(Energie_VS_nh)[:anzahl,0],color='red',label=r' $\frac{J_h}{~ \gamma_h \omega_h}$')
plt.plot(np.asarray(nh_list3)[:anzahl],np.asarray(Energie_VS_nh)[:anzahl,1],color='green',label=r' $\frac{J_c}{~\gamma_h \omega_h}$')
plt.plot(np.asarray(nh_list3)[:anzahl],np.asarray(Energie_VS_nh)[:anzahl,2],color='orange',label=r' $\frac{J_{cav}}{~\gamma_h \omega_{h}}$')
plt.plot(np.asarray(nh_list3)[:anzahl],np.asarray(PundJ)[:anzahl],'.',color='black',label=r'$~ \frac{P+J_c+J_h+J_{cav}}{~ \gamma_h \omega_{h}}$')
plt.plot(np.asarray(nh_list3)[:anzahl],np.asarray(P_list4)[:anzahl],'--',label=r'$ \frac{P}{~ \gamma_h \omega_{h}}$')
plt.plot(np.asarray(nh_list3)[:anzahl],np.asarray(Energie_VS_nh4)[:anzahl,0],'--',color='red',label=r' $\frac{J_h}{~ \gamma_h \omega_h}$')
plt.plot(np.asarray(nh_list3)[:anzahl],np.asarray(Energie_VS_nh4)[:anzahl,1],'--',color='green',label=r' $\frac{J_c}{~\gamma_h \omega_h}$')
plt.plot(np.asarray(nh_list3)[:anzahl],np.asarray(Energie_VS_nh4)[:anzahl,2],'--',color='orange',label=r' $\frac{J_{cav}}{~\gamma_h \omega_{h}}$')
#plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(n_list)[:anzahl,0],'*',color='black',label=r'numerical solved EqoM')

legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('white')


    


#legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
#legend.get_frame().set_facecolor('C0')
#Linien in plt
'''plt.axvline(x=2.6)
plt.axvline(x=2.6)
plt.axvline(x=5.5)
plt.axvline(x=0.17)
plt.axvline(x=20)
plt.axvline(x=1.7)

'''

#Just J_cav



step=1
nh_list4=[]
for i in range(200):
    nh_list4.append(nh)
    nh+=step
    

nh=0
g=14*kappa
anzahl=200
omega_f=30
omega_h=100
Energie_VS_nh=[]
for i in range(anzahl):
    nh=nh+step
    list_temp=[]
    list_temp=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2)
      #Erstellt eine Liste mit Wären von g 
    Energie_VS_nh.append(list_temp)
    Energie_VS_nh=Energie_VS_nh

P_list3,nh_list3=Diverse_Loups.P5(g,H_free, Trans_12, Trans_13, Trans_23, a,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2, anzahl,step)

print(P_list3,nh_list3)

f=0
nh=0
Energie_VS_nh2=[]
for i in range(anzahl):
    nh=nh+step
    list_temp2=[]
    list_temp2=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2)
      #Erstellt eine Liste mit Wären von g 
    Energie_VS_nh2.append(list_temp2)
    Energie_VS_nh2=Energie_VS_nh2



fig, ax = plt.subplots()
ax.set_xlabel(r' $n_h$', fontsize=23)
ax.set_ylabel(r' Heat current', fontsize=15)
plt.title(r' current $J_{cav}$ VS $n_h$ Doble threshold behaviour')
plt.plot(np.asarray(nh_list4)[:anzahl],np.asarray(Energie_VS_nh)[:anzahl,2],label=r' $\frac{J_{cav}}{~\gamma_h \omega_{h}}$')
plt.plot(np.asarray(nh_list4)[:anzahl],np.asarray(Energie_VS_nh2)[:anzahl,2],'--',label=r' $\frac{J_{cav}}{~\gamma_h \omega_{h}}$')



print('Dichtematrix von nh=0',DichteMatrix(0,0.02,0.02,Hdilde))
plt.show()
 

"""

#Projector occupation probability

Delta1=Delta2=0
gamma_h = gamma_c = 1

nc = ncav = 0.0
g=14*kappa

Delta1=0
Delta2=0
anzahl=100
nh=0.1
nc=ncav=nf=0.1
n_list=[]
nh_list=[]
f=0
Photonnumber_list=[]
nh2 = np.linspace(0, 70, 100)

nh_list=[]
Trace_list=[]
nh=0.1 #set nh again to zero
step=0.7

H_free=omega_1*proj_1+h*omega_2*proj_2+h*omega_3*proj_3+h*omega_f*a.dag()*a

H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

V=f*a.dag()+f*a #das got glaub nid

H=H_free+H_int


print(a, a.dag(), a*a*a.dag()-a*a.dag()*a)

#Hdilde=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 

Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a) 
n_list=((Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , f , nh, ncav , nc, gamma_c, gamma_h, g , kappa, anzahl,step,'nh')))  
for i in range(anzahl):
    
    #if (isinstance(n_list[i], complex) or n_list[i]<n_list[i-1]-10):
    #    n_list[i]=n_list[i-1]
    nh_list.append(nh)
    
    Trace_list_temp=Diverse_Loups.ProjectorP_mit_faktor(Trans_12,h,g,nh,proj_1,proj_2,proj_3,V ,omega_2,omega_1,omega_d,omega_f,a ,nc,nf,gamma_h,gamma_c,kappa,A1,A2,A3,A4,A5,A6)
    Trace_list.append(Trace_list_temp)
    ist_temp=[]
    list_temp=Diverse_Loups.Photonnumber(nh,a,proj_1,proj_2,proj_3,Hdilde,nc,ncav,gamma_h,gamma_c,kappa,A1,A2,A3,A4,A5,A6,omega_d,omega_f,omega_1,omega_2,H_int,f)
    #g_list.append(i/100)  #Erstellt eine Liste mit Wären von g 
    Photonnumber_list.append(list_temp)
    #print(n_list[i],Photonnumber_list[i])
    nh=nh+step



fig2, ax = plt.subplots()
ax.set_xlabel(r' $n_h$', fontsize=FZ)
#ax.set_ylabel('probability')
#plt.title('stationary atomic population')
    
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(n_list[3])[:anzahl],'--',color='blue')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(n_list[1])[:anzahl],'--',color='green')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(n_list[2])[:anzahl],'--',color='orange')
plt.plot(np.asarray(nh_list)[:100],np.asarray(Trace_list)[:100,0],color='green',label=r'$P_1$')
plt.plot(np.asarray(nh_list)[:100],np.asarray(Trace_list)[:100,1],color='orange',label=r'$P_2$')
plt.plot(np.asarray(nh_list)[:100],np.asarray(Trace_list)[:100,2],color='blue',label=r'$P_3$')

legend = ax.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
fig2.set_figheight(9)
fig2.set_figwidth(13)







#plt.show()

##############################################################################################
f=0

anzahl =100
step=0.5/anzahl
nc=ncav=nf=0
nh=11

Cav_Power_mit_nh=[]
Cav_Power_ohne_nh=[]
f_list=[]
Energy_VS_f=[]
Energy_VS_f_nh=[]
J_h_f_list=[]
J_c_f_list=[]
J_cav_f_list=[]
J_tot_f_list=[]
P_ana=[]
Power=[]
Tot=[]
f_Line=np.linspace(0, 0.5, anzahl)
JH_f=Diverse_Loups.High_f_Approx(kappa,f_Line,gamma_c,gamma_h,nc,nh,g,nf)[0]
JC_f=omega_c/omega_h*Diverse_Loups.High_f_Approx(kappa,f_Line,gamma_c,gamma_h,nc,nh,g,nf)[1]
JCAV_f=omega_f/omega_h*Diverse_Loups.High_f_Approx(kappa,f_Line,gamma_c,gamma_h,nc,nh,g,nf)[2]
POWER_f=omega_f/omega_h*Diverse_Loups.High_f_Approx(kappa,f_Line,gamma_c,gamma_h,nc,nh,g,nf)[5]
n_f_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , 0 , nh, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"f"))
f_Line =np.linspace(0,0.5,100)
P_Pre2=Diverse_Loups.High_f_Approx(kappa,f_Line,gammac,gammah,nc,nh,g,ncav)[0]*(-1)*omega_f/omega_h
f_LineZ =np.linspace(0,step*30,100)
P_PreZ=Diverse_Loups.High_f_Approx(kappa,f_LineZ,gammac,gammah,nc,nh,g,ncav)[0]*(-1)*omega_f/omega_h
for  i in range(anzahl):


    P_ana.append(-30/omega_h*2*f*(n_f_list[4][i]))
    J_cav_f_list.append(30/omega_h*2*kappa*(nf-n_f_list[0][i]))
    J_h_f_list.append(150/omega_h*(nh*n_f_list[1][i]-(nh+1)*n_f_list[3][i]))
    J_c_f_list.append(120/omega_h*(nc*n_f_list[2][i]-(nc+1)*n_f_list[3][i]))
    J_tot_f_list.append(J_c_f_list[i]+J_h_f_list[i]+J_cav_f_list[i]+P_ana[i])
    list_temp=[]
    list_temp=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h)
    Energy_VS_f.append(list_temp)
    list_temp_nh=[]
    list_temp_nh=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, 200 ,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h)
    Energy_VS_f_nh.append(list_temp_nh)
    Power.append(Diverse_Loups.P2(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,f,omega_2,g)/omega_h)
    Tot.append(list_temp[0]+list_temp[1]+list_temp[2]+Power[i])
    Cav_Power_mit_nh.append(list_temp[2]+Power[i])
    #Cav_Power_ohne_nh.append(list_temp[2])
    f_list.append(f)
    f+=step




fig1, ax = plt.subplots()
ax.set_xlabel(r' $\frac{f}{\gamma}$', fontsize=FZ)
#ax.set_ylabel(r' $\langle J \rangle $', fontsize=21)
plt.title('')


plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Power)[:anzahl],'-',color='green',linewidth=3,label=r' $\frac{\langle P\rangle}{~\gamma_h \omega_{h}}$')

ax.plot(np.asarray(f_list)[:anzahl],np.asarray(J_h_f_list)[:anzahl],'--',color='red')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Tot)[:anzahl],'-',linewidth=3,color='black',label=r'$\frac{\Sigma_i \langle J_{i}\rangle +\langle P \rangle}{~\gamma_h \omega_{h}}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(J_c_f_list)[:anzahl],'--',color='blue')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(J_cav_f_list)[:anzahl],'--',color='orange')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(P_ana)[:anzahl],'--',color='green')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energy_VS_f)[:anzahl,1],'-',color='blue',linewidth=3,label=r' $\frac{\langle J_{c}\rangle }{~\gamma_h \omega_{h}}$')
ax.plot(np.asarray(f_list)[:anzahl],np.asarray(Energy_VS_f)[:anzahl,0],'-',color='red',linewidth=3,label=r' $\frac{\langle J_{h}\rangle }{~\gamma_h \omega_{h}}$')
plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energy_VS_f)[:anzahl,2],'-',color='orange',linewidth=3,label=r' $\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}}$')
#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energy_VS_f_nh)[:anzahl,1],'*',color='blue',linewidth=3)
#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energy_VS_f_nh)[:anzahl,0],'*',color='red',linewidth=3)
#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(Energy_VS_f_nh)[:anzahl,2],'*',color='orange',linewidth=3)
plt.plot(np.asarray(f_list)[:anzahl],JC_f,'-',color='blue')#,label=r'$\frac{\langle J_{c}\rangle }{~\gamma_h \omega_{h}}$ Eq 32')
plt.plot(np.asarray(f_list)[:anzahl],JH_f,'-',color='red')#,label=r'$\frac{\langle J_{h}\rangle }{~\gamma_h \omega_{h}} $ Eq 31')
plt.plot(f_Line,JCAV_f,color='orange')#,label=r'$\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}} $ ')
plt.plot(f_Line,POWER_f,color='green')#,label=r'$\frac{\langle P \rangle }{~\gamma_h \omega_{h}} $ ')

""""
axins= ax.inset_axes([0.1,0.75,0.5,0.25])
axins.plot(np.asarray(f_list)[:anzahl], np.asarray(Energy_VS_f)[:anzahl,0],color='red')
axins.plot(np.asarray(f_list)[:anzahl], np.asarray(J_h_f_list)[:anzahl],'--',color='red')
mark_inset(ax, axins, loc1=2, loc2=4, fc="None", ec="grey",color='grey')
"""
legend = ax.legend(loc='upper right', shadow=True, fontsize=23)
legend.get_frame().set_facecolor('white')
fig1.set_figheight(9)
fig1.set_figwidth(13)
plt.axvline(x=0)
plt.axvline(x=0.1)
plt.axvline(x=0.3)
plt.axvline(x=0.5)
###############################################################################################
plt.show()
###########################################################################################################################


fig5, ax = plt.subplots()
ax.set_xlabel(r' $\frac{f}{\gamma}$', fontsize=FZ)
#ax.set_ylabel(r' $\langle J \rangle $', fontsize=21)
plt.title('')

ax.plot(np.asarray(f_list)[:anzahl], np.asarray(Cav_Power_mit_nh)[:anzahl],'-',color='purple',label=r'$\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}}+\frac{\langle P \rangle }{~\gamma_h \omega_{h}}$ ')
ax.plot(f_Line,P_Pre2,label=r' $\frac{\langle P \rangle }{~\gamma_h \omega_{h}}$  ')
ax.plot(np.asarray(f_list)[:anzahl],np.asarray(Energy_VS_f)[:anzahl,2],'-',color='orange',linewidth=3,label=r' $\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}}$')
ax.plot(np.asarray(f_list)[:anzahl],np.asarray(Power)[:anzahl],'-',color='green',linewidth=3,label=r' $\frac{\langle P\rangle}{~\gamma_h \omega_{h}}$')

axins= ax.inset_axes([0.53,0.53,0.47,0.47])
axins.plot(np.asarray(f_list)[:30], np.asarray(Cav_Power_mit_nh)[:30],color='purple')
axins.plot(f_LineZ, P_PreZ)
axins.plot(np.asarray(f_list)[:30], np.asarray(Energy_VS_f)[:30,2],color='orange')
axins.plot(np.asarray(f_list)[:30], np.asarray(Power)[:30],color='green')
#axins.plot(np.asarrayf_Line, P_Pre2)

mark_inset(ax, axins, loc1=2, loc2=4, fc="None", ec="red",color='red')


legend = ax.legend(loc='lower left', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
fig5.set_figheight(9)
fig5.set_figwidth(13)
plt.draw()
plt.show()
Delta1=Delta2=0
gamma_h = gamma_c = 1

nc = nf=ncav = 0.0

kappa=0.2

step=1
Delta1=0
Delta2=0
anzahl=200
g=2.8
nc=ncav=nf=0

#g=0
n_list=[]
nh_list=[]
f=0.2
Photonnumber_list=[]
nh2 = np.linspace(0, 70, 100)

nh3_list=[]
nh_list=[]
Trace_list=[]
nh=0.0 #set nh again to zero
Anal=[]
J_h_list=[]
J_c_list=[]
J_cav_list=[]
J_tot_list=[]

J_h_f_list=[]
J_c_f_list=[]
J_cav_f_list=[]  
J_tot_f_list=[]

J_h_f2_list=[]
J_c_f2_list=[]
J_cav_f2_list=[]  
J_tot_f2_list=[]

J_h_f3_list=[]
J_c_f3_list=[]
J_cav_f3_list=[]  
J_tot_f3_list=[]

nh_Line=np.linspace(0, 200, 70)
f_Line=np.linspace(0, 1, 70)
JH=Diverse_Loups.High_f_Approx(kappa,0.2,gamma_c,gamma_h,nc,nh_Line,g,nf)[0]
JC=omega_c/omega_h*Diverse_Loups.High_f_Approx(kappa,0.2,gamma_c,gamma_h,nc,nh_Line,g,nf)[1]
JCAV=omega_f/omega_h*Diverse_Loups.High_f_Approx(kappa,0.2,gamma_c,gamma_h,nc,nh_Line,g,nf)[2]
P=omega_f/omega_h*Diverse_Loups.High_f_Approx(kappa,0.2,gamma_c,gamma_h,nc,nh_Line,g,nf)[5]
JH2=Diverse_Loups.High_f_Approx(kappa,28,gamma_c,gamma_h,nc,nh_Line,g,nf)[0]
JC2=omega_c/omega_h*Diverse_Loups.High_f_Approx(kappa,28,gamma_c,gamma_h,nc,nh_Line,g,nf)[1]
P2=omega_f/omega_h*Diverse_Loups.High_f_Approx(kappa,28,gamma_c,gamma_h,nc,nh_Line,g,nf)[5]
JCAV2=omega_f/omega_h*Diverse_Loups.High_f_Approx(kappa,28,gamma_c,gamma_h,nc,nh_Line,g,nf)[2]
JH3=Diverse_Loups.High_f_Approx(kappa,0.7,gamma_c,gamma_h,nc,nh_Line,g,nf)[0]
JC3=omega_c/omega_h*Diverse_Loups.High_f_Approx(kappa,0.7,gamma_c,gamma_h,nc,nh_Line,g,nf)[1]
P3=omega_f/omega_h*Diverse_Loups.High_f_Approx(kappa,0.7,gamma_c,gamma_h,nc,nh_Line,g,nf)[5]
JCAV3=omega_f/omega_h*Diverse_Loups.High_f_Approx(kappa,0.7,gamma_c,gamma_h,nc,nh_Line,g,nf)[2]

Energie_vs_nh=[]
Cav_Power_mit_f=[]
Cav_Power_ohne_f=[]
Total=[]
Energie_vs_nh_f=[]
n2_list=[]
nh=0.0
Power_ohne_f=[]
Power_mit_f=[]
Total=[]
P_ana=[]
P_ana_f=[]
n_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , 0 , nh, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"nh"))
n_f_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , f , nh, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"nh"))
n_f2_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , 0.7 , nh, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"nh"))
n_f3_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , 20 , nh, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"nh"))
for i in range(anzahl):
   
    
    list_temp=[]
    list_temp=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,0,omega_f,omega_2,omega_h)
    
    Energie_vs_nh.append(list_temp)
    Power_ohne_f.append(Diverse_Loups.P2(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,0,omega_2,g)/omega_h)
    
    
    
    
    
   
    
    list_temp2=[]
    list_temp2=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h)
    
    Energie_vs_nh_f.append(list_temp2)
    Power_mit_f.append(Diverse_Loups.P2(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,f,omega_2,g)/omega_h)
    Cav_Power_mit_f.append(Power_mit_f[i]+ list_temp2[2])
    Cav_Power_ohne_f.append(Power_ohne_f[i]+ list_temp[2])
    #Total.append(list_temp[0]+list_temp[1]+list_temp[2]+Power_ohne_f[i])
    Total.append(list_temp2[0]+list_temp2[1]+list_temp2[2]+Power_mit_f[i])
    nh_list.append(nh)
    nh3_list.append(nh)
    
    J_cav_list.append(30/omega_h*2*kappa*(nf-n_list[0][i]))
    J_h_list.append(150/omega_h*(nh*n_list[1][i]-(nh+1)*n_list[3][i]))
    J_c_list.append(120/omega_h*(nc*n_list[2][i]-(nc+1)*n_list[3][i]))
    J_tot_list.append(J_c_list[i]+J_h_list[i]+J_cav_list[i])
    P_ana.append(-2*30/omega_h*f*(n_list[4][i]))
    J_cav_f_list.append(30/omega_h*2*kappa*(nf-n_f_list[0][i]))
    J_h_f_list.append(150/omega_h*(nh*n_f_list[1][i]-(nh+1)*n_f_list[3][i]))
    J_c_f_list.append(120/omega_h*(nc*n_f_list[2][i]-(nc+1)*n_f_list[3][i]))
    J_tot_f_list.append(J_c_list[i]+J_h_list[i]+J_cav_list[i])
    P_ana_f.append(-2*30/omega_h*f*(n_f_list[4][i]))
    
    J_cav_f2_list.append(30/omega_h*2*kappa*(nf-n_f2_list[0][i]))
    J_h_f2_list.append(150/omega_h*(nh*n_f2_list[1][i]-(nh+1)*n_f2_list[3][i]))
    J_c_f2_list.append(120/omega_h*(nc*n_f2_list[2][i]-(nc+1)*n_f2_list[3][i]))
    
    J_cav_f3_list.append(30/omega_h*2*kappa*(nf-n_f2_list[0][i]))
    J_h_f3_list.append(150/omega_h*(nh*n_f3_list[1][i]-(nh+1)*n_f3_list[3][i]))
    J_c_f3_list.append(120/omega_h*(nc*n_f3_list[2][i]-(nc+1)*n_f3_list[3][i]))
    
    

    nh+=step




fig2, (ax1,ax2) = plt.subplots(1,2)
ax1.tick_params( labelsize=tz)
ax2.tick_params( labelsize=tz)
ax1.set_xlabel(r' $n_h$', fontsize=FZ)
ax2.set_xlabel(r' $n_h$', fontsize=FZ)
#ax.set_ylabel('heat  current')
ax1.set_title('a)')
ax2.set_title('b)')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_h_list)[:anzahl],'--',color='red')
#plt.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_tot_list)[:anzahl],'--',color='black',label=r'$\frac{J_{tot}}{~\gamma_h \omega_{h}}$analytisch')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_c_list)[:anzahl],'--',color='blue')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_cav_list)[:anzahl],'--',color='orange')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_h_f_list)[:anzahl],'--',color='red',alpha=0.5)
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_tot_f_list)[:anzahl],'--',color='black',alpha=0.5)
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_c_f_list)[:anzahl],'--',color='blue',alpha=0.5)
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_cav_f_list)[:anzahl],'--',color='orange',alpha=0.5)
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh)[:anzahl,1],'-',color='blue',linewidth=3,label=r' $\frac{\langle J_{c}\rangle }{~\gamma_h \omega_{h}}$')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh)[:anzahl,0],'-',color='red',linewidth=3,label=r' $\frac{\langle J_{h}\rangle }{~\gamma_h \omega_{h}}$')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh)[:anzahl,2],'-',color='orange',linewidth=3,label=r' $\frac{\langle J_{cav} \rangle  }{~\gamma_h \omega_{h}}$')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh_f)[:anzahl,0],'-',color='red',linewidth=3,alpha=0.5,)
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh_f)[:anzahl,2],'-',color='orange',linewidth=3,alpha=0.5)
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh_f)[:anzahl,1],'-',color='blue',linewidth=3,alpha=0.5)
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Total)[:anzahl],'-',color='black',label=r'$\frac{\Sigma_i \langle J_{i}\rangle +\langle P \rangle }{~\gamma_h \omega_{h}}$')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Power_ohne_f)[:anzahl],'-',color='green',linewidth=3,label=r' $\frac{\langle P \rangle }{~\gamma_h \omega_{h}}$')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Power_mit_f)[:anzahl],'-',color='green',linewidth=3,alpha=0.5)
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(P_ana)[:anzahl],'--',color='green')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(P_ana_f)[:anzahl],'--',color='green',alpha=0.5)
ax2.plot(nh_Line,JC,'-',color='blue',label=r'$\frac{\langle J_{c}\rangle }{~\gamma_h \omega_{h}} $ ')
ax2.plot(nh_Line,JH,'-',color='red',label=r'$\frac{\langle J_{h}\rangle }{~\gamma_h \omega_{h}} $ ')
ax2.plot(nh_Line,JC2,'-',color='blue',alpha=0.3)
ax2.plot(nh_Line,JH2,'-',color='red',alpha=0.3)
ax2.plot(nh_Line,JC3,'-',color='blue',alpha=0.5)
ax2.plot(nh_Line,JH3,'-',color='red',alpha=0.5)

"""""
ax2.plot(nh_Line,P,'-',color='red',label=r'$\frac{\langle P \rangle }{~\gamma_h \omega_{h}} $ Eq 32')
ax2.plot(nh_Line,JCAV,'-',color='blue',label=r'$\frac{\langle J_{cav} \rangle }{~\gamma_h \omega_{h}} $ Eq 32')
ax2.plot(nh_Line,P2,'-',color='red',alpha=0.3)
ax2.plot(nh_Line,JCAV2,'-',color='blue',alpha=0.3)
ax2.plot(nh_Line,P3,'-',color='red',alpha=0.5)
ax2.plot(nh_Line,JCAV3,'-',color='red',alpha=0.5)
"""
ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_h_f_list)[:anzahl],'--',color='red',alpha=1)
ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_c_f_list)[:anzahl],'--',color='blue',alpha=1)
#ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_cav_f_list)[:anzahl],'--',color='orange',alpha=1)

ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_h_f2_list)[:anzahl],'--',color='red',alpha=0.5)
ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_c_f2_list)[:anzahl],'--',color='blue',alpha=0.5)
#ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_cav_f2_list)[:anzahl],'--',color='orange',alpha=0.5)

ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_h_f3_list)[:anzahl],'--',color='red',alpha=0.3)
ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_c_f3_list)[:anzahl],'--',color='blue',alpha=0.3)
#ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_cav_f3_list)[:anzahl],'--',color='orange',alpha=0.3)



legend = ax1.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
fig2.set_figheight(9)
fig2.set_figwidth(13)
legend = ax2.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
fig2.set_figheight(9)
fig2.set_figwidth(13)

ax1.axvline(x=0.1)
ax1.axvline(x=0.5)
ax1.axvline(x=11)
ax1.axvline(x=195)













f_Line =np.linspace(0,0.5,100)
#P_Pre=Diverse_Loups.High_f_Approx(kappa,f_Line,gammac,gammah,nc,nh,g,ncav)[0]*(-1)*omega_f/omega_h

nh_Line =np.linspace(0,200,100)
P_Pre_vs_nh=-1*Diverse_Loups.High_f_Approx(kappa,f,gammac,gammah,nc,nh_Line,g,ncav)[0]*omega_f/omega_h
FigJP, (ax,ax1) = plt.subplots(1,2)
ax.set_xlabel(r' $\frac{f}{\gamma}$', fontsize=FZ)
ax.set_title('a)')
ax1.set_title('b)')
#ax.set_ylabel(r' $\langle J \rangle $', fontsize=21)


ax.plot(np.asarray(f_list)[:anzahl], np.asarray(Cav_Power_mit_nh)[:anzahl],'-',color='purple',label=r'$\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}}+\frac{\langle P \rangle }{~\gamma_h \omega_{h}}$ ')
ax.plot(f_Line,P_Pre2,label=r' $\frac{\langle P \rangle }{~\gamma_h \omega_{h}}$  ')
ax.plot(np.asarray(f_list)[:anzahl],np.asarray(Energy_VS_f)[:anzahl,2],'-',color='orange',linewidth=3,label=r' $\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}}$')
ax.plot(np.asarray(f_list)[:anzahl],np.asarray(Power)[:anzahl],'-',color='green',linewidth=3,label=r' $\frac{\langle P\rangle}{~\gamma_h \omega_{h}}$')

axins= ax.inset_axes([0.60,0.60,0.40,0.40])
axins.plot(np.asarray(f_list)[:30], np.asarray(Cav_Power_mit_nh)[:30],color='purple')
axins.plot(f_LineZ, P_PreZ)
axins.plot(np.asarray(f_list)[:30], np.asarray(Energy_VS_f)[:30,2],color='orange')
axins.plot(np.asarray(f_list)[:30], np.asarray(Power)[:30],color='green')
#axins.plot(np.asarrayf_Line, P_Pre2)

mark_inset(ax, axins, loc1=2, loc2=4, fc="None", ec="grey",color='red')


#legend = ax.legend(loc='lower left', shadow=True, fontsize=fz)
#legend.get_frame().set_facecolor('white')
FigJP.set_figheight(9)
FigJP.set_figwidth(13)


ax1.set_xlabel(r' $n_h$', fontsize=FZ)
#ax.set_ylabel('heat  current')

ax1.plot(nh_Line,P_Pre_vs_nh,label=r' $\frac{\langle P \rangle }{~\gamma_h \omega_{h}}$ ')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Power_ohne_f)[:anzahl],'-',color='green',label=r' $\frac{\langle P \rangle }{~\gamma_h \omega_{h}}$')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh_f)[:anzahl,2],'-',color='orange',alpha=0.5)
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Cav_Power_mit_f)[:anzahl],'-',color='purple',alpha=0.5)
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Cav_Power_ohne_f)[:anzahl],'-',color='purple',alpha=0.5,label=r'$\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}}+\frac{\langle P \rangle }{~\gamma_h \omega_{h}}$ ')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Power_mit_f)[:anzahl],'-',color='green',alpha=0.5)
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh)[:anzahl,2],'-',color='orange',label=r' $\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}}$')
legend = ax1.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
FigJP.set_figheight(9)
FigJP.set_figwidth(13)
plt.show()

"""""
fig2, ax = plt.subplots()
ax.set_xlabel(r' $n_h$', fontsize=21)
ax.set_ylabel('heat  current')
plt.title('')
plt.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh_f)[:anzahl,0],'-',color='red',alpha=0.5,)
plt.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh_f)[:anzahl,2],'-',color='orange',alpha=0.5)
plt.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh_f)[:anzahl,1],'-',color='green',alpha=0.5)
plt.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Power_mit_f)[:anzahl],'-',color='blue',alpha=0.5)
plt.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh)[:anzahl,1],'-',color='green',label=r' $\frac{J_{c}}{~\gamma_h \omega_{h}}$')
plt.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh)[:anzahl,0],'-',color='red',label=r' $\frac{J_{h}}{~\gamma_h \omega_{h}}$')
plt.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Energie_vs_nh)[:anzahl,2],'-',color='orange',label=r' $\frac{J_{cav}}{~\gamma_h \omega_{h}}$')
plt.plot(np.asarray(nh3_list)[:anzahl],np.asarray(Power_ohne_f)[:anzahl],'-',color='blue',label=r' $\frac{P_{cav}}{~\gamma_h \omega_{h}}$')

legend = ax.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
fig3.set_figheight(9)
fig3.set_figwidth(13)

"""
#################################################################################################################

Delta1=Delta2=0
gamma_h = gamma_c = 1
g = 0
nc = ncav = 0.0





anzahl=100
step=2.8/anzahl
Delta1=0
Delta2=0


nc=nf=0

f=0.2
nh=5
Energie_vs_g=[]
Energie_vs_g_f=[]

n_list=[]
g_list=[]
J_h_list=[]

J_c_list=[]
J_cav_list=[]
J_tot_list=[]


J_h_f_list=[]

J_c_f_list=[]
J_cav_f_list=[]
J_tot_f_list=[]

Energie_vs_nh=[]

Tot=[]
Tot2=[]
P_list_f=[]
P_ana=[]
P_ana_f=[]
P_list=Diverse_Loups.P(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,0,anzahl,step,omega_h)
P_list_f=Diverse_Loups.P(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,f,anzahl,step,omega_h)
n_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , 0 , nh, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"g"))
n_f_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , f , nh, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"g"))

g_Line=np.linspace(0, 2.8, 100)

JH_g=Diverse_Loups.High_f_Approx(kappa,0.2,gamma_c,gamma_h,nc,nh,g_Line,nf)[0]
JC_g=omega_c/omega_h*Diverse_Loups.High_f_Approx(kappa,f,gamma_c,gamma_h,nc,nh,g_Line,nf)[1]
JH2_g=Diverse_Loups.High_f_Approx(kappa,0.2,gamma_c,gamma_h,nc,nh,g_Line,nf)[0]
JC2_g=omega_c/omega_h*Diverse_Loups.High_f_Approx(kappa,f,gamma_c,gamma_h,nc,nh,g_Line,nf)[1]
JCAV_g=omega_f/omega_h*Diverse_Loups.High_f_Approx(kappa,f,gamma_c,gamma_h,nc,nh,g_Line,nf)[2]
POWER_g=(omega_f/omega_h)*Diverse_Loups.High_f_Approx(kappa,f,gamma_c,gamma_h,nc,nh,g_Line,nf)[5]
TOT=JH_g+JC_g
TOT2=JH2_g+JC2_g+JCAV_g+POWER_g
for i in range(anzahl):
   
    
    list_temp=[]
    list_temp=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,0,omega_f,omega_2,omega_h)
    Energie_vs_g.append(list_temp)
    
    list_temp2=[]
    list_temp2=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h)
    
    Energie_vs_g_f.append(list_temp2)
    
    Tot.append(list_temp[0]+list_temp[1]+list_temp[2]+P_list[i])
    Tot2.append(list_temp2[0]+list_temp2[1]+list_temp2[2]+P_list_f[i])
    
    
    J_cav_list.append(30/omega_h*2*kappa*(nf-n_list[0][i]))
    J_h_list.append(150/omega_h*(nh*n_list[1][i]-(nh+1)*n_list[3][i]))
    J_c_list.append(120/omega_h*(nc*n_list[2][i]-(nc+1)*n_list[3][i]))
    P_ana.append(-2*30/omega_h*0*(n_f_list[4][i]))
    J_tot_list.append(J_c_list[i]+J_h_list[i]+J_cav_list[i])
    
    J_cav_f_list.append(30/omega_h*2*kappa*(nf-n_f_list[0][i]))
    J_h_f_list.append(150/omega_h*(nh*n_f_list[1][i]-(nh+1)*n_f_list[3][i]))
    J_c_f_list.append(120/omega_h*(nc*n_f_list[2][i]-(nc+1)*n_f_list[3][i]))
    P_ana_f.append(-2*30/omega_h*f*(n_f_list[4][i]))
    J_tot_f_list.append(J_c_f_list[i]+J_h_f_list[i]+J_cav_f_list[i]+P_ana_f[i])
    
    
    
    g=g+step
    g_list.append(g)


Gplot, (ax1, ax2) = plt.subplots(1,2)
ax1.tick_params( labelsize=tz)
ax2.tick_params( labelsize=tz)
ax2.set_xlabel(r' $\frac{g}{\gamma}$', fontsize=FZ)
ax1.set_xlabel(r' $\frac{g}{\gamma}$', fontsize=FZ)
ax1.set_title('a)')
ax2.set_title('b)')
#ax2.set_ylabel('heat current')


ax1.plot(np.asarray(g_list)[:anzahl],np.asarray(P_ana)[:anzahl],'--',color='green')
ax2.plot(np.asarray(g_list)[:anzahl],np.asarray(P_ana_f)[:anzahl],'--',color='green')
ax1.plot(np.asarray(g_list)[:anzahl],np.asarray(Energie_vs_g)[:anzahl,1],'-',color='blue',linewidth=3,label=r' $\frac{\langle J_{c}\rangle }{~\gamma_h \omega_{h}}$')
ax1.plot(np.asarray(g_list)[:anzahl],np.asarray(Energie_vs_g)[:anzahl,0],'-',color='red',linewidth=3,label=r' $\frac{\langle J_{h}\rangle }{~\gamma_h \omega_{h}}$')
ax1.plot(np.asarray(g_list)[:anzahl],np.asarray(Energie_vs_g)[:anzahl,2],'-',color='orange',linewidth=3,label=r' $\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}}$')
ax2.plot(np.asarray(g_list)[:anzahl],np.asarray(Energie_vs_g_f)[:anzahl,0],'-',color='red',linewidth=3,label=r' $\frac{\langle J_{h}\rangle }{~\gamma_h \omega_{h}}$')
ax2.plot(np.asarray(g_list)[:anzahl],np.asarray(Energie_vs_g_f)[:anzahl,2],'-',color='orange',linewidth=3,label=r' $\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}}$')
ax2.plot(np.asarray(g_list)[:anzahl],np.asarray(Energie_vs_g_f)[:anzahl,1],'-',color='blue',linewidth=3,label=r' $\frac{\langle J_{c}\rangle }{~\gamma_h \omega_{h}}$')
#ax1.plot(np.asarray(g_list)[:anzahl],np.asarray(J_tot_list)[:anzahl],'--',color='black',label=r'$\frac{\Sigma_i \langle J_{i}\rangle +\langle P \rangle}{~\gamma_h \omega_{h}}$')
ax1.plot(np.asarray(g_list)[:anzahl],np.asarray(P_list)[:anzahl],'-',color='green',linewidth=3,label=r' $\frac{\langle P \rangle }{~\gamma_h \omega_{h}}$')
ax2.plot(np.asarray(g_list)[:anzahl],np.asarray(P_list_f)[:anzahl],'-',color='green',linewidth=3,label=r' $\frac{\langle P \rangle }{~\gamma_h \omega_{h}}$')
ax1.plot(np.asarray(g_list)[:anzahl],np.asarray(J_c_list)[:anzahl],'--',color='blue')
ax1.plot(np.asarray(g_list)[:anzahl],np.asarray(J_h_list)[:anzahl],'--',color='red')
ax1.plot(np.asarray(g_list)[:anzahl],np.asarray(J_cav_list)[:anzahl],'--',color='orange')

#ax2.plot(np.asarray(g_list)[:anzahl],np.asarray(J_tot_f_list)[:anzahl],'--',color='black')

ax2.plot(np.asarray(g_list)[:anzahl],np.asarray(J_c_f_list)[:anzahl],'--',color='blue')
ax2.plot(np.asarray(g_list)[:anzahl],np.asarray(J_h_f_list)[:anzahl],'--',color='red')
ax2.plot(np.asarray(g_list)[:anzahl],np.asarray(J_cav_f_list)[:anzahl],'--',color='orange')
ax1.plot(np.asarray(g_list)[:anzahl],np.asarray(Tot)[:anzahl],'-',color='black',linewidth=3,label=r'$\frac{\Sigma_i \langle J_{i}\rangle +\langle P \rangle}{~\gamma_h \omega_{h}}$')
ax2.plot(np.asarray(g_list)[:anzahl],np.asarray(Tot2)[:anzahl],'-',color='black',linewidth=3,label=r'$\frac{\Sigma_i \langle J_{i}\rangle +\langle P \rangle}{~\gamma_h \omega_{h}}$')
ax1.plot(g_Line,JC_g,'-',color='blue')#,label=r'$\frac{\langle J_{c}\rangle }{~\gamma_h \omega_{h}} $ Eq 31')
ax1.plot(g_Line,JH_g,'-',color='red')#,label=r'$\frac{\langle J_{h}\rangle }{~\gamma_h \omega_{h}} $ Eq 32')
ax2.plot(g_Line,JC2_g,'-',color='blue')
ax2.plot(g_Line,JH2_g,'-',color='red')
#ax2.plot(g_Line,JCAV_g,'-',color='orange')
#ax2.plot(g_Line,POWER_g,'-',color='green')
#ax2.plot(g_Line,TOT2,'-',color='black',label=r'$\frac{\Sigma_i \langle J_{i}\rangle +\langle P \rangle}{~\gamma_h \omega_{h}}$')
#ax1.plot(g_Line,TOT,'-',color='black')
#plt.plot(np.asarray(g_list)[:anzahl],np.asarray(Tot)[:anzahl],'-',color='black',alpha=0.2,label=r' $\frac{J_{tot}}{~\gamma_h \omega_{h}}$')

legend = ax2.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
Gplot.set_figheight(9)
Gplot.set_figwidth(13)

plt.show()

"""
############################################################################################
#Entropy. production
def T(omega,n):
    T=h*omega/(kb*(np.log((1/n)+1)))
    return T



nh_list=[]
Trace_list=[]
nh=0.1 #set nh again to zero
nc=0.01
nf=0.01
nh2=0.1
nh_list2=[]
Entropy=[]
Entropy2=[]
g=14*kappa
for i in range(100):
    list_temp=[]
    list_temp=Diverse_Loups.Entropy(nh2,Trans_12,a, kb,h,g,H_free,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,0)
    list_temp2=[]
    list_temp2=Diverse_Loups.Entropy(nh2,Trans_12,a, kb,h,g,H_free,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,0.5)
    #g_list.append(i/100)  #Erstellt eine Liste mit Wären von g 
    Entropy.append(list_temp)
    Entropy2.append(list_temp2)
    nh2=nh2+10
    nh_list2.append(nh2)

#Liste von Stings in floats konvertieren
#float_list2=list(np.float_(Energie_VS_g))
print(Entropy) 

#result=mesolve(H, rho0, tlist)
#print(D(c_op_list,rho)[3])


print("Die Temperatur des warmen Bades ist: ",T(omega_h,nh))
print("Die Temperatur des kalten Bades ist: ",T(omega_c,nc))
print(Trace_list_temp)

fig3, ax = plt.subplots()

ax.set_xlabel(r' $n_h$', fontsize=23)
ax.set_ylabel('Entropy production rate')
plt.title(r' Entropy Production  rate vs $n_h$ ')
plt.plot(np.asarray(nh_list2)[:100],np.asarray(Entropy)[:100,0],label=r' $\frac{J_h}{T_h}$',color='red')
plt.plot(np.asarray(nh_list2)[:100],np.asarray(Entropy)[:100,1],label=r' $\frac{J_c}{T_h}$',color='green')
plt.plot(np.asarray(nh_list2)[:100],np.asarray(Entropy)[:100,2],label=r' $\frac{J_{cav}}{T_c}$',color='pink')
#plt.plot(np.asarray(nh_list2)[:100],np.asarray(Entropy)[:100,3],label=r' $\frac{J_{cav}}{T_{cav}}$',color='orange')

legend = ax.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
fig3.set_figheight(9)
fig3.set_figwidth(13)
#Linien in plt
"""
"""plt.axvline(x=2.6)
plt.axvline(x=2.6)
plt.axvline(x=5.5)
plt.axvline(x=0.17)
plt.axvline(x=20)
plt.axvline(x=1.7)"""












Delta1=Delta2=0
gamma_h = gamma_c = 1

 

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
from Loup_for_different_coupling import Diverse_Loups as Diverse_Loups
import multiprocessing as mp
import csv
#Konstante Grössen
########################################################################################################
qutip.settings.has_mkl = False            

step=1
Delta1=Delta2=0
gamma_h = gamma_c = 1
g = 6*kappa
nc = ncav = 0.0

Delta1=0
Delta2=0
anzahl=100
nh=0
nc=nf=0

f=0.5
ladder_list=[]
ladder_list_f=[]
nh_list=[]
sigma_list=[]
sigma_list_f=[]
a_Ana=[]
a_Ana2=[]
n_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , 0 , nh, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"nh"))
n_f_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , f , nh, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"nh"))

for i in range(anzahl):
    
    ladder_list.append(Diverse_Loups.LadderOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,0))
    sigma_list.append((-Diverse_Loups.SigmaOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,0)*g-0)/kappa)
    sigma_list_f.append((-Diverse_Loups.SigmaOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,f)*g-f)/kappa)
    ladder_list_f.append(Diverse_Loups.LadderOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,f))
    a_Ana.append(((-n_f_list[6][i])*g-0.55)/kappa)
    a_Ana2.append(n_f_list[4][i])
    nh_list.append(nh)
    nh+=step
    
    


fig2, ax = plt.subplots()
ax.set_xlabel(r' $nh$', fontsize=FZ)
ax.set_ylabel(r'$<a>$',fontsize=FZ)
plt.title('')
    

#plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(a_Ana2)[:anzahl],'*',color='green',label=r'$ \langle a\rangle, (eqm), f=0.5$')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(ladder_list)[:anzahl],'-',color='black',label=r'$ a , f=0$')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(ladder_list_f)[:anzahl],'+',color='red',label=r' $a, f=0.5 $')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(sigma_list)[:anzahl],'--',color='black',label=r' $(\langle \sigma_{12}\rangle g-f)/\kappa, f=0$')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(sigma_list_f)[:anzahl],'-',color='red',alpha=0.4,label=r' $(\langle \sigma_{12}\rangle g-f)/\kappa, f=0.5$')
plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(a_Ana)[:anzahl],'-',color='green',label=r' $(\langle \sigma_{12}\rangle g-f)/\kappa, (eqm), f=0.5$')
#plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(Power_mit_f)[:anzahl],'-',color='black',alpha=0.4,label=r' $Power$')

legend = ax.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
fig3.set_figheight(9)
fig3.set_figwidth(13)





""""
#High  f approximation
f=3

nh_Lischt = np.arange(0., 100., 0.5)

# red dashes, blue squares and green triangles
plt.plot(nh_Lischt, nh_Lischt*Diverse_Loups.OccP1_Analytic(gamma_c,gamma_h,kappa,nc,nh_Lischt, ncav,f)-(nh_Lischt+1)*(1-Diverse_Loups.OccP2_Analytic(gamma_c,gamma_h,kappa,nc,nh_Lischt, ncav,f)-Diverse_Loups.OccP1_Analytic(gamma_c,gamma_h,kappa,nc,nh_Lischt, ncav,f)), 'r--')
plt.show()


"""




Delta1=Delta2=0
gamma_h = gamma_c = 1
kappa=0.07
g = 3*kappa
#g=2.8
nc = ncav = 0.0

Delta1=0
Delta2=0
anzahl=20
nh=0
nc=nf=0
step=0.03

f=0.000
ladder_list=[]
ladder_list_f=[]
f_list=[]
sigma_list=[]
sigma_list_f=[]
a_Ana=[]
a_Ana2=[]
f_list3=[]
warm=0

#n_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , f , 0, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"f"))
#n_f_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , f , warm, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"f"))

for i in range(anzahl):
   
    
    
    ladder_list.append(-Diverse_Loups.LadderOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, 0,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,f))
    #sigma_list.append((Diverse_Loups.SigmaOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, 0,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,f)*g-f)/kappa)
     
    
    #sigma_list_f.append((-Diverse_Loups.SigmaOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, warm,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,f)*g-f)/kappa)
    ladder_list_f.append(-Diverse_Loups.LadderOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, warm,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,f))
    print(ladder_list_f[i])
    #a_Ana.append(((-n_f_list[6][i])*g-f)/kappa)
    #a_Ana2.append(n_f_list[4][i])
    f_list.append(f)
    f_list3.append(f**3)
    f+=step
    
    
import pylab
import matplotlib.pyplot as plt
#g = 2.8
#nc = ncav = 0.0
#kappa=10
Delta1=0
Delta2=0
nh=0
nc=nf=0

f=0
ladder_list2=[]
#ladder_list_f=[]
f_list=[]
sigma_list=[]
sigma_list_f=[]
a_Ana=[]
a_Ana2=[]
n_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , f , 0, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"f"))
n_f_list=(Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , f , 5, ncav , nc, gamma_c, gamma_h, g ,kappa,anzahl,step,"f"))

F=np.linspace(0,step*anzahl,anzahl)
Average_A=np.imag((-1j*F*(2*F**2*np.sqrt(kappa) + g**2*np.sqrt(kappa) - g**3*np.sqrt(kappa/g**2)))/((2*F**2 + g**2)*kappa**1.5))
for i in range(anzahl):

    
    
    ladder_list2.append(Diverse_Loups.LadderOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, 0,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,f))
    sigma_list.append((-Diverse_Loups.SigmaOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, 0,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,f)*g-f)/kappa)
    
    sigma_list_f.append((-Diverse_Loups.SigmaOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, 5,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,f)*g-f)/kappa)
    
    a_Ana.append((-(n_f_list[6][i])*g-f)/kappa)
    a_Ana2.append(n_f_list[4][i])
    f_list.append(f)
    f+=step
    
    


fig2, (ax,ax1) = plt.subplots(1,2)
ax1.tick_params( labelsize=tz)
ax.tick_params( labelsize=tz)
ax1.set_xlabel(r' $\frac{f}{\gamma}$', fontsize=FZ)
ax1.set_ylabel(r'$<a>$',fontsize=FZ)
ax.set_title('a)')
ax1.set_title('b)')

ax.plot(F,Average_A,'--',color='black',label=r'$\langle a\rangle$')
#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(a_Ana2)[:anzahl],'*',color='green',label=r'$\langle a\rangle, (eqm), nh=5$')
ax.plot(np.asarray(f_list)[:anzahl],np.asarray(ladder_list2)[:anzahl],'-',color='black',label=r'$ \langle a\rangle $')
#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(ladder_list_f)[:anzahl],'-',color='red',label=r' $\langle a\rangle, n_h=5 $')
#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(sigma_list)[:anzahl],'*',color='black',label=r' $(\langle \sigma_{12}\rangle g-f)\kappa, n_h=0$')
#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(sigma_list_f)[:anzahl],'-',color='red',alpha=0.4,label=r' $(\langle \sigma_{12}\rangle g-f)/\kappa , n_h=5$')
#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(a_Ana)[:anzahl],'-',color='green',label=r' $(\langle \sigma_{12}\rangle g-f)/\kappa, eqm,  n_h=5$')
#plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(Power_mit_f)[:anzahl],'-',color='black',alpha=0.4,label=r' $Power$')

legend = ax.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')



plt.yscale('log',base=3) 
plt.xscale('log',base=3)
ax.set_xlabel(r' $\frac{f}{\gamma}$', fontsize=FZ)
ax.set_ylabel(r'$<a>$',fontsize=FZ)

    

#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(a_Ana2)[:anzahl],'*',color='green',label=r'$\langle a\rangle, (eqm), nh=5$')
#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(ladder_list)[:anzahl],'+',color='black',label=r'$ \langle a\rangle , nh=0$')
ax1.plot(np.asarray(f_list)[:anzahl],np.asarray(ladder_list)[:anzahl],'-',color='black',label=r' $\langle a\rangle$')
ax1.plot(np.asarray(f_list)[:anzahl],np.asarray(f_list3)[:anzahl],'-',color='red',label=r' $f^3/\gamma^3 $')
#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(sigma_list)[:anzahl],'*',color='black',label=r' $(\langle \sigma_{12}\rangle g-f)\kappa, n_h=0$')
#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(sigma_list_f)[:anzahl],'-',color='red',alpha=0.4,label=r' $(\langle \sigma_{12}\rangle g-f)/\kappa , n_h=5$')
#plt.plot(np.asarray(f_list)[:anzahl],np.asarray(a_Ana)[:anzahl],'-',color='green',label=r' $(\langle \sigma_{12}\rangle g-f)/\kappa, eqm  n_h=5$')
#plt.plot(np.asarray(nh_list)[:anzahl],np.asarray(Power_mit_f)[:anzahl],'-',color='black',alpha=0.4,label=r' $Power$')

legend = ax1.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
fig2.set_figheight(9)
fig2.set_figwidth(13)

plt.show()


Delta1=Delta2=0
gamma_h = gamma_c = 1

nc = ncav = 0.0
kappa = 0.2
g = 14*kappa
step=1
f=20
anzahl=200
nh=0

n_list=[]
nh_list=[]

Photonnumber_list=[]
nh2 = np.linspace(0, 70, 100)

nh3_list=[]
nh_list=[]
Trace_list=[]
nh=0 #set nh again to zero
Anal=[]
J_h_list=[]
J_h2_list=[]
J_c_list=[]
J_cav_list=[]
J_tot_list=[]
Energie_vs_nh=[]
P_ana=[]

n_list=Diverse_Loups.EquationOfMotion3(Delta1 , Delta2 , f , nh, ncav , nc, gamma_c, gamma_h, g , kappa, anzahl,step,'nh')

#n2_list=[]

for i in range(anzahl):
   
    
    #list_temp=[]
    #list_temp=Diverse_Loups.EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,1/3*kappa,omega_d,proj_2,f,omega_f,omega_2)
    nh_list.append(nh)
    #Energie_vs_nh.append(list_temp)
    #Trace_list_temp=Diverse_Loups.ProjectorP(nh,proj_1,proj_2,proj_3,Hdilde,nc,ncav,gamma_h,gamma_c,kappa,A1,A2,A3,A4,A5,A6)
    #Trace_list.append(Trace_list_temp)
    ist_temp=[]
    #list_temp=Diverse_Loups.Photonnumber(nh,a,proj_1,proj_2,proj_3,Hdilde,nc,ncav,gamma_h,gamma_c,kappa,A1,A2,A3,A4,A5,A6,omega_d,omega_f,omega_1,omega_2,H_int,f)
    #g_list.append(i/100)  #Erstellt eine Liste mit Wären von g 
    #Photonnumber_list.append(list_temp)
    
    
    J_cav_list.append(2*omega_f*kappa*(nf-n_list[0][i])/omega_h)
    J_h_list.append(omega_h*((nh*n_list[1][i]-(nh+1)*n_list[3][i]))/omega_h)
    J_c_list.append(omega_c*(nc*n_list[2][i]-(nc+1)*n_list[3][i])/omega_h)
  
    P_ana.append(-2*30*f*(n_list[4][i])/omega_h)
    J_tot_list.append(J_c_list[i]+J_h_list[i]+J_cav_list[i]+P_ana[i]/omega_h)
    
    nh3_list.append(nh)
    nh=nh+step
        
        
        
        
        
        
        
nh_Line=np.linspace(0, 200, 70)
f_Line=np.linspace(0, 1, 70)
JH=Diverse_Loups.High_f_Approx(kappa,28,gamma_c,gamma_h,nc,nh_Line,g,nf)[0]
JC=omega_c/omega_h*Diverse_Loups.High_f_Approx(kappa,28,gamma_c,gamma_h,nc,nh_Line,g,nf)[1]
JCAV=2*omega_f/omega_h*kappa*(ncav-(f/kappa)**2)+0*nh_Line
POWER=2*omega_f/omega_h*(f**2/kappa)+0*nh_Line
#JH_f=Diverse_Loups.High_f_Approx(kappa,f_Line,gamma_c,gamma_h,nc,5,g,nf)[0]
#JC_f=omega_c/omega_h*Diverse_Loups.High_f_Approx(kappa,f_Line,gamma_c,gamma_h,nc,5,g,nf)[1]
#JCAV=Diverse_Loups.High_f_Approx(kappa,f,gamma_c,gamma_h,nc,nh_Line,g,nf)[2]
#JH=(4*f**2*g**2*gammac*gammah*(-nc + nh_Line))/(gammac*gammah*kappa**2*(gammac*nc + gammah*nh_Line)*(nc + nh_Line + 3*nc*nh_Line) + 4*f**2*g**2*(gammac*(2 + 3*nc) + gammah*(2 + 3*nh_Line)))





#ax1.plot(f_Line,JC_f,color='red',label=r'$\langle J_h \rangle $ high $f$ approx')
#ax1.plot(f_Line,JH_f,color='red',label=r'$\langle J_c \rangle $ high $f$ approx')

#plt.plot(nh_Line,JCAV,color='pink',label=r'$\langle a\rangle$ "from paper"')
#plt.show()        

fig2, (ax1,ax2) = plt.subplots(1,2)
ax1.tick_params( labelsize=tz)
ax2.tick_params( labelsize=tz)

ax1.set_xlabel(r' $n_h$', fontsize=21)

#ax1.set_ylabel(r' $ \langle J \rangle $' ,fontsize=25)
#plt.title(r'$n_h$ VS Power and $J_{cav}$ ')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(P_ana)[:anzahl],'--',alpha=1,color='green')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_h_list)[:anzahl],'--',color='red')#,label=r'$\frac{\langle J_{h} \rangle }{~\gamma_h \omega_{h}}$')
#ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_tot_list)[:anzahl],'--',color='black',label=r'$\frac{\langle J_{tot} \rangle }{~\gamma_h \omega_{h}}$')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_c_list)[:anzahl],'--',color='blue')#label=r' $\frac{\langle J_{c} \rangle }{~\gamma_h \omega_{h}}$ ')
ax1.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_cav_list)[:anzahl],'--',color='orange')#label=r' $\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}}$ S')
ax1.plot(nh_Line,JCAV,color='orange',label=r'$\frac{\langle J_{cav}\rangle }{~\gamma_h \omega_{h}} $ ')
ax1.plot(nh_Line,POWER,color='green',label=r'$\frac{\langle P \rangle }{~\gamma_h \omega_{h}} $ ')
legend = ax1.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')


ax2.set_xlabel(r' $n_h$', fontsize=FZ)

#ax2.set_ylabel(r' $ \langle J \rangle $' ,fontsize=25)
#plt.title(r'$n_h$ VS $J_c$ and $J_{h}$ ')
#ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(P_ana)[:anzahl],'--',alpha=1,color='blue')
ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_h_list)[:anzahl],'--',color='red')#,label=r'$\frac{\lange J_{h}\rangle }{~\gamma_h \omega_{h}}$')
#ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_tot_list)[:anzahl],'--',color='black',label=r'$\frac{J_{tot}}{~\gamma_h \omega_{h}}$analytisch')
ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_c_list)[:anzahl],'--',color='blue')#,label=r' $\frac{\lanlge J_{c}\rangle }{~\gamma_h \omega_{h}}$ ')
#ax2.plot(np.asarray(nh3_list)[:anzahl],np.asarray(J_cav_list)[:anzahl],'--',color='orange',label=r' $\frac{J_{cav}}{~\gamma_h \omega_{h}}$ analytisch')
ax2.plot(nh_Line,JC,color='blue',label=r'$\frac{\langle J_{h}\rangle }{~\gamma_h \omega_{h}} $ ')
ax2.plot(nh_Line,JH,color='red',label=r'$\frac{\langle J_{c}\rangle }{~\gamma_h \omega_{h}} $ ')


legend = ax2.legend(loc='upper right', shadow=True, fontsize=fz)
legend.get_frame().set_facecolor('white')
fig3.set_figheight(9)
fig3.set_figwidth(13)
plt.show()


