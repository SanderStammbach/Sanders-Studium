# /bin/env/python
from calendar import c
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
from qutip import ptrace 
#from qutip import eigenstates as eigenstates


class Diverse_Loups():
    def EnergieCalculator(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list):

        

        H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

        H=H_free+H_int
        
        
        rho = steadystate(H, c_op_list) ######## Are you sure its with only photons H_free?
        rho_f=rho.ptrace(1)

        def D(c_op_list,rho):
            D=[]
            for i in range(6):
                D.append(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho-rho*c_op_list[i].dag()*c_op_list[i]))
            return D

        Liste_von_Q=[] # ExpectValue for Thermal Energy

        Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1])))
        Liste_von_Q.append(-1*np.trace(H_free*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3])))
        Liste_von_Q.append(-1*np.trace(H_free*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5])))
        #Liste_von_Q.append(g)  g in der liste anfügen

        float_list= list(np.float_(Liste_von_Q))
        print(float_list)    
        Liste_von_Q=float_list

        return(Liste_von_Q)




    def Funktion(nh,proj_1,proj_2,proj_3,H,nc,nf,gamma_h,gamma_c,kappa,A1,A2,A3,A4,A5,A6):
        nh_list=[]
        Trace_list=[]
        
        

        gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
        gamma_2=(nh)*gamma_h
        gamma_3=(nc+1)*gamma_c
        gamma_4=(nc)*gamma_c
        kappa_5=(nf+1)*kappa ####goes to zero
        kappa_6=(nf)*kappa

        c_op_list=[]    
        c_op_list.append(np.sqrt(gamma_1)*A1)
        c_op_list.append(np.sqrt(gamma_2)*A2)
        c_op_list.append(np.sqrt(gamma_3)*A3)
        c_op_list.append(np.sqrt(gamma_4)*A4)
        c_op_list.append(np.sqrt(kappa_5)*A5)
        c_op_list.append(np.sqrt(kappa_6)*A6)
            

        rho = steadystate(H, c_op_list)
        Trace_list.append(np.trace(proj_1*rho))
        Trace_list.append(np.trace(proj_2*rho))
        Trace_list.append(np.trace(proj_3*rho))
        

        float_list2= list(np.float_(Trace_list))
        #print(float_list)    
        Trace_list=float_list2

        return Trace_list
            

        #A1=Trans_13
        #A2=Trans_13.dag()
        #A3=Trans_23
        #A4=Trans_23.dag()
        #A5=a
        #A6=a.dag()
        #c_op_list=[]
        



        #gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
        #gamma_2=(nh)*gamma_h
        #gamma_3=(nc+1)*gamma_c
        #gamma_4=(nc)*gamma_c
        #gamma_5=(nf+1)*gamma_f ####goes to zero
        #gamma_6=(nf)*gamma_f


        #c_op_list.append(np.sqrt(gamma_1)*A1)
        #c_op_list.append(np.sqrt(gamma_2)*A2)
        #c_op_list.append(np.sqrt(gamma_3)*A3)
        #c_op_list.append(np.sqrt(gamma_4)*A4)
        #c_op_list.append(np.sqrt(gamma_5)*A5)
        #c_op_list.append(np.sqrt(gamma_6)*A6)

       
    def Entropy(nh,Trans_12,a, kb,h,g,H,H_free,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_c,omega_h,omega_f):
        A1=Trans_13
        A2=Trans_13.dag()
        A3=Trans_23
        A4=Trans_23.dag()
        A5=a
        A6=a.dag()
        
        H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

        H=H_free+H_int
        
        
        
        
        gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
        gamma_2=(nh)*gamma_h
        gamma_3=(nc+1)*gamma_c
        gamma_4=(nc)*gamma_c
        kappa_5=(nf+1)*kappa ####goes to zero
        kappa_6=(nf)*kappa
        
        c_op_list=[]    
        c_op_list.append(np.sqrt(gamma_1)*A1)
        c_op_list.append(np.sqrt(gamma_2)*A2)
        c_op_list.append(np.sqrt(gamma_3)*A3)
        c_op_list.append(np.sqrt(gamma_4)*A4)
        c_op_list.append(np.sqrt(kappa_5)*A5)
        c_op_list.append(np.sqrt(kappa_6)*A6)
        def T(n,omega):
            
            T=h*omega/(kb*(np.log((1/n)+1)))
            return T

        rho = steadystate(H, c_op_list) ######## Are you sure its with only photons H_free?
        rho_f=rho.ptrace(1)
        def D(c_op_list,rho):
            D=[]
            for i in range(6):
                D.append(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho+rho*c_op_list[i].dag()*c_op_list[i]))
            return D

        Liste_von_Q=[] # ExpectValue for Thermal Energy
        Liste_von_Q_f=[]
        Liste_von_Q_c=[]
        Liste_von_Q_h=[]
        #Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1]))/(T(nh,omega_h))+np.trace(H_free*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3]))/T(nc,omega_c)+np.trace(H_free*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5]))/T(nf,omega_f))
        Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1]))+np.trace(H_free*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3]))+np.trace(H_free*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5])))
        Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1])))#/(T(nh,omega_h)))#liste_von_Q_h
        Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3])))#/(T(nc,omega_c)))
        Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5])))#/(T(nf,omega_f)))
            #Liste_von_Q.append(g)  g in der liste anfügen

        float_list= list(np.float_(Liste_von_Q))
        print(float_list)    
        Liste_von_Entropy=float_list

        return(Liste_von_Entropy)

    def Photonnumber(nh2,a,proj_1,proj_2,proj_3,H,nc,nf,gamma_h,gamma_c,kappa,A1,A2,A3,A4,A5,A6):
        
        Trace_list=[]
        
        

        gamma_1=(nh2+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
        gamma_2=(nh2)*gamma_h
        gamma_3=(nc+1)*gamma_c
        gamma_4=(nc)*gamma_c
        kappa_5=(nf+1)*kappa ####goes to zero
        kappa_6=(nf)*kappa

        c_op_list=[]    
        c_op_list.append(np.sqrt(gamma_1)*A1)
        c_op_list.append(np.sqrt(gamma_2)*A2)
        c_op_list.append(np.sqrt(gamma_3)*A3)
        c_op_list.append(np.sqrt(gamma_4)*A4)
        c_op_list.append(np.sqrt(kappa_5)*A5)
        c_op_list.append(np.sqrt(kappa_6)*A6)
            

        rho = steadystate(H, c_op_list)
        Trace_list.append

        n=np.trace(a.dag()*a*rho)
        float_list2= list(np.float_(Trace_list))
        #print(float_list)    
        Trace_list=float_list2



        return n 
    
    def EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2):

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

        omega_1=0
        V=f*a.dag()+f*a
        H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())
        
        H=H_free+H_int -omega_d*(a.dag()*a+proj_2) #+ f*(a+a.dag()) 
        H_free1=H_free 
        Hdilde=H_int+V +(omega_2-(30+omega_d))*(a.dag()*a)+(omega_f-omega_d)*proj_2   
        rho = steadystate(Hdilde, c_op_list) ######## Are you sure its with only photons H_free?
        rho_f=rho.ptrace(1)

        def D(c_op_list,rho):
            D=[]
            for i in range(6):
                D.append(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho+rho*c_op_list[i].dag()*c_op_list[i]))
            return D

        Liste_von_Q=[] # ExpectValue for Thermal Energy

        Liste_von_Q.append(np.trace(H_free1*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1])))
        Liste_von_Q.append(np.trace(H_free1*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3])))
        Liste_von_Q.append(np.trace(H_free1*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5])))
        #Liste_von_Q.append(g)  g in der liste anfügen

        float_list= list(np.float_(Liste_von_Q))
        print(float_list)    
        Liste_von_Q=float_list

        return(Liste_von_Q)
    





    def P(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,f,omega_2):
      
        def Ptr(H_free,Hdilde,rho):
                Power=0
                Power=-1j*np.trace(H_free*Hdilde*rho-H_free*rho*Hdilde)
                return Power
                #dt_rho=dt_rho+(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho-rho*c_op_list[i].dag()*c_op_list[i]))
        P_list=[]
        g=0
        for i in range(200):
            g=g+i/100
            H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

            H=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 
        
            V=f*a.dag()+f*a
            Hdilde=H_int+V +(30-(30+omega_d))*(a.dag()*a)+(omega_f-omega_d)*proj_2  
            rho = steadystate(Hdilde, c_op_list) ######## Are you sure its with only photons H_free?
        
            P_list.append(Ptr(H_free,Hdilde,rho))

        return(P_list)

    
    def P2(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,f,omega_2,g):

          
            H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

            H=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 
        
            V=f*a.dag()+f*a
            Hdilde=H_int+V +(30-(30+omega_d))*(a.dag()*a)+(omega_f-omega_d)*proj_2  
            rho = steadystate(Hdilde, c_op_list) ######## Are you sure its with only photons H_free?
            Power=-1j*np.trace((H_free*Hdilde-Hdilde*H_free)*rho)
            

            return(Power)
    

    def Fullcounting(c_op_list):
        x=c_op_list.qutip.eigenstates()

        return x

    def current(H_free, Trans_12, a, h,c_op_list,omega_d,omega_f ,proj_2,g,f,anzahl):
      
        def D(c_op_list,rho):
            D=[]
            for i in range(6):
                D.append(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho+rho*c_op_list[i].dag()*c_op_list[i]))
            return D
        
        
        float_list=[]
        for i in range(anzahl):
            f=f+1/80
            H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

            H=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 
        
            V=f*a.dag()+f*a
            Hdilde=H_int+V +(30-(30+omega_d))*(a.dag()*a)+(omega_f-omega_d)*proj_2  
            rho = steadystate(Hdilde, c_op_list)
            Liste_von_Q=[]
            Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1])))
            Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3])))
            Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5])))
            float_list.append(list(np.float_(Liste_von_Q)))
            print(Liste_von_Q)
        
        #Liste_von_Q.append(g)  g in der liste anfügen

        
           
        Liste_von_Q=float_list

        return(Liste_von_Q)
    


    def P3(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,omega_2,g,f,anzahl):
      
        def Ptr(H_free,Hdilde,rho):
                Power=0
                Power=-1j*np.trace(H_free*Hdilde*rho-H_free*rho*Hdilde)
                return Power
                #dt_rho=dt_rho+(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho-rho*c_op_list[i].dag()*c_op_list[i]))
        P_list=[]
        
        for i in range(anzahl):
            f=f+1/80
            H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

            H=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 
        
            
            Hdilde=H_int+f*(a+a.dag()) +30*(a.dag()*a)#+(omega_f-omega_d)*proj_2  
            rho = steadystate(Hdilde, c_op_list) 
        
            P_list.append(Ptr(H_free,Hdilde,rho))

        return(P_list)
    


    
    def P4(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,omega_2,g,f,anzahl):
      
        def Ptr(H_free,Hdilde,rho):
                Power=0
                Power=-1j*np.trace(H_free*Hdilde*rho-H_free*rho*Hdilde)
                return Power
                #dt_rho=dt_rho+(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho-rho*c_op_list[i].dag()*c_op_list[i]))
        P_list=[]
        
        for i in range(anzahl):
            g=g+1/120
            H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

            H=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 
        
            
            Hdilde=H_int+f*(a+a.dag()) +30*(a.dag()*a)#+(omega_f-omega_d)*proj_2  
            rho = steadystate(Hdilde, c_op_list) 
        
            P_list.append(Ptr(H_free,Hdilde,rho))

        return(P_list)
    
    
    
