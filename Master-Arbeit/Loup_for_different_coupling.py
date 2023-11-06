# /bin/env/python
from ast import mod
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
from math import gamma, gcd
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
from sympy import symbols, Eq, solve,nsolve
from scipy.optimize import fsolve
from math import exp




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
    



    def EnergieCalculator_mit_faktor(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2,omega_h):
        
        gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
        gamma_2=(nh)*gamma_h
        gamma_3=(nc+1)*gamma_c
        gamma_4=(nc)*gamma_c
        kappa_5=(nf+1)*2*kappa ####goes to zero
        kappa_6=(nf)*2*kappa/2

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
    


    def ProjectorP(nh,proj_1,proj_2,proj_3,Hdilde,nc,nf,gamma_h,gamma_c,kappa,A1,A2,A3,A4,A5,A6):
        nh_list=[]
        Trace_list=[]
        
        

        gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
        gamma_2=(nh)*gamma_h
        gamma_3=(nc+1)*gamma_c
        gamma_4=(nc)*gamma_c
        kappa_5=(nf+1)*2*kappa ####goes to zero
        kappa_6=(nf)*2*kappa

        c_op_list=[]    
        c_op_list.append(np.sqrt(gamma_1)*A1)
        c_op_list.append(np.sqrt(gamma_2)*A2)
        c_op_list.append(np.sqrt(gamma_3)*A3)
        c_op_list.append(np.sqrt(gamma_4)*A4)
        c_op_list.append(np.sqrt(kappa_5)*A5)
        c_op_list.append(np.sqrt(kappa_6)*A6)
            

        rho = steadystate(Hdilde, c_op_list)
        Trace_list.append(np.trace(proj_1*rho))
        Trace_list.append(np.trace(proj_2*rho))
        Trace_list.append(np.trace(proj_3*rho))
        

        float_list2= list(np.float_(Trace_list))
        #print(float_list)    
        Trace_list=float_list2

        return Trace_list
    

    
    
    def ProjectorP_mit_faktor(Trans_12,h,g,nh,proj_1,proj_2,proj_3,V ,omega_2,omega_1,omega_d,omega_f,a ,nc,nf,gamma_h,gamma_c,kappa,A1,A2,A3,A4,A5,A6):
        
        H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())
        Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a) 
        nh_list=[]
        Trace_list=[]
        
        

        gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
        gamma_2=(nh)*gamma_h
        gamma_3=(nc+1)*gamma_c
        gamma_4=(nc)*gamma_c
        kappa_5=(nf+1)*2*kappa ####goes to zero
        kappa_6=(nf)*2*kappa

        c_op_list=[]    
        c_op_list.append(np.sqrt(gamma_1)*A1)
        c_op_list.append(np.sqrt(gamma_2)*A2)
        c_op_list.append(np.sqrt(gamma_3)*A3)
        c_op_list.append(np.sqrt(gamma_4)*A4)
        c_op_list.append(np.sqrt(kappa_5)*A5)
        c_op_list.append(np.sqrt(kappa_6)*A6)
            

        rho = steadystate(Hdilde, c_op_list)
        Trace_list.append(np.trace(proj_1*rho))
        Trace_list.append(np.trace(proj_2*rho))
        Trace_list.append(np.trace(proj_3*rho))
        

        float_list2= list(np.float_(Trace_list))
        #print(float_list)    
        Trace_list=float_list2

        return Trace_list
            

    def Photonnumber(nh,a,proj_1,proj_2,proj_3,Hdilde,nc,ncav,gamma_h,gamma_c,kappa,A1,A2,A3,A4,A5,A6,omega_d,omega_f,omega_1,omega_2,H_int,f):
        
        Trace_list=[]
        
        

        gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
        gamma_2=(nh)*gamma_h
        gamma_3=(nc+1)*gamma_c
        gamma_4=(nc)*gamma_c
        kappa_5=(ncav+1)*2*kappa ####goes to zero
        kappa_6=(ncav)*2*kappa

        c_op_list=[]    
        c_op_list.append(np.sqrt(gamma_1)*A1)
        c_op_list.append(np.sqrt(gamma_2)*A2)
        c_op_list.append(np.sqrt(gamma_3)*A3)
        c_op_list.append(np.sqrt(gamma_4)*A4)
        c_op_list.append(np.sqrt(kappa_5)*A5)
        c_op_list.append(np.sqrt(kappa_6)*A6)

        
        V=f*(a+a.dag())  
        Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)   
        rho = steadystate(Hdilde, c_op_list)
        Trace_list.append

        n=np.trace(a.dag()*a*rho)
        float_list2= list(np.float_(Trace_list))
        #print(float_list)    
        Trace_list=float_list2



        return n 


    def Photonnumber2(nh,a,proj_1,proj_2,proj_3,Trans_12,nc,ncav,gamma_h,gamma_c,kappa,A1,A2,A3,A4,A5,A6,omega_d,omega_f,omega_1,omega_2,f,g):
        
        Trace_list=[]
        
        

        gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
        gamma_2=(nh)*gamma_h
        gamma_3=(nc+1)*gamma_c
        gamma_4=(nc)*gamma_c
        kappa_5=(ncav+1)*2*kappa ####goes to zero
        kappa_6=(ncav)*2*kappa

        c_op_list=[]    
        c_op_list.append(np.sqrt(gamma_1)*A1)
        c_op_list.append(np.sqrt(gamma_2)*A2)
        c_op_list.append(np.sqrt(gamma_3)*A3)
        c_op_list.append(np.sqrt(gamma_4)*A4)
        c_op_list.append(np.sqrt(kappa_5)*A5)
        c_op_list.append(np.sqrt(kappa_6)*A6)

        H_int=g*(Trans_12*a.dag()+a*Trans_12.dag())
        V=f*(a+a.dag())  
        Hdilde=H_int+V +(0)*(proj_2)+(0)*(a.dag()*a)
        rho = steadystate(Hdilde, c_op_list)
        Trace_list.append

        n=np.trace(a.dag()*a*rho)
        float_list2= list(np.float_(Trace_list))
        #print(float_list)    
        Trace_list=float_list2



        return n 

       
    def Entropy(nh,Trans_12,a, kb,h,g,H_free,nc,nf,gamma_h,gamma_c,kappa,Trans_13,Trans_23,omega_f,omega_d,omega_1,omega_2,proj_2,f):
        
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
     
        Liste_von_Q.append((np.trace(-H_free1*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1])))/(T(nh,omega_h)*omega_h))
        Liste_von_Q.append((np.trace(-H_free1*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3])))/(T(nc,omega_c)*omega_h))
        Liste_von_Q.append((np.trace(-H_free1*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5])))/(T(nf,omega_f)*omega_h))
        Liste_von_Q.append((np.trace(-H_free1*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1])))/(T(nh,omega_h)*omega_h)+(np.trace(-H_free1*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3])))/(T(nc,omega_c)*omega_h)+(np.trace(-H_free1*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5])))/(T(nf,omega_f)*omega_h))
        #Liste_von_Q.append(g)  g in der liste anfügen

        float_list= list(np.float_(Liste_von_Q))
            
        Liste_von_Q=float_list

        return(Liste_von_Q)





    def P(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,f,anzahl,step,omega_h):
        
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
        omega_1=0
        omega_2=30
        



        def Ptr(H_free,Hdilde,rho):
                Power=0
                Power=-1j*np.trace(H_free*Hdilde*rho-H_free*rho*Hdilde)
                return Power
                #dt_rho=dt_rho+(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho-rho*c_op_list[i].dag()*c_op_list[i]))
        P_list=[]
        g=0
        for i in range(anzahl):
            g=g+step
            H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

            H=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 
        
            V=f*a.dag()+f*a
            Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)    
            rho = steadystate(Hdilde, c_op_list) ######## Are you sure its with only photons H_free?
        
            P_list.append(Ptr(H_free,Hdilde,rho)/omega_h)

        return(P_list)

    
    def P2(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,f,omega_2,g):
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
          
            H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

            H=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 
            omega_1=0
            V=f*a.dag()+f*a
            Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)   
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
            Hdilde=H_int+V +0*(a.dag()*a)+(omega_f-omega_d)*proj_2  
            rho = steadystate(Hdilde, c_op_list)
            Liste_von_Q=[]
            Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[0]+D(c_op_list,rho)[1])))
            Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[2]+D(c_op_list,rho)[3])))
            Liste_von_Q.append(np.trace(H_free*(D(c_op_list,rho)[4]+D(c_op_list,rho)[5])))
            float_list.append(list(np.float_(Liste_von_Q)))
            
        
        #Liste_von_Q.append(g)  g in der liste anfügen

        
           
        Liste_von_Q=float_list

        return(Liste_von_Q)
    


    def P3(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,omega_2,g,f,anzahl,step):
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

        def Ptr(H_free,Hdilde,rho):
                Power=0
                Power=-1j*np.trace(H_free*Hdilde*rho-H_free*rho*Hdilde)
                return Power
                #dt_rho=dt_rho+(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho-rho*c_op_list[i].dag()*c_op_list[i]))
        P_list=[]
        
        for i in range(anzahl):
            f=f+step
            H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

            H=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 
        
            V=f*(a+a.dag())
            omega_1=0
            Hdilde=H_int+V +(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a)   
            rho = steadystate(Hdilde, c_op_list) 
        
            P_list.append(Ptr(H_free,Hdilde,rho))

        return(P_list)
    


    
    def P4(H_free, Trans_12, Trans_13, Trans_23, a, nh,nf,nc,h,kb,gamma_h,gamma_c,kappa,c_op_list,omega_d,omega_f ,proj_2,omega_2,g,f,anzahl):
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
        def Ptr(H_free,Hdilde,rho):
                Power=0
                Power=-1j*np.trace(H_free*Hdilde*rho-H_free*rho*Hdilde)
                return Power
                #dt_rho=dt_rho+(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho-rho*c_op_list[i].dag()*c_op_list[i]))
        P_list=[]
        
        for i in range(anzahl):
            g=g+1/12
            H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

            H=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 
        
            
            Hdilde=H_int+f*(a+a.dag()) +30*(a.dag()*a)+(omega_f-omega_d)*proj_2  
            rho = steadystate(Hdilde, c_op_list) 
        
            P_list.append(Ptr(H_free,Hdilde,rho))

        return(P_list)
    


    def P5(g,H_free, Trans_12, Trans_13, Trans_23, a,nf,nc,h,kb,gamma_h,gamma_c,kappa,omega_d,proj_2,f,omega_f,omega_2, anzahl,step):
        

        
        def Ptr(H_free,Hdilde,rho):
                Power=0
                Power=-1j*np.trace(H_free*Hdilde*rho-H_free*rho*Hdilde)
                return Power
                #dt_rho=dt_rho+(c_op_list[i]*rho*c_op_list[i].dag()-1/2*(c_op_list[i].dag()*c_op_list[i]*rho-rho*c_op_list[i].dag()*c_op_list[i]))
        P_list=[]


        nh=0
        nh_List=[]
        for i in range(anzahl):

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
            nh=nh+step
            nh_List.append(nh)
            H_int=h*g*(Trans_12*a.dag()+a*Trans_12.dag())

            H=H_free+H_int -omega_d*(a.dag()*a+proj_2) + f*(a+a.dag()) 
        
            
            Hdilde=H_int+f*(a+a.dag()) +0*(a.dag()*a)+(omega_f-omega_d)*proj_2  
            rho = steadystate(Hdilde, c_op_list) 
        
            P_list.append(Ptr(H_free,Hdilde,rho))


        return(P_list,nh_List)

    def N_Analytic(gammac,gammah,kappa,g,nh,f,nc):#The new one
        return (        (4*f**2*g**2*(2*gammac + gammah*(2 + 3*nh)) - kappa*(gammac*gammah*kappa*nh*(2*kappa + gammah*nh) + 2*g**2*(-(gammac*gammah*nh) + 2*kappa*(gammac + gammah + 2*gammah*nh))) + np.sqrt(64*g**2*kappa**2*(2*gammac + gammah*(2 + 3*nh))*(2*g**2*gammac*gammah*kappa*nh + f**2*(gammac*gammah*nh*(2*kappa + gammah*nh) + 4*g**2*(gammac + gammah + 2*gammah*nh))) + 4*(-4*f**2*g**2*(2*gammac + gammah*(2 + 3*nh)) + kappa*(gammac*gammah*kappa*nh*(2*kappa + gammah*nh) + g**2*(-2*gammac*gammah*nh + 4*kappa*(gammac + gammah + 2*gammah*nh))))**2)/2.)/(8.*g**2*kappa**2*(2*gammac + gammah*(2 + 3*nh))))

        
    def N_Analytic2(gamma,kappa,g,nh,ncav,nc):
          
        return (-((2*gamma + gamma*nc + 2*gamma*nh)*(gamma**2*kappa*(2*kappa + gamma*nc + gamma*nh)*(nc + nh + 3*nc*nh) + 4*g**2*(-2*gamma*kappa*(-1 + 2*ncav) + gamma*nc*(gamma + kappa - 3*kappa*ncav) - gamma*(gamma - 2*kappa + 3*kappa*ncav)*nh))) + np.sqrt((2*gamma + gamma*nc + 2*gamma*nh)**2*(gamma**4*kappa**2*(2*kappa + gamma*nc + gamma*nh)**2*(nc + nh + 3*nc*nh)**2 + 8*g**2*gamma**2*kappa*(2*kappa + gamma*nc + gamma*nh)*(nc + nh + 3*nc*nh)*(2*gamma*kappa*(1 + 2*ncav) + gamma*nc*(gamma + kappa + 3*kappa*ncav) + gamma*(-gamma + 2*kappa + 3*kappa*ncav)*nh) + 16*g**4*(gamma**4*(nc - nh)**2 + kappa**2*(gamma + gamma*(1 + nc + 2*ncav + 3*nc*ncav) + 2*gamma*nh + gamma*ncav*(2 + 3*nh))**2 +  2*gamma**2*kappa*(gamma*nc**2*(1 - 3*ncav + 6*nh) + nh*(2*gamma*(3 + 2*ncav) + gamma*(4 + 3*ncav)*nh) + nc*(2*gamma - 2*gamma*ncav + 3*gamma*(3 + ncav)*nh + 6*gamma*nh*(1 + nh) - gamma*ncav*(2 + 3*nh)))))))/(8.*g**2*kappa*(2*gamma + gamma*nc + 2*gamma*nh)*(gamma*(2 + 3*nc) + gamma*(2 + 3*nh)))
    

    def OccP1_Analytic(gammac,gammah,kappa,nc,nh, ncav,f):

        return(-4*f**2*(gammac + gammah + 2*gammac*nc + gammah*nh) + kappa*(gammac*gammah*nc*(1 + nh) - 2*kappa*ncav*(gammac + gammah + 2*gammac*nc + gammah*nh)))/(gammac*gammah*kappa*(nc + nh + 3*nc*nh))

    def OccP2_Analytic(gammac,gammah,kappa,nc,nh,ncav ,f):
        return((4*f**2*(gammac + gammah + gammac*nc + 2*gammah*nh) + kappa*(gammac*gammah*(1 + nc)*nh + 2*kappa*ncav*(gammac + gammah + gammac*nc + 2*gammah*nh)))/(gammac*gammah*kappa*(nc + nh + 3*nc*nh)))

    def EquationOfMotion(Delta1 , Delta2 , f , nh, ncav , nc, gammac, gammah, g , kappa):
         
        import sympy as sym
        ai, ar, Trans12R,Trans12I,ATrans12I,ATrans12R,p1,p2,n = sym.symbols('ai ar Trans12R Trans12I ATrans12I ATrans12R p1 p2 n')
        
        eq1=sym.Eq(-(ar*Delta1) - f - ai*kappa - g*Trans12R,0)
        eq2=sym.Eq(ai*Delta1 - ar*kappa + g*Trans12I,0)
        eq3=sym.Eq(-2*ATrans12I*g + gammac*((1 + nc)*(1 - p1 - p2) - nc*p2),0)
        eq4=sym.Eq(2*ATrans12I*g + gammah*(-(nh*p1) + (1 + nh)*(1 - p1 - p2)),0)
        eq5=sym.Eq(ATrans12R*Delta1 - ATrans12R*Delta2 - ATrans12I*kappa - (ATrans12I*gammac*nc)/2 -(ATrans12I*gammah*nh)/2 + g*(p2 + n*(-p1 + p2)) - f*Trans12R,0)
        eq6=sym.Eq(-(ATrans12I*Delta1) + ATrans12I*Delta2 - ATrans12R*kappa - (ATrans12R*gammac*nc)/2- (ATrans12R*gammah*nh)/2 + f*Trans12I,0)
        eq7=sym.Eq(-(ar*g*(p1 - p2)) - (gammac*nc*Trans12I)/2. - (gammah*nh*Trans12I)/2,0)
        eq8=sym.Eq(ai*g*(p1 - p2) - (gammac*nc*Trans12R)/2. - (gammah*nh*Trans12R)/2,0)     
        eq9=sym.Eq(-2*ai*f + 2*ATrans12I*g + kappa*(-n + ncav),0)
        

        sol = sym.nsolve((eq1,eq2,eq3,eq4,eq5,eq6,eq6,eq8,eq9),(ai, ar, Trans12R, Trans12I, ATrans12I, ATrans12R, p1, p2 ,n))


        return sol


    """from scipy.optimize import fsolve
    def EqO(Delta1 , Delta2 , f , nh, ncav , nc, gammac, gammah, g , kappa,ai, ar, Trans12R,Trans12I,ATrans12I,ATrans12R,p1,p2,n):
         def eqs(Delta1 , Delta2 , f , nh, ncav , nc, gammac, gammah, g , kappa,ai, ar, Trans12R,Trans12I,ATrans12I,ATrans12R,p1,p2,n):
              return(-(ar*Delta1) - f - ai*kappa - g*Trans12R,ai*Delta1 - ar*kappa + g*Trans12I,-2*ATrans12I*g + gammac*((1 + nc)*(1 - p1 - p2) - nc*p2),2*ATrans12I*g + gammah*(-(nh*p1) + (1 + nh)*(1 - p1 - p2)),ATrans12R*Delta1 - ATrans12R*Delta2 - ATrans12I*kappa - (ATrans12I*gammac*nc)/2 -(ATrans12I*gammah*nh)/2 + g*(p2 + n*(-p1 + p2)) - f*Trans12R,-(ATrans12I*Delta1) + ATrans12I*Delta2 - ATrans12R*kappa - (ATrans12R*gammac*nc)/2- (ATrans12R*gammah*nh)/2 + f*Trans12I,-(ar*g*(p1 - p2)) - (gammac*nc*Trans12I)/2. - (gammah*nh*Trans12I)/2,ai*g*(p1 - p2) - (gammac*nc*Trans12R)/2. - (gammah*nh*Trans12R)/2,-2*ai*f + 2*ATrans12I*g + kappa*(-n + ncav))
         

         losung=fsolve(eqs,[0,0,0,0,0,0,0,0,0])

         return losung"""
    



 



    def EquationOfMotion2(Delta1 , Delta2 , f , nh, ncav , nc, gammac, gammah, g , kappa):
        
        def equations(vars):
            ai ,ar, Trans12R,Trans12I,ATrans12I,ATrans12R,p1,p2,n=vars
        
         
        
            eq1=-(ar*Delta1) -f- ai*kappa - g*Trans12R
            eq2=(ai*Delta1 - ar*kappa + g*Trans12I)
            eq3=-2*ATrans12I*g + gammac*((1 + nc)*(1 - p1 - p2) - nc*p2)
            eq4=(2*ATrans12I*g + gammah*(-(nh*p1) + (1 + nh)*(1 - p1 - p2)))
            eq5=(ATrans12R*Delta1 - ATrans12R*Delta2 - ATrans12I*kappa - (ATrans12I*gammac*nc)/2 -(ATrans12I*gammah*nh)/2 + g*(p2 + n*(-p1 + p2)) - f*Trans12R)
            eq6=(-(ATrans12I*Delta1) + ATrans12I*Delta2 - ATrans12R*kappa - (ATrans12R*gammac*nc)/2- (ATrans12R*gammah*nh)/2 + f*Trans12I)
            eq7=(-(ar*g*(p1 - p2)) - (gammac*nc*Trans12I)/2. - (gammah*nh*Trans12I)/2)
            eq8=(ai*g*(p1 - p2) - (gammac*nc*Trans12R)/2. - (gammah*nh*Trans12R)/2)     
            eq9=(+2*ai*f + 2*ATrans12I*g + 2*kappa*(-n + ncav))

            return[eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9]
        
        ai,ar, Trans12R,Trans12I,ATrans12I,ATrans12R,p1,p2,n=fsolve(equations,(-1,0,0,0,0,0,1,1,1))

        
        return [n,p1,p2,(1-p1-p2),ai]
    

    def EquationOfMotion3(Delta1 , Delta2 , f , nh, ncav , nc, gammac, gammah, g , kappa, anzahl,step,variable):
        if variable=="nh":
            g=g
            nh=0
            f=f

        elif variable=="f":
            f=0
            g=g
            nh=nh

        else:
            g=0
            nh=nh
            f=f

        n_list=[]
        p1_list=[]
        p2_list=[]
        ai_list=[]
        p3_list=[]
        ar_list=[]
        sigma_12_list=[]
        nS=0.1
        p1S=0.8
        p2S=0.07
        aiS=-1
        arS=0
        Trans12RS=0
        Trans12IS=0
        ATrans12IS=0.2
        ATrans12RS=0
        
        
        
        
        for i in range(anzahl):
            def equations(vars):
                ai ,ar, Trans12R,Trans12I,ATrans12I,ATrans12R,p1,p2,n=vars
        
         
        
                eq1=-(ar*Delta1) -f- ai*kappa - g*Trans12R
                eq2=(ai*Delta1 - ar*kappa + g*Trans12I)
                eq3=-2*ATrans12I*g + gammac*((1 + nc)*(1 - p1 - p2) - nc*p2)
                eq4=(2*ATrans12I*g + gammah*(-(nh*p1) + (1 + nh)*(1 - p1 - p2)))
                eq5=(ATrans12R*Delta1 - ATrans12R*Delta2 - ATrans12I*kappa - (ATrans12I*gammac*nc)/2 -(ATrans12I*gammah*nh)/2 + g*(p2 + n*(-p1 + p2)) - f*Trans12R)
                eq6=(-(ATrans12I*Delta1) + ATrans12I*Delta2 - ATrans12R*kappa - (ATrans12R*gammac*nc)/2- (ATrans12R*gammah*nh)/2 + f*Trans12I)
                eq7=(-(ar*g*(p1 - p2)) - (gammac*nc*Trans12I)/2. - (gammah*nh*Trans12I)/2)
                eq8=(ai*g*(p1 - p2) - (gammac*nc*Trans12R)/2 - (gammah*nh*Trans12R)/2)     
                eq9=(-2*ai*f + 2*ATrans12I*g + 2*kappa*(-n + ncav))

                return[eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9]
        
            ai,ar, Trans12R,Trans12I,ATrans12I,ATrans12R,p1,p2,n=fsolve(equations,(aiS,arS,Trans12RS,Trans12IS,ATrans12IS,ATrans12RS,p1S,p2S,nS))
            n_list.append(n)
            p1_list.append(p1)
            p2_list.append(p2)
            ai_list.append(ai)
            ar_list.append(ar)
            sigma_12_list.append(Trans12R)
            p3=1-(p1+p2)
            p3_list.append(p3)
            nS=n
            print('fsddsffsdasdfasdffddasf',p2S)
            p1S=p1
            p2S=p2
            aiS=ai
            arS=ar
            Trans12RS=Trans12R
            Trans12IS=Trans12I
            ATrans12IS=ATrans12I
            ATrans12RS=ATrans12R
            if variable=="g":
                g+=step
            elif variable=="nh":
                nh+=step

            elif variable=="f":
                f+=step
            else:
                print('error')
        return n_list,p1_list,p2_list, p3_list,ai_list,ar_list,sigma_12_list
    


    def LadderOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,f):
        
        Trace_list=[]
        

        

        gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
        gamma_2=(nh)*gamma_h
        gamma_3=(nc+1)*gamma_c
        gamma_4=(nc)*gamma_c
        kappa_5=(ncav+1)*2*kappa ####goes to zero
        kappa_6=(ncav)*2*kappa

        A1=Trans_13
        A2=Trans_13.dag()
        A3=Trans_23
        A4=Trans_23.dag()
        A5=a
        A6=a.dag()

        c_op_list=[]    
        c_op_list.append(np.sqrt(gamma_1)*A1)
        c_op_list.append(np.sqrt(gamma_2)*A2)
        c_op_list.append(np.sqrt(gamma_3)*A3)
        c_op_list.append(np.sqrt(gamma_4)*A4)
        c_op_list.append(np.sqrt(kappa_5)*A5)
        c_op_list.append(np.sqrt(kappa_6)*A6)

        H_int=g*(Trans_12*a.dag()+a*Trans_12.dag())
        V=f*(a+a.dag())  
        Hdilde=H_int+f*(a+a.dag()) #+(omega_2-(omega_1+omega_d))*(proj_2)+(omega_f-omega_d)*(a.dag()*a) 
        rho = steadystate(Hdilde, c_op_list)
        
        rho_f=rho.ptrace(1)  ### State in the cavity

        

        A=np.imag(np.trace(a*rho))
        



        return A 
        

        

    def SigmaOperator(g,H_free, Trans_12, Trans_13, Trans_23, a, nh,ncav,nc,h,kb,gamma_h,gamma_c,kappa,proj_2,f):
        
            Trace_list=[]
        
        

            gamma_1=(nh+1)*gamma_h #### unsicher wegen vorfaktor 1/2 
            gamma_2=(nh)*gamma_h
            gamma_3=(nc+1)*gamma_c
            gamma_4=(nc)*gamma_c
            kappa_5=(ncav+1)*2*kappa ####goes to zero
            kappa_6=(ncav)*2*kappa

            A1=Trans_13
            A2=Trans_13.dag()
            A3=Trans_23
            A4=Trans_23.dag()
            A5=a
            A6=a.dag()

            c_op_list=[]    
            c_op_list.append(np.sqrt(gamma_1)*A1)
            c_op_list.append(np.sqrt(gamma_2)*A2)
            c_op_list.append(np.sqrt(gamma_3)*A3)
            c_op_list.append(np.sqrt(gamma_4)*A4)
            c_op_list.append(np.sqrt(kappa_5)*A5)
            c_op_list.append(np.sqrt(kappa_6)*A6)

            H_int=g*(Trans_12*a.dag()+a*Trans_12.dag())
            
            Hdilde=H_int+f*(a+a.dag())   +(0)*(proj_2)+(0)*(a.dag()*a)
            rho = steadystate(Hdilde, c_op_list)
            Trace_list.append
            rho_f=rho.ptrace(1)  ### State in the cavity

        

            sigma=np.real(np.trace(Trans_12*rho))
        



            return sigma 