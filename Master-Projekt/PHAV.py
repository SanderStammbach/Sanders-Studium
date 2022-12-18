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
from qutip import *
from qutip import ptrace 
from Loup_for_different_coupling import Diverse_Loups as Diverse_Loups
import multiprocessing as mp
import csv

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
norm=8200
#a=4 bestes ergebnis
a=3
nph=30
b_fock_alt=qutip.states.fock(nph,0)

rho_alt=b_fock_alt*b_fock_alt.dag()

rho_alt_2=rho_alt






y=0

for n in range(nph):
    fock=qutip.states.fock(nph,n)
    rho_PHAV=((np.exp(-1*(a**2))+(a**(2*n))/(np.math.factorial(n)))/8200)*(fock*fock.dag())
   
    x=(np.exp(-a*a)+(a**(2*n))/(np.math.factorial(n)))/8200
    print(x)
    y=x+y
    rho_alt=rho_PHAV+rho_alt


RHO=rho_alt-rho_alt_2
qutip.plot_wigner_fock_distribution(RHO,colorbar='colorbar')
print(RHO,y)
plt.show()