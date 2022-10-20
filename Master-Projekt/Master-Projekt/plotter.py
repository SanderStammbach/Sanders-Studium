from lzma import FILTER_ARM, FILTER_SPARC
from unittest import skip
from numpy import arange
from pandas import read_csv
from scipy.optimize import curve_fit
from matplotlib import pyplot
import numpy as np
import matplotlib.pyplot as plt

from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
# define the true objective function

#def VanDerWallsGL(v,a1,b1,r):
#    return r*50/(v-b1)-a1/(v**2)
    

data1 = np.loadtxt('Speicherort.csv',delimiter=',',dtype=str )#skiprows=1# ))
T,g = data1[:,0], data1[:,3]





# curve fit
fig, ax = plt.subplots()
#Fit,_ = curve_fit(VanDerWallsGL, v, p,[0.00078,0.00008])
#fitParameter = Fit[0]
#fehlerMatrix = Fit[1]
a1,b1 = Fit
plt.scatter(T,g )
vFunk = np.linspace(0.8, 4, 20)
pFunk = VanDerWallsGL(vFunk, a1, b1)
#for i in range(len(fitParameter)):
 #   print('Parameter',i+1,':',fitParameter[i],'+-', np.sqrt(fehlerMatrix[i][i]))



plt.plot(vFunk, pFunk, '--', color='gray')

#verschönerung

plt.title('')
ax.set_xlabel('T')
ax.set_ylabel('g')
plt.show()
print('a= ',a1,'b= ', b1)
