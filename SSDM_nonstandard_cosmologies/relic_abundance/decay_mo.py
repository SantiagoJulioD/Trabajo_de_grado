#!/usr/bin/env python

#import pandas as pd
import numpy as np
#from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm
from numpy import *
#from scipy.integrate  import quad
import cmath
import pandas as pd

import subprocess
import sys
import os


os.chdir('/home/santiago/MicrOmegas/micromegas_6.2.3/micromegas_6.2.3/SingletDM')
def writeinputf(file,dictionary):
    ''' write dictionarys parameters in file'''
    data1 = open(file,'w')
    for items in dictionary.items():
        data1.write("%s %s\n" % items)

    data1.close()


newlist1=[ ]

datapar={"Mh": 125.2, "Mdm1": 50., "laSH": 0.1}

sltns=0

MH   =  125.2
Mdm1  = 1.
laSH  = 1e-11   

datapar['Mh']    = MH
datapar['Mdm1']   = Mdm1
datapar['laSH']  = laSH

sltns=0

def obtain_Gamma(M,lam):    
    datapar['Mdm1'] = M
    datapar['laSH'] = lam
    writeinputf('data.dat',datapar)

    subprocess.getoutput("./main data.dat > output.dat")
    aa=subprocess.getoutput("grep 'h -> ~x1,~x1' output.dat | awk -F'h' '{print $1}'")
    bb=subprocess.getoutput("grep 'total width' output.dat | awk -F'=' '{print $2}'")
    return float(aa)*float(bb)    

obtain_Gamma = np.vectorize(obtain_Gamma)

Ms = np.logspace(0,np.log10(MH/2))
la = 1e-11 # las = np.logspace(-13,-10,10)

# MM, LL = np.meshgrid(Ms,las)
OO = obtain_Gamma(Ms,la)

result = pd.DataFrame()
result['M'] = Ms
result['Gamma'] = OO

result.to_csv('Gamma_mo.csv')

fig, ax = plt.subplots(figsize=(8,6))
ax.loglog(Ms,OO)
#CS = ax.contourf(np.log10(MM),np.log10(LL),np.log10(OO),cmap='inferno',)
#CSl = ax.contour(CS,levels=[np.log10(0.12),],colors='red',linestyles='solid',linewidths=2)
#ax.vlines(np.log10(m_h/2),-13,-10)

# p = CSl.collections[0].get_paths()[0]
# v = p.vertices
# x = v[:,0]
# y = v[:,1]
# params = pd.DataFrame()
# params['M'] = x
# params['lambda'] = y
# params.to_csv('params_micromegas.csv')

plt.show()





# for irun in range(0,10001,1):
#     if (irun%1000==0):
#         print ('irun=',irun)
#     omegaS = -1.0
#     xf     = -1.0
#     Ms2  = 10**( (log10(1e4)-log10(5e1))*np.random.uniform(0,1)+log10(5e1))
#     laS  = 10**( (log10(1e0)-log10(1e-4))*np.random.uniform(0,1)+log10(1e-4))
#     # writing data.dat input file for micromegas
#     datapar['Ms2']   = Ms2
#     datapar['laS']  = laS
#     writeinputf('data.dat',datapar)
#     # running micromegas and extracting the relic density (omega)
#     subprocess.getoutput("./main data.dat > output.dat")
#     aa=subprocess.getoutput("grep 'Omega' output.dat | awk -F'=' '{print $3}'")
#     bb=subprocess.getoutput("grep 'Omega' output.dat | awk -F'=' '{print $2}' | awk -F' ' '{print $1}'")
#     omegaS = float(aa)
#     xf     = float(bb)
#     if (omegaS > 0.0):
#         sltns=sltns+1
#         newlist1.append([Mp1,Ms2,omegaS,laS,ys,yp,xf])
# print ('sltns=',sltns)                                                                                               
# datos = np.asarray(newlist1)
# np.savetxt('nuevosdatos.txt',datos)