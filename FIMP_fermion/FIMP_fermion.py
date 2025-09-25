import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import kn
from scipy.interpolate import interp1d
import astropy.constants as ct
import pandas as pd
import subprocess
import sys
import os

# Global variables

MP = 2.4e18
T0 = 2.725*ct.k_B.value/ct.e.value*1e-9
gstars0 = 3.91
gstar0 = 3.38
s0 = 2*np.pi**2/45*gstars0*T0**3
gS = 1

# Classes

class FIMP_fermion:
    '''
    Class for fermionic FIMP dark matter with Z4 symmetry
    https://arxiv.org/pdf/2308.05249

    We consider case 1: S is thermal and can decay into two fermions
    '''
    def __init__(self,M_fer,M_bos,MH,laS,ys,yp):
        self.M_fer = M_fer
        self.M_bos = M_bos
        self.MH = MH
        self.laS = laS
        self.ys = ys
        self.yp = yp
        self.obs_omega = 0.12

    def relic_abundance(self,small_coupling,micromegas):
        if small_coupling:
            Omegah2 = 2.744*(self.M_fer/1)*(1000/self.M_bos)*(self.ys*1e10)**2
        else:
            rdof = pd.read_csv('/home/santiago/Trabajo_de_grado/FIMP_fermion/rdof.csv')
            #rdof = pd.read_csv('https://raw.githubusercontent.com/SantiagoJulioD/Trabajo_de_grado/refs/heads/main/FIMP_fermion/rdof.csv')
            gstar_s = interp1d(rdof['Temp'][::-1],rdof['g_s'][::-1],bounds_error=False,fill_value=(rdof['g_s'].values[-1],rdof['g_s'].values[0]))
            gstar_rho = interp1d(rdof['Temp'][::-1],rdof['g_rho'][::-1],bounds_error=False,fill_value=(rdof['g_rho'].values[-1],rdof['g_rho'].values[0]))
            alpha = self.ys**2*(1-4*self.M_fer**2/self.M_bos**2)**1.5+self.yp**2*(1-4*self.M_fer**2/self.M_bos**2)**0.5
            integrand = lambda x: kn(1,x)*x**3/((gstar_rho(self.M_bos/x))**0.5*gstar_s(self.M_bos/x))
            I = quad(integrand,self.M_bos*1e-12,np.inf)[0]
            const = gS*45*MP*alpha/((np.pi**2/90)**0.5*16*np.pi**5*self.M_bos)  
            Y = const*I
            Omegah2 = s0*self.M_fer*Y/(3*MP**2*(2.13e-42)**2)
        if micromegas:
            os.chdir('/home/santiago/MicrOmegas/micromegas_6.2.3/micromegas_6.2.3/Z4fermion')
            def writeinputf(file,dictionary):
                ''' write dictionarys parameters in file'''
                data1 = open(file,'w')
                for items in dictionary.items():
                    data1.write("%s %s\n" % items)

                data1.close()
            if isinstance(self.M_fer,np.ndarray):
                Omegah2_micro = np.zeros_like(self.M_fer)
                for i in range(len(self.M_fer)):
                    if i%10==0:
                        print(i)
                    datapar = {"MH": self.MH, "Mp1": self.M_fer[i], "Ms2": self.M_bos, "laS": self.laS, "ys": self.ys, "yp": self.yp}          
                    writeinputf('data.dat',datapar)

                    subprocess.getoutput("./main data.dat > output.dat")
                    aa=subprocess.getoutput("grep 'omega freeze-in' output.dat | awk -F'=' '{print $2}'")

                    Omegah2_micro[i] = float(aa)

            else:
                datapar = {"MH": self.MH, "Mp1": self.M_fer, "Ms2": self.M_bos, "laS": self.laS, "ys": self.ys, "yp": self.yp}          
                writeinputf('data.dat',datapar)

                subprocess.getoutput("./main data.dat > output.dat")
                aa=subprocess.getoutput("grep 'omega freeze-in' output.dat | awk -F'=' '{print $2}'")

                Omegah2_micro[i] = float(aa)          

            return {'Omegah2_ODE':Omegah2,'Omegah2_MO':Omegah2_micro}
        else:
            return {'Omegah2_ODE':Omegah2}
        
    def plot_omega(self,small_coupling,micromegas,low_m,up_m,N_m,low_y,up_y,N_y,colors):

        masses = np.logspace(low_m,up_m,N_m)

        ys = np.logspace(low_y,up_y,N_y)

        self.M_fer = masses
        # print(ys)
        for i,y in enumerate(ys):
            self.ys = y
            self.yp = y
            Om = self.relic_abundance(small_coupling,micromegas)
            plt.loglog(masses,Om['Omegah2_ODE'],label=r'$\log y_{s,p}=$%.2f'%np.log10(self.ys),color=colors[i])
            if micromegas:
               plt.loglog(masses,Om['Omegah2_MO'],label=r'$\log y_{s,p}=$%.2f (micrOMEGAs)'%np.log10(self.ys),color=colors[len(ys)+i],linestyle=(0,(5,5)))
            

        plt.hlines(0.12,1e-3,1e3,colors='k')
        plt.grid()
        plt.xlabel(r'$M_\psi$')
        plt.ylabel(r'$\Omega h^2$')
        plt.ylim(1e-4,1e3)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        plt.show()

    def plot_y(self,low_mp,up_mp,N_mp,low_ms,up_ms,N_ms,colors):
        Mps = np.logspace(low_mp,up_mp,N_mp)
        MSs = np.logspace(low_ms,up_ms,N_ms)

        self.M_fer = Mps

        rdof = pd.read_csv('/home/santiago/Trabajo_de_grado/FIMP_fermion/rdof.csv')
        #rdof = pd.read_csv('https://raw.githubusercontent.com/SantiagoJulioD/Trabajo_de_grado/refs/heads/main/FIMP_fermion/rdof.csv')
        gstar_s = interp1d(rdof['Temp'][::-1],rdof['g_s'][::-1],bounds_error=False,fill_value=(rdof['g_s'].values[-1],rdof['g_s'].values[0]))
        gstar_rho = interp1d(rdof['Temp'][::-1],rdof['g_rho'][::-1],bounds_error=False,fill_value=(rdof['g_rho'].values[-1],rdof['g_rho'].values[0]))
        for i,M in enumerate(MSs):
            self.M_bos = M
            integrand = lambda x: kn(1,x)*x**3/((gstar_rho(self.M_bos/x))**0.5*gstar_s(self.M_bos/x))
            I = quad(integrand,self.M_bos*1e-12,np.inf)[0]
            const = gS*45*MP/((np.pi**2/90)**0.5*16*np.pi**5*self.M_bos)  
            yss = (self.obs_omega/(2.744e8*Mps*self.M_bos**3*(1-4*Mps**2/self.M_bos**2)**1.5*gS*const*I))**0.5
            yps = (self.obs_omega/(2.744e8*Mps*self.M_bos**3*(1-4*Mps**2/self.M_bos**2)**0.5*gS*const*I))**0.5
            plt.loglog(Mps,yss,color=colors[i])
            plt.loglog(Mps,yps,color=colors[i],linestyle=(0,(5,5)))
        plt.grid()
        plt.xlabel(r'$M_\psi$')
        plt.ylabel(r'$y_{p,s}$')
        #plt.ylim(1e-4,1e3)
        #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        plt.show()

                     

# subprocess.getoutput("./main data.dat > output.dat")
# aa=subprocess.getoutput("grep 'omega freeze-in' output.dat | awk -F'=' '{print $2}'")
# print(aa)