import numpy as np
import matplotlib.pyplot as plt

import PlotFuncs_DarkPhoton as PF

import DarkAxionPortal
from Units import *

import os

rootdir = os.path.dirname(os.path.abspath(__file__)) + "/"

DAP = DarkAxionPortal.Model(PQ_Phi = 1.0, e_D = 0.1, D_psi = 3)

f_DP = 0.1
    
def DAP_projections(f_DP=1.0):
    A_Si=28.0855 #Atomic mass of Silicon
    N_a=6.02214076e23 #Avogadro
    density=f_DP*4.5e8 #local energy-density of dark matter= 0.3 GeV/cm^3
    
    #Cross section in Mbarns
    PECross = np.loadtxt(rootdir + '../data/Si_PhotoelectricAbsorptionCrossSection.txt')
    PECross[:,1]=PECross[:,1]*10**-15*N_a/A_Si #To put it in cm^2 kg^-1
    
    y2 = ax.get_ylim()[1]
    #plt.fill_between(PECross[:,0],np.sqrt((3/((density/PECross[:,0])*PECross[:,1]*2.5902e15*5))),y2=y2,edgecolor='k',linewidth=2.5, linestyle="-.", facecolor="red",zorder=0.01, alpha=0.75)
    plt.fill_between(PECross[:,0],np.sqrt((0.1/((density/PECross[:,0])*PECross[:,1]*2.5902e15*5))),y2=y2,edgecolor='k',linewidth=2.5, linestyle="-.", facecolor="sienna",zorder=0.001, alpha=0.75)
    
    #plt.text(1.0,4e-16,r'{\bf LBC}'+'\n'+'(This Work)',fontsize=21,color="red",rotation=0,rotation_mode='anchor',ha='center',va='center') 
    plt.text(2.5,1e-17/np.sqrt(f_DP),r'{\bf DAMIC-M}'+'\n'+'(This Work)',fontsize=21,color="sienna",rotation=0,rotation_mode='anchor',ha='center',va='center')
    #plt.plot([3e0,8],[6e-16,8e-16],'k-',lw=2.5,color="red")
    plt.plot([7e0,9],[3e-17/np.sqrt(f_DP),5e-17/np.sqrt(f_DP)],'k-',lw=2.5,color="sienna")

    
    
#--------------------------------------------------------------------------------------

fig,ax = PF.FigSetup(Shape="Square", chi_max = 1e-8, m_min = 1e-2, lfs=40, tfs=35, upper_xlabel=r"Easter Egg")

# DPDM
#DarkMatter(ax)

# Axion haloscopes
#Haloscopes(ax)

# LSW/Helioscopes
#LSW(ax)
PF.CAST(ax, text_on=False)
plt.text(1e3*(1-0.01),1e-9*(1+0.08),r'{\bf CAST}',fontsize=27,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')   
plt.text(1e3,1e-9,r'{\bf CAST}',fontsize=27,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')

PF.SHIPS(ax)

# Tests of coulomb law
#Coulomb(ax)

# Reactor neutrinos
#TEXONO(ax)



# DPDM searches
PF.Xenon(ax, f_DP=f_DP, text_on=False)
plt.text(8e2,1e-16,r'{\bf XENON}',fontsize=22,color='crimson',rotation=0,rotation_mode='anchor',ha='center',va='center')

PF.DAMIC(ax, text_on = False, f_DP=f_DP)
plt.text(16, 0.4e-13/np.sqrt(f_DP), "DAMIC", fontsize=20, rotation=-90)
PF.SENSEI(ax, f_DP=f_DP, text_on = False)

plt.text(1e0,3e-14, "SENSEI", fontsize=20, rotation=0, color="Firebrick")
plt.plot([5e0,8e0],[3.8e-14,5.3e-14],'k-',lw=2.5,color="Firebrick")



PF.FUNK(ax, text_on = False, f_DP=f_DP)
plt.text(2, 3e-12/np.sqrt(f_DP), "FUNK", fontsize=20, rotation=-90)
#Tokyo(ax)
#SHUKET(ax)
#DarkEfield(ax)
#WISPDMX(ax)
#SQuAD(ax)
#DMPathfinder(ax)

# Astrophysical bounds
PF.StellarBounds(ax, Higgsed=True, e_D=DAP.e_D)
#COBEFIRAS(ax)
#Jupiter(ax)
#Earth(ax)
#Crab(ax)
#IGM(ax)
#LeoT(ax)

#KineticMixing_production()
#Decay()
DAP_projections(f_DP)

plt.gcf().text(0.89,0.13,r'$f_{\gamma^\prime} = '+ str(int(f_DP*100)) + r'\%; \,\,\rho_0 = 0.45$ GeV cm$^{-3}$',fontsize=25,ha='right',va='bottom')

PF.MySaveFig(fig,'DarkPhoton_KineticMixing_Projections_fDP_'+str(int(f_DP*100)) + 'pct')