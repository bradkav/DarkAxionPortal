import numpy as np
import matplotlib.pyplot as plt

import PlotFuncs_DarkPhoton as PF

import DarkAxionPortal
from Units import *

DAP = DarkAxionPortal.Model(PQ_Phi = 1.0, e_D = 0.1, D_psi = 3)

def KineticMixing_production():

    m_dp=np.linspace(1e-3,1e5,10000)

    #Production through resonant conversion
    #-------------------------------------
    eps100=DAP.calc_epsilon_Conversion(m_dp*eV, f_dp=1.0)
    #chi10=DAP.calc_epsilon_Conversion(m_dp*eV, 1e-1)
    eps1=DAP.calc_epsilon_Conversion(m_dp*eV, f_dp=1e-2)
    #chi01=DAP.calc_epsilon_Conversion(m_dp*eV, 1e-3)
    
    
    plt.plot(m_dp, eps100, color='white', linestyle=':', lw=2.5)
    #plt.plot(m_dp, chi10, color='black', linestyle=':')
    plt.plot(m_dp, eps1, color='white', linestyle=':', lw=2.5)
    #[m_dp < 3e4]
    plt.plot(m_dp, eps1, color='grey', linestyle='-', lw=2.5, zorder=-1)
    plt.plot(m_dp, eps100, color='grey', linestyle='-', lw=2.5, zorder=-1)
    
    x_adjust = 500
    plt.text(x_adjust*0.20*(1+1e-2), (x_adjust**(-3/2))*5e-6*(1-1e-2), r"$f_{\gamma'} = 100\%$ ($\gamma \rightarrow \gamma'$)", ha = "left", fontsize=18, rotation=-38.5, color='black', rotation_mode='anchor')
    plt.text(x_adjust*0.20, (x_adjust**(-3/2))*5e-6, r"$f_{\gamma'} = 100\%$ ($\gamma \rightarrow \gamma'$)", ha = "left", fontsize=18, rotation=-38.5, color='white', rotation_mode='anchor')
    #plt.text(1, 1e-7, r"$f_{\gamma'} = 10\%$", ha = "left")
    plt.text(x_adjust*0.13*(1+1e-2), (x_adjust**(-3/2))*1.5e-7*(1 - 1e-2), r"$f_{\gamma'} = 1\%$", ha = "left", fontsize=18, rotation=-38.5, color='black')
    plt.text(x_adjust*0.13, (x_adjust**(-3/2))*1.5e-7, r"$f_{\gamma'} = 1\%$", ha = "left", fontsize=18, rotation=-38.5, color='white')
    
    
    #Production through Dark Higgs-strahlung
    #-------------------------------------
    eps100_DH = DAP.calc_epsilon_DarkHiggsstrahlung(m_dp*eV, f_dp=1.0)
    eps1_DH = DAP.calc_epsilon_DarkHiggsstrahlung(m_dp*eV, f_dp=1e-2)
    
    plt.plot(m_dp, eps100_DH, color='white', linestyle='--', lw=2.5)
    #plt.plot(m_dp, chi10, color='black', linestyle=':')
    plt.plot(m_dp, eps1_DH, color='white', linestyle='--', lw=2.5)
    
    x_adjust = 220
    rot = -14.3
    plt.text(x_adjust*0.5*(1+1e-2), (x_adjust**(-1/2))*7.5e-7*(1-1e-2), r"$f_{\gamma'} = 100\%$ (Dark Higgs-strahlung)", ha = "left", fontsize=18, rotation=rot, color='black', rotation_mode='anchor')
    plt.text(x_adjust*0.5, (x_adjust**(-1/2))*7.5e-7, r"$f_{\gamma'} = 100\%$ (Dark Higgs-strahlung)", ha = "left", fontsize=18, rotation=rot, color='white', rotation_mode='anchor')
    #plt.text(1, 1e-7, r"$f_{\gamma'} = 10\%$", ha = "left")
    plt.text(x_adjust*0.4*(1+1e-2), (x_adjust**(-1/2))*4e-8*(1 - 1e-2), r"$f_{\gamma'} = 1\%$", ha = "left", fontsize=18, rotation=rot, color='black')
    plt.text(x_adjust*0.4, (x_adjust**(-1/2))*4e-8, r"$f_{\gamma'} = 1\%$", ha = "left", fontsize=18, rotation=rot, color='white')

def Decay():

    m_dp=np.linspace(1e-3,1e5,10000)
    plt.fill_between(m_dp, 5.8e8*(m_dp/10)**(-9/2),1e2, color='grey', linestyle='-', lw=2.5, zorder=5)
    plt.plot(m_dp, 5.8e8*(m_dp/10)**(-9/2), color='black', linestyle='-', lw=2.5, zorder=5)
    
    #plt.text(8e3*(1+1e-2), 1e-6*(1-1e-2), r"$\tau_{\gamma^\prime} < \tau_\mathrm{U}$", ha = "left", fontsize=18, rotation=-45, color='black', rotation_mode='anchor', zorder=10)
    plt.text(3e4, 1e-6, r"$\tau_{\gamma^\prime} < \tau_\mathrm{U}$", ha = "left", fontsize=18, rotation=-63, color='black', rotation_mode='anchor', zorder=10)
    
    
#--------------------------------------------------------------------------------------

fig,ax = PF.FigSetup(Shape="Square", chi_max = 1e-5, m_min = 1e-2)

# DPDM
#DarkMatter(ax)

# Axion haloscopes
#Haloscopes(ax)

# LSW/Helioscopes
#LSW(ax)
PF.CAST(ax)
PF.SHIPS(ax)

# Tests of coulomb law
#Coulomb(ax)

# Reactor neutrinos
#TEXONO(ax)

# DPDM searches
PF.Xenon(ax)
PF.DAMIC(ax, text_on = False)
plt.text(16, 1e-13, "DAMIC", fontsize=20, rotation=-90)
PF.SENSEI(ax)
PF.FUNK(ax, text_on = False)
plt.text(2, 3e-12, "FUNK", fontsize=20, rotation=-90)
#Tokyo(ax)
#SHUKET(ax)
#DarkEfield(ax)
#WISPDMX(ax)
#SQuAD(ax)
#DMPathfinder(ax)

# Astrophysical bounds
PF.StellarBounds(ax, Higgsed=True, e_D = DAP.e_D)
#COBEFIRAS(ax)
#Jupiter(ax)
#Earth(ax)
#Crab(ax)
#IGM(ax)
#LeoT(ax)

KineticMixing_production()
Decay()
#DAP_projections()

plt.gcf().text(0.89,0.13,r'$\rho_0 = 0.45$ GeV cm$^{-3}$',fontsize=25,ha='right', va='bottom')

PF.MySaveFig(fig,'DarkPhoton_KineticMixing_Production')
