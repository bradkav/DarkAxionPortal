import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc_file("matplotlibrc")
matplotlib.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}',
    r'\usepackage{amssymb}']


tfs = 14.0
#lfs = 20.0
nfs = 13.0
plt.rc('font', family='serif',size=tfs)
plt.rcParams['axes.linewidth'] = 2.0

import DarkAxionPortal
from Units import *


#------------------------------
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f_dp", "--f_dp",help = "DM fraction in Dark Photons: f_dp = [0, 1]", type=float, default=1.0)
args = parser.parse_args()
f_dp = args.f_dp

AXION_LIMITS = 2

DAP = DarkAxionPortal.Model(PQ_Phi = 1.0, e_D = 0.1, D_psi = 3)

print("> Plotting for f_dp = ", f_dp)
#------------------------------

def add_T_RH_label(T_RH, first=False):
    T_RH_str = "10^{" + str(int(np.log10(T_RH/GeV))) + "}\,\mathrm{GeV}$"
    if (first):
        x0 = 3e2
        T_RH_str = r"$T_\mathrm{RH} = " + T_RH_str
    else:
        x0 = 9e2
        T_RH_str = r"$" + T_RH_str
    y0 = 1.11*DAP.calc_fa_DPDensity_GluonFusion(x0*eV, f_dp, T_RH)/GeV
    if (y0 < 5e11):
        ax.text(x0,y0,T_RH_str,rotation=33, rotation_mode='anchor', fontsize=12, zorder=0)

T_RH_list = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 1000])*1e9*GeV

m_dp_list = np.geomspace(1, 1e9, 100)*eV

rescale = 0.85

#-------------------------------

fig = plt.figure(figsize=(6*rescale,5*rescale))
ax = plt.gca()

ax.set_xscale('log')
ax.set_yscale('log')

for i, T_RH in enumerate(T_RH_list):
    if ((i > 0) and (i < 9)):
        lw = 0.25
    else:
        lw = 1.5
    #ax.plot(m_dp_list/eV, DAP.calc_fa_DPDensity(m_dp_list, f_dp, T_RH)/GeV, color='k', linewidth=lw)
    fa_list = DAP.calc_fa_DPDensity_GluonFusion(m_dp_list, f_dp, T_RH)
    mask1 = T_RH < fa_list #Post-inflationary 
    mask2 = T_RH > fa_list #Pre-inflationary
    ax.plot(m_dp_list[mask1]/eV, fa_list[mask1]/GeV, color='k', linestyle='-', linewidth=lw)
    ax.plot(m_dp_list[mask2]/eV, fa_list[mask2]/GeV, color='k', linestyle='--', linewidth=lw)
    
#TO-DO: Implement as a loop
add_T_RH_label(1e9*GeV, first=True)
add_T_RH_label(1e10*GeV)
add_T_RH_label(1e11*GeV)
add_T_RH_label(1e12*GeV)


col_lepto = "mediumpurple"
ax.fill_between(m_dp_list/eV,0,DAP.calc_fa_DPDensity_GluonFusion(m_dp_list, f_dp, 1e9)/GeV, color=col_lepto,zorder=0, alpha=0.5)
ax.text(3e8,1.4e9,r"$T_\mathrm{RH} \lesssim 10^9$"+' GeV disfavoured by theory', zorder=10, fontsize=12.0, weight='bold', ha='right')

#Adding relic density lines for axions
f_a_min = DAP.calc_fa_AxionDensity(f_ax=(1-f_dp), theta_i = np.pi/DAP.N_ax)
f_a_med = DAP.calc_fa_AxionDensity(f_ax=(1-f_dp), theta_i = np.pi/(np.sqrt(3)*DAP.N_ax))
#f_a_BR = DAP.calc_fa_AxionDensity_full(f_ax=(1-f_dp), theta_i = 0.0, H_I=1e6*GeV)
#print(f_a_BR/(1e9*GeV))
#print(DAP.calc_Omega_Axion(f_a_BR, theta_i = 0.0, gamma = 1, H_I = 1e6*GeV)/DAP.Omegah2_dm)
#f_a_max = DAP.calc_fa_AxionDensity(f_ax=(1-f_dp), theta_i = np.pi/4.0)
#print(f_a_min/f_a_min_full)

col_a = "darkolivegreen"
if (AXION_LIMITS == 1):
    
    #plt.fill_between(m_dp_list/eV, f_a_min/GeV, 1e13, color='red', alpha=0.5)
    plt.axhline(f_a_min/GeV, linestyle='-', color=col_a)
    if (1e9 < f_a_min/GeV < 1e12):
        #plt.text(8e5*(1-2e-2), 1.58*f_a_min/GeV*(1+5e-3), r"$f_\mathrm{ax} = " + str((1-f_dp)) + "$ " +  r"$(\theta_i = \pi)$", color='white', fontsize=12, va='top')
        plt.text(8e5, 1.58*f_a_min/GeV, r"$f_\mathrm{ax} = " + str((1-f_dp)) + "$ " +  r"$(\theta_i = \pi)$", color=col_a, fontsize=nfs, va='top')

    plt.axhline(f_a_med/GeV, linestyle='-', color=col_a)
    if (1e9 < f_a_med/GeV < 1e12):
        props = dict(boxstyle='square,pad=0.18',facecolor='white', alpha=0.5, edgecolor='none')    
        plt.text(8e5, 1.58*f_a_med/GeV, r"$f_\mathrm{ax} = " + str((1-f_dp)) + "$ " +  r"$(\theta_i = \frac{\pi}{\sqrt{3}})$", color=col_a, fontsize=nfs, va='top')#, bbox=props)

    #plt.axhline(f_a_BR/GeV, linestyle='-', color=col_a)
    #if (1e9 < f_a_BR/GeV < 1e12):
        #props = dict(boxstyle='square,pad=0.18',facecolor='white', alpha=0.5, edgecolor='none')
    #    plt.text(8e5, 1.58*f_a_BR/GeV, r"$f_\mathrm{ax} = " + str((1-f_dp)) + "$ " +  r"$(\theta_i = 0, BR)$", color=col_a, fontsize=12, va='top')#, bbox=props)

elif (AXION_LIMITS == 2):
    col_a = "darkolivegreen"
    f1, f2 = DAP.calc_fa_AxionDensity_Buschmann(f_ax = (1-f_dp))
    plt.fill_between(m_dp_list/eV, f1/GeV, f2/GeV, color=col_a, alpha=0.7, zorder=4)
    if (1e9 < f1/GeV < 1e12):
        plt.text(7e8, np.sqrt(f1*f2)/GeV, r"$f_\mathrm{ax} = " + str((1-f_dp)) + "$ (post-inflationary)", fontsize=nfs, va='center', zorder=10, color='w', ha='right')
    
    #for fac in np.linspace(0, 0.3, 10):
    
    plt.axhline(f_a_min/GeV, linestyle='-', color=col_a, alpha=0.7, lw=2.5)
        
    #gradient_fill(m_dp_list/eV, m_dp_list*0.0 + f_a_min/GeV, color=col_a)
        
    if (1e9 < f_a_min/GeV < 1e12):
            plt.text(7e8, 0.9*f_a_min/GeV, r"$f_\mathrm{ax} = " + str((1-f_dp)) + "$ (pre-inflationary)", color='k', fontsize=nfs, va='top', ha='right')
            plt.annotate("", xytext = (5e6, f_a_min/GeV), xy =   (5e6, 1.5*f_a_min/GeV), arrowprops=dict(width=2.5, headlength=5.0,headwidth=9.0, shrink=0, facecolor=col_a, lw=0.0, shrinkA = 0, shrinkB = 0, alpha=0.7))
    
    
#Adding the 'not Freeze-In' region
m_dp_min = DAP.calc_m_dp_min(f_dp)
ax.axvline(m_dp_min/eV,linestyle='-',color="black",zorder=1000)
ax.fill_betweenx([1e9,1e12],m_dp_min/eV, facecolor='salmon', zorder=5)
if (m_dp_min/eV > 4.0):
    ax.text(0.75*m_dp_min/eV,3e9,"Freeze-in not achieved", weight='bold',rotation=90, rotation_mode='anchor',zorder=10, fontsize=14.0)


plt.xlabel(r"$m_{\gamma'}$ [eV]")
plt.ylabel(r"$f_a$ [GeV]")

plt.xlim(1, 1e9)
plt.ylim(1e9, 1e12)

locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=50)
locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=100)
ax.xaxis.set_major_locator(locmaj)
ax.xaxis.set_minor_locator(locmin)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

plt.tick_params(axis='x', which = 'both', top=False)
plt.tick_params(axis='y', which = 'both', right=False)

secax = ax.secondary_yaxis('right', functions=(lambda x: DAP.fa_to_ma(x)/eV, lambda x: DAP.ma_to_fa(x)))
secax.set_ylabel(r'$m_a$ [eV]')

plt.title(r"$f_{\gamma'} = \Omega_{\gamma'}/\Omega_\mathrm{DM} = " + str(float(f_dp)) + " $", pad=10)

plt.savefig("plots/FreezeInRegion_fdp=" + str(int(f_dp*100)) + "pct.pdf", bbox_inches='tight')
#plt.show()