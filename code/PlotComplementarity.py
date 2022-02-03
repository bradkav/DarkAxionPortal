import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
rootdir = os.path.dirname(os.path.abspath(__file__)) + "/"

matplotlib.rc_file(rootdir + "matplotlibrc")
matplotlib.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}',
    r'\usepackage{amssymb}']

tfs = 16.0
lfs = 18.0
nfs = 12.0
plt.rc('font', family='serif',size=tfs)
plt.rcParams['axes.linewidth'] = 2.0


import DarkAxionPortal 
from Units import *





#------------------------------
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f_dp", "--f_dp",help = "DM fraction in Dark Photons: f_dp = [0, 1]", type=float, default=0.1)
parser.add_argument("-e_D", "--e_D",help = "Coupling constant of Dark U(1), e_D", type=float, default=0.1)
parser.add_argument("-show_e_D", "--show_e_D", help = "Include value of e_D in plot", type=int, default=0)
args = parser.parse_args()
f_dp = args.f_dp
e_D  = args.e_D
SHOW_e_D = False
if (args.show_e_D > 0):
    SHOW_e_D = True

DAP = DarkAxionPortal.Model(PQ_Phi = 1.0, e_D = e_D, D_psi = 3)

print("g_agammagamma:", DAP.calc_G_app(1))

AXION_LIMITS = 2

print("> Plotting for f_dp = ", f_dp, "; e_D = ", e_D)
#------------------------------

T_RH_list = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 1000])*1e9*GeV

m_dp_list = np.geomspace(1, 1e9, 200)*eV

rescale = 0.95

fig = plt.figure(figsize=(6*rescale,5*rescale))
ax = plt.gca()

ax.set_xscale('log')
ax.set_yscale('log')

ax.tick_params(which='major',width=1.5)
ax.tick_params(which='minor',width=1)

for loc in ["bottom", "top", "left", "right"]:
    ax.spines[loc].set_zorder(1e10)

    
#Adding the 'not Freeze-In' region
m_dp_min = DAP.calc_m_dp_min(f_dp)
PLOT_FREEZE_IN = False
if (PLOT_FREEZE_IN):
    ax.axvline(m_dp_min/eV,linestyle='-',color="black",zorder=1000)
    ax.fill_betweenx([1e-6,1e-3],m_dp_min/eV, facecolor='salmon', zorder=5)
    if (m_dp_min/eV > 4.0):
        ax.text(0.3*m_dp_min/eV,3e-6,"Freeze-in not achieved", weight='bold',rotation=90, rotation_mode='anchor',zorder=10, fontsize=nfs)


#Adding the 'DAMIC-M' region
ax.fill_between([10, 200],1e-10, 1, facecolor='yellow', zorder=5, alpha=0.7)
ax.axvline(10,linestyle='-',color="black",zorder=5.1)
ax.axvline(200,linestyle='-',color="black",zorder=5.1)
#if (m_dp_min/eV > 4.0):
ax.text(18, 4.5e-6,"DAMIC-M\nprojected", weight='bold',rotation=90, rotation_mode='anchor',zorder=10, fontsize=nfs, va='center')

#Add ADMX region
#ax.fill_between(m_dp_list/eV, 2e-6, 4e-6, facecolor=[0.8, 0.0, 0.0], zorder=4)
ax.fill_between(m_dp_list/eV, 2e-6, 4e-6, facecolor=[0.8, 0.0, 0.0], zorder=100002)
ax.axhline(2e-6,linestyle='-',color="black",zorder=100002)
ax.axhline(4e-6,linestyle='-',color="black",zorder=100002)
ax.text(1e3, 2.7e-6,"ADMX excluded", weight='bold',rotation=0, rotation_mode='anchor',zorder=100002, fontsize=nfs, va='center')

#Add Axion projected region
ax.fill_between(m_dp_list/eV, 7e-7, 400e-6, facecolor=[0.8, 0.0, 0.0], zorder=3, alpha=0.7)
ax.axhline(7e-7,linestyle='-',color="black",zorder=3.1)
ax.axhline(400e-6,linestyle='-',color="black",zorder=3.1)
ax.text(1.2, 30e-6,"Haloscopes\nprojected", weight='bold',rotation=0, rotation_mode='anchor',zorder=10, fontsize=nfs, va='center')

#ax.axhline(,linestyle='-.',color="black",zorder=1000)
#ax.axhline(200,linestyle='-.',color="black",zorder=1000)
#if (m_dp_min/eV > 4.0):
#ax.text(150, 3e-6,"DAMIC-M sensitive", weight='bold',rotation=90, rotation_mode='anchor',zorder=10, fontsize=14.0)

#Add Freeze-in region
#def calc_fa_max(m_dp):
        

fa_list1 = DAP.calc_fa_DPDensity_GluonFusion(m_dp_list, f_dp, T_RH=1e9*GeV)
fa_list2 = DAP.calc_fa_DPDensity_GluonFusion(m_dp_list, f_dp, T_RH=1)**(4)
ma_list1 = DAP.fa_to_ma(fa_list1)
ma_list2 = DAP.fa_to_ma(fa_list2)
#mask2 = 1e13*GeV < fa_list2 #Pre-inflationary

col_DP = [0, 0, 0.9]

inds = np.arange(len(m_dp_list))[m_dp_list > m_dp_min]
ax.fill_between(m_dp_list[inds]/eV, ma_list1[inds]/eV, ma_list2[inds]/eV, facecolor=col_DP, zorder=4, alpha=0.7)
xDP = 3e3
yDP = 3e-5
if (e_D < 0.05):
    xDP = 6e3
    yDP = 2e-4
plt.text(xDP, yDP, "Dark Photon\nproduction", color='white', fontsize=nfs, va='center', zorder=100)

plt.plot(m_dp_list[inds]/eV, ma_list1[inds]/eV, color='black', zorder=4.1)
plt.plot(m_dp_list[inds]/eV, ma_list2[inds]/eV, color='black', zorder=4.1)
plt.plot([m_dp_list[inds[0]]/eV, m_dp_list[inds[0]]/eV], [ma_list1[inds[0]]/eV, ma_list2[inds[0]]/eV], color='black', zorder=4.1)

#Add Axion production
#col_a = "blue"
#col_a = [ 0, 0, 0.9]
col_a = "cornflowerblue"
#col_a = "darkolivegreen"
#f1, f2 = DAP.calc_fa_AxionDensity_Buschmann(f_ax = (1-f_dp))
#plt.fill_between(m_dp_list/eV, f1/GeV, f2/GeV, color=col_a, alpha=0.7, zorder=4)
#if (1e9 < f1/GeV < 1e12):
#    plt.text(3e2, np.sqrt(f1*f2)/GeV, r"$f_\mathrm{ax} = " + str((1-f_dp)) + "$ (post-inflationary, Buschmann et al.)", fontsize=10, va='center', zorder=10, color='w')

#for fac in np.linspace(0, 0.3, 10):
f_a_min = DAP.calc_fa_AxionDensity(f_ax=(1-f_dp), theta_i = np.pi/DAP.N_ax)
plt.axhline(DAP.fa_to_ma(f_a_min)/eV, linestyle='-', color=col_a, lw=2.5, zorder=25)
    
#gradient_fill(m_dp_list/eV, m_dp_list*0.0 + f_a_min/GeV, color=col_a)
    
if (1e9 < f_a_min/GeV < 1e12):
        plt.text(4.5e3, 1.1*DAP.fa_to_ma(0.67*f_a_min)/eV, r"Axion production", color=col_a, fontsize=nfs*0.9, va='top', zorder=10)
        plt.annotate("", xytext = (1e3, DAP.fa_to_ma(f_a_min)/eV), xy =   (1e3, DAP.fa_to_ma(1.8*f_a_min)/eV), arrowprops=dict(width=2.5, headlength=5.0,headwidth=9.0, shrink=0, facecolor=col_a, lw=0.0, shrinkA = 0, shrinkB = 0, alpha=1.0), zorder=25)
    

#Intersecting region:
inds = np.arange(len(m_dp_list))[(m_dp_list <  200*eV) & (m_dp_list > m_dp_min) & (m_dp_list > 10*eV)]
#print(inds)
lower = 1.15*np.maximum(ma_list2, np.ones_like(ma_list2)*7e-7*eV)
upper = 0.9*np.ones_like(m_dp_list)*400e-6*eV


if (np.any(lower[inds] < upper[inds])):
    #print(m_dp_list[inds[-1]])
    mnew = m_dp_list[inds]
    mnew[0] *= 1.05
    mnew[-1] *= 0.9
    mlist2 = np.concatenate((mnew , [mnew[-1]], mnew[::-1], [mnew[0]]))
    ylist = np.concatenate((lower[inds], [upper[inds[-1]]],upper[inds[::-1]], [lower[inds[0]]]))

    #plt.fill_between(mnew/eV, lower[inds]/eV, upper[inds]/eV, color='white', zorder=50)

    #plt.plot(mlist2/eV, ylist/eV, lw=3, zorder=100000, alpha=1.0, color='white')
    plt.plot(mlist2/eV, ylist/eV, lw=2, zorder=100001, alpha=1.0, color='black', linestyle='--')



plt.xlabel(r"Dark Photon mass $m_{\gamma^\prime}$ [eV]",fontsize=lfs)
plt.ylabel(r"Axion mass $m_a$ [eV]",fontsize=lfs)

plt.xlim(1, 1e5)
plt.ylim(3e-7, 3e-3)

locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=50)
locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=100)
ax.xaxis.set_major_locator(locmaj)
ax.xaxis.set_minor_locator(locmin)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

plt.tick_params(axis='x', which = 'both', top=True)
plt.tick_params(axis='y', which = 'both', right=True)


plottext = r"$f_{\mathrm{ax}} = " + str(int(100*float(1-f_dp))) + "\%; \,\, f_{\gamma'} = " + str(int(100*float(f_dp))) + "\% $"
if (SHOW_e_D):
    plottext += r"; $e^\prime = " + str(float(e_D)) + "$"

props = dict(boxstyle='round,pad=0.18',facecolor='white', alpha=0.9, edgecolor='none')
plt.text(8.2e4, 1.9e-3, plottext, ha='right',va='center', fontsize=12.0, bbox=props, zorder=50)
outfile = rootdir + "../plots/Complementarity_" + str(int(f_dp*100)) + "pct"
if (SHOW_e_D):
    outfile += "_eD_" + str(e_D)

#plt.title(r"$f_{\mathrm{ax}} = " + str(int(100*float(1-f_dp))) + "\%; \,\, f_{\gamma'} = " + str(int(100*float(f_dp))) + "\% $", pad=10, fontsize=14.0)

plt.savefig(outfile + ".pdf", bbox_inches='tight')
#plt.show()