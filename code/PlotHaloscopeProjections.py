"""
Originally adapted from https://github.com/cajohare/AxionLimits/
"""

from numpy import *
import matplotlib.pyplot as plt

import sys
from Dirs import axionlimits_dir



from PlotFuncs import BlackHoleSpins, AxionPhoton, MySaveFig, UpperFrequencyAxis,col_alpha,FigSetup

fig,ax = FigSetup(Shape='Square',xlab = r'Axion mass $m_a$ [eV]',ylab='$|g_{a\gamma\gamma}|$ [GeV$^{-1}$]',mathpazo=True,\
                 m_min=1e-6,m_max=1e-3,g_min=2e-16,g_max=1e-9,lfs=40, tfs=35,xtick_rotation=0,FrequencyAxis=True,N_Hz=1e6,upper_xlabel=r"Frequency [MHz]",xlabel_pad=15)

y2 = 1e-9
#xscale = 1.03*(log10(3e-4) - log10(1e-6))/(log10(1e-3) - log10(1e-6))


f_ax = 0.9
scale = (f_ax)**-0.5

ADMX_col = 'crimson'


# CAST
AxionPhoton.Helioscopes(ax,projection=False,text_on=False)
plt.gcf().text(0.2*(1-0.005),0.805*(1+0.005),r'{\bf CAST}',color='k',fontsize=40,rotation=0,rotation_mode='anchor',ha='center',va='center')
plt.gcf().text(0.2,0.805,r'{\bf CAST}',color=col_alpha(ADMX_col,0.0),fontsize=40,rotation=0,rotation_mode='anchor',ha='center',va='center')

# QCD axion
AxionPhoton.QCDAxion(ax,text_on=False,C_center=abs(5/3-1.92)*(44/3-1.92)/2,C_width=0.7)
def g_x(C_ag,m_a):
        return 2e-10*C_ag*m_a
mvals = array([1e-6,1e-3])
plt.plot(mvals,g_x(44/3-1.95,mvals),'--',lw=3,color='brown',zorder=0.0,alpha=0.2)
plt.plot(mvals,g_x(abs(5/3-1.92),mvals),'--',lw=3,color='brown',zorder=0.1,alpha=0.2)
#plt.gcf().text(0.89,0.53,r'$E/N = 44/3$',fontsize=20,color='brown',rotation=13,rotation_mode='anchor',ha='right',alpha=0.9)
#plt.gcf().text(0.89,0.305,r'$E/N = 5/3$',fontsize=20,color='brown',rotation=13,rotation_mode='anchor',ha='right',alpha=0.9)
plt.gcf().text(0.89,0.49,r'KSVZ',fontsize=20,color='darkred',rotation=20,rotation_mode='anchor',ha='right',alpha=1)
plt.gcf().text(0.89,0.413,r'DFSZ',fontsize=20,color='darkred',rotation=20,rotation_mode='anchor',ha='right',alpha=1)



# RBF/UF
UF_col = [0.8,0,0]
RBF_col = 'Darkred'
dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/RBF_UF_Haloscopes.txt")
plt.fill_between(dat[:,0],scale*dat[:,1],y2=y2,edgecolor='k',facecolor='darkred',zorder=0.1)
dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/UF.txt")
plt.fill_between(dat[:,0],scale*dat[:,1],y2=y2,edgecolor='k',facecolor=[0.8,0,0],zorder=0.1)
plt.gcf().text(0.315*(1-0.003),0.32*(1+0.003),r'{\bf UF}',color='k',fontsize=30,rotation=0,rotation_mode='anchor',ha='center',va='center')
plt.gcf().text(0.315,0.32,r'{\bf UF}',color=col_alpha(UF_col,1),fontsize=30,rotation=0,rotation_mode='anchor',ha='center',va='center')
plt.gcf().text(0.355*(1-0.003),0.55*(1+0.003),r'{\bf RBF}',color='k',fontsize=30,rotation=0,rotation_mode='anchor',ha='center',va='center')
plt.gcf().text(0.355,0.55,r'{\bf RBF}',color=col_alpha(RBF_col,0.1),fontsize=30,rotation=0,rotation_mode='anchor',ha='center',va='center')


# ADMX
dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX.txt")
plt.fill_between(dat[:,0],scale*dat[:,1],y2=y2,edgecolor='k',facecolor=ADMX_col,zorder=0.1)
dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX2018.txt")
plt.fill_between(dat[:,0],scale*dat[:,1],y2=y2,edgecolor='k',facecolor=ADMX_col,zorder=0.1)
dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX2019_1.txt")
plt.fill_between(dat[:,0],scale*dat[:,1],y2=y2,edgecolor='k',facecolor=ADMX_col,zorder=0.1)
dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX2019_2.txt")
plt.fill_between(dat[:,0],scale*dat[:,1],y2=y2,edgecolor='k',facecolor=ADMX_col,zorder=0.1)
dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX2021.txt")
plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='k',facecolor=ADMX_col,zorder=0.1)
plt.gcf().text(0.22*(1-0.005),0.4*(1+0.005),r'{\bf ADMX}',color='k',fontsize=35,rotation=90,rotation_mode='anchor',ha='center',va='center')
plt.gcf().text(0.22,0.4,r'{\bf ADMX}',color=col_alpha(ADMX_col,0.2),fontsize=35,rotation=90,rotation_mode='anchor',ha='center',va='center')

# HAYSTAC
HAYSTAC_col = 'orangered'
dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/HAYSTAC_highres.txt")
dat2 = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/HAYSTAC_2020_highres.txt")
plt.fill_between(dat[:,0],scale*dat[:,1],y2=y2,edgecolor='k',facecolor=HAYSTAC_col,zorder=0.1)
plt.fill_between(dat2[:,0],scale*dat2[:,1],y2=y2,edgecolor='k',facecolor=HAYSTAC_col,zorder=0.1)
plt.gcf().text(0.436,0.357,r'{\bf HAYSTAC}',color=col_alpha(HAYSTAC_col,1),fontsize=17,rotation=90,rotation_mode='anchor',ha='center',va='center')
plt.gcf().text(0.473,0.42,r'{\bf HAYSTAC}',color=col_alpha(HAYSTAC_col,1),fontsize=17,rotation=90,rotation_mode='anchor',ha='center',va='center')


# CAPP
CAPP_col = [1, 0.1, 0.37]
dat1 = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/CAPP-1.txt")
dat2 = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/CAPP-2.txt")
dat3 = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/CAPP-3.txt")
plt.fill_between(dat1[:,0],scale*dat1[:,1],y2=y2,edgecolor='k',facecolor=CAPP_col,zorder=0.1)
plt.fill_between(dat2[:,0],scale*dat2[:,1],y2=y2,edgecolor='k',facecolor=CAPP_col,zorder=0.1)
plt.fill_between(dat3[:,0],scale*dat3[:,1],y2=y2,edgecolor='k',facecolor=CAPP_col,zorder=0.1)
plt.gcf().text(0.385,0.31,r'{\bf CAPP}',color=col_alpha(CAPP_col,1),fontsize=20,rotation=90,rotation_mode='anchor',ha='center',va='center')
plt.gcf().text(0.425,0.45,r'{\bf CAPP}',color=col_alpha(CAPP_col,1),fontsize=15,rotation=-90,rotation_mode='anchor',ha='center',va='center')
plt.gcf().text(0.347,0.35,r'{\bf CAPP}',color=col_alpha(CAPP_col,1),fontsize=15,rotation=-90,rotation_mode='anchor',ha='center',va='center')

# ORGAN
ORGAN_col = [0.7,0.2,0.05]
AxionPhoton.ORGAN(ax,text_on=False,col=ORGAN_col, f_ax = f_ax)
plt.gcf().text(0.645,0.615,r'{\bf ORGAN}',color=col_alpha(ORGAN_col,1),fontsize=20,rotation=90,rotation_mode='anchor',ha='center',va='center')

# GrAHal
GrAHal_col = '#b53e5a'
AxionPhoton.GrAHal(ax,text_on=False,col=GrAHal_col, f_ax = f_ax)
plt.gcf().text(0.497,0.51,r'{\bf GrAHal}',color=col_alpha(GrAHal_col,1),fontsize=15,rotation=-90,rotation_mode='anchor',ha='center',va='center')


# ADMX Sidecar
dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX_Sidecar.txt")
plt.fill_between(dat[:,0],scale*dat[:,1],y2=y2,edgecolor='k',facecolor=ADMX_col,zorder=0.1)
dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX_Sidecar_JTWPA.txt")
plt.fill_between(dat[:,0],scale*dat[:,1],y2=y2,edgecolor='k',facecolor=ADMX_col,zorder=0.1)
plt.gcf().text(0.471,0.55,r'{\bf ADMX}',color=ADMX_col,fontsize=18,rotation=90,rotation_mode='anchor',ha='center',va='center')


# QUAX
QUAX_col = [0.6,0.1,0.3]
dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/QUAX.txt")
dat2 = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/QUAX2.txt")
plt.plot([dat[0,0],dat[0,0]],[scale*dat[0,1],y2],color=QUAX_col,lw=2,zorder=0.05)
plt.plot([dat2[0,0],dat2[0,0]],[scale*dat2[0,1],y2],color=QUAX_col,lw=2,zorder=0.05)
plt.gcf().text(0.538,0.479,r'{\bf QUAX}-$a\gamma$',color=col_alpha(QUAX_col,1),fontsize=20,rotation=90,rotation_mode='anchor',ha='center',va='center')


# RADES
RADES_col = [0.95,0.1,0.1]
dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/RADES.txt")
plt.plot([dat[0,0],dat[0,0]],[scale*dat[0,1],y2],color=RADES_col,lw=2,zorder=0.05)
plt.gcf().text(0.517,0.53,r'{\bf RADES}',color=col_alpha(RADES_col,1),fontsize=15,rotation=90,rotation_mode='anchor',ha='center',va='center')


# Neutron stars
dat = loadtxt(axionlimits_dir +'limit_data/AxionPhoton/NeutronStars_GreenBank.txt')
plt.fill_between(dat[:,0],scale*dat[:,1],y2=y2,edgecolor='k',facecolor='ForestGreen',zorder=0.11)
#dat = loadtxt(axionlimits_dir +'limit_data/AxionPhoton/NeutronStars_VLA.txt')
#plt.fill_between(dat[:,0],scale*dat[:,1],y2=y2,edgecolor='k',facecolor='SeaGreen',zorder=0.11)
dat = loadtxt(axionlimits_dir +'limit_data/AxionPhoton/NeutronStars_Battye.txt') 
plt.plot(dat[:,0],dat[:,1],color='k',lw=1,zorder=0.2)
plt.fill_between(dat[:,0],scale*dat[:,1],y2=1e0,zorder=0.2,color='#4ad46a')
plt.gcf().text(0.50*(1-0.003),0.81*(1+0.003),r'{\bf Neutron Stars}',color='k',fontsize=30,rotation=0,rotation_mode='anchor',ha='center',va='center')
plt.gcf().text(0.50,0.81,r'{\bf Neutron Stars}',color=col_alpha('Green',0.4),fontsize=30,rotation=0,rotation_mode='anchor',ha='center',va='center')

plt.gcf().text(0.34*(1-0.001),0.76*(1+0.001),r'{\bf Foster et al.}',color='k',fontsize=25,rotation=0,rotation_mode='anchor',ha='center',va='center',multialignment='center')
plt.gcf().text(0.34,0.76,r'Foster et al.',color=col_alpha('ForestGreen',0.7),fontsize=25,rotation=0,rotation_mode='anchor',ha='center',va='center',multialignment='center')

#plt.gcf().text(0.47*(1-0.001),0.76*(1+0.001),r'{\bf Darling}',color='k',fontsize=25,rotation=0,rotation_mode='anchor',ha='center',va='center',multialignment='center')
#plt.gcf().text(0.47,0.76,r'Darling',color=col_alpha('SeaGreen',0.7),fontsize=25,rotation=0,rotation_mode='anchor',ha='center',va='center',multialignment='center')

plt.gcf().text(0.58*(1-0.001),0.76*(1+0.001),r'{\bf Battye et al.}',color='k',fontsize=25,rotation=0,rotation_mode='anchor',ha='center',va='center',multialignment='center')
plt.gcf().text(0.58,0.76,r'Battye et al.',color=col_alpha('#20c742',0.7),fontsize=25,rotation=0,rotation_mode='anchor',ha='center',va='center',multialignment='center')

#AxionPhoton.Haloscopes(ax,projection=True,fs=20,text_on=True,BASE_arrow_on=True)

text_on = True
AxionPhoton.ALPHA(ax,text_on=text_on, text_shift=[0.95, 1], f_ax = f_ax)
AxionPhoton.MADMAX(ax,text_on=text_on, text_shift=[1.4,1], f_ax = f_ax)
AxionPhoton.ORGAN(ax,projection=True,text_on=text_on, f_ax = f_ax)
#AxionPhoton.ADMX(ax,projection=True,text_on=text_on, f_ax = f_ax)

dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/Projections/ADMX_Projected.txt")
plt.plot(dat[:,0],(f_ax**-0.5)*dat[:,1],'-',linewidth=1.5,color=ADMX_col,zorder=0)
plt.fill_between(dat[:,0],(f_ax**-0.5)*dat[:,1],y2=y2,edgecolor=None,facecolor=ADMX_col,zorder=0,alpha=0.1)

plt.text(3.5e-5,1e-15,r'{\bf ADMX}',fontsize=20,color=ADMX_col,rotation=0,ha='left',va='top',clip_on=True)
plt.plot([3e-5,2e-5],[1e-15,0.2e-14],'k-',lw=1.5)

INCLUDE_CADEX = False
if (INCLUDE_CADEX):
    #CADEx
    CADEx_col = [0.95,0.5,0.5]
    #CADEx_col = 'k'
    mvals = geomspace(330e-6, 460e-6)
    gvals = (f_ax**-0.5)*1.6e-13*mvals/370e-6
    plt.fill_between(mvals,gvals,y2=y2,edgecolor='k',facecolor=CADEx_col,zorder=0.1, linewidth=3.0, alpha=0.4)
    xvals = concatenate((mvals[0:1], mvals, mvals[-2:-1]))
    yvals = concatenate(([y2], gvals, [y2]))
    plt.plot(xvals, yvals, color='k', zorder=0.11, linewidth=3.0)
    #plt.fill_between(mvals,4e-13*mvals/400e-6,y2=y2,edgecolor='k',facecolor=None,zorder=0.1, linewidth=3.0)
    #plt.plot([330, 460],[dat[0,1],y2],color=CADEx_col,lw=2,zorder=0.05)
    plt.gcf().text(0.764,0.615,r'{\bf CADEx}',color='k',fontsize=30,rotation=90,rotation_mode='anchor',ha='center',va='center')

plt.gca().xaxis.labelpad = 0

plt.gcf().text(0.89,0.13,r'$f_\mathrm{ax} = 90\%; \,\,\rho_0 = 0.45$ GeV cm$^{-3}$',fontsize=25,ha='right', va='bottom')

#plt.title('0',color='w',pad=10)
MySaveFig(fig,'AxionPhoton_RadioFreqCloseup_withProjections', pngsave=False)