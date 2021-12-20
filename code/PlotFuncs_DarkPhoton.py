"""

Originally adapted from https://github.com/cajohare/AxionLimits/
"""

from numpy import *
from numpy.random import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from matplotlib import colors
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from scipy.stats import norm

from Dirs import axionlimits_dir

import os

rootdir = os.path.dirname(os.path.abspath(__file__)) + "/"
pltdir = rootdir + '../plots/'
pltdir_png = rootdir + '../plots/'

def MySaveFig(fig,pltname,pngsave=False):
    fig.savefig(pltdir+pltname+'.pdf',bbox_inches='tight')
    if pngsave:
        fig.savefig(pltdir_png+pltname+'.png',bbox_inches='tight')

def col_alpha(col,alpha=0.1):
    rgb = colors.colorConverter.to_rgb(col)
    bg_rgb = [1,1,1]
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(rgb, bg_rgb)]

from matplotlib import patches
from matplotlib import text as mtext
import numpy as np
import math

class CurvedText(mtext.Text):
    """
    A text object that follows an arbitrary curve.
    """
    def __init__(self, x, y, text, axes, **kwargs):
        super(CurvedText, self).__init__(x[0],y[0],' ', **kwargs)

        axes.add_artist(self)

        ##saving the curve:
        self.__x = x
        self.__y = y
        self.__zorder = self.get_zorder()

        ##creating the text objects
        self.__Characters = []
        for c in text:
            if c == ' ':
                ##make this an invisible 'a':
                t = mtext.Text(0,0,'a')
                t.set_alpha(0.0)
            else:
                t = mtext.Text(0,0,c, **kwargs)

            #resetting unnecessary arguments
            t.set_ha('center')
            t.set_rotation(0)
            t.set_zorder(self.__zorder +1)

            self.__Characters.append((c,t))
            axes.add_artist(t)


    ##overloading some member functions, to assure correct functionality
    ##on update
    def set_zorder(self, zorder):
        super(CurvedText, self).set_zorder(zorder)
        self.__zorder = self.get_zorder()
        for c,t in self.__Characters:
            t.set_zorder(self.__zorder+1)

    def draw(self, renderer, *args, **kwargs):
        """
        Overload of the Text.draw() function. Do not do
        do any drawing, but update the positions and rotation
        angles of self.__Characters.
        """
        self.update_positions(renderer)

    def update_positions(self,renderer):
        """
        Update positions and rotations of the individual text elements.
        """

        #preparations

        ##determining the aspect ratio:
        ##from https://stackoverflow.com/a/42014041/2454357

        ##data limits
        xlim = self.axes.get_xlim()
        ylim = self.axes.get_ylim()
        ## Axis size on figure
        figW, figH = self.axes.get_figure().get_size_inches()
        ## Ratio of display units
        _, _, w, h = self.axes.get_position().bounds
        ##final aspect ratio
        aspect = ((figW * w)/(figH * h))*(ylim[1]-ylim[0])/(xlim[1]-xlim[0])

        #points of the curve in figure coordinates:
        x_fig,y_fig = (
            np.array(l) for l in zip(*self.axes.transData.transform([
            (i,j) for i,j in zip(self.__x,self.__y)
            ]))
        )

        #point distances in figure coordinates
        x_fig_dist = (x_fig[1:]-x_fig[:-1])
        y_fig_dist = (y_fig[1:]-y_fig[:-1])
        r_fig_dist = np.sqrt(x_fig_dist**2+y_fig_dist**2)

        #arc length in figure coordinates
        l_fig = np.insert(np.cumsum(r_fig_dist),0,0)

        #angles in figure coordinates
        rads = np.arctan2((y_fig[1:] - y_fig[:-1]),(x_fig[1:] - x_fig[:-1]))
        degs = np.rad2deg(rads)


        rel_pos = 10
        for c,t in self.__Characters:
            #finding the width of c:
            t.set_rotation(0)
            t.set_va('center')
            bbox1  = t.get_window_extent(renderer=renderer)
            w = bbox1.width
            h = bbox1.height

            #ignore all letters that don't fit:
            if rel_pos+w/2 > l_fig[-1]:
                t.set_alpha(0.0)
                rel_pos += w
                continue

            elif c != ' ':
                t.set_alpha(1.0)

            #finding the two data points between which the horizontal
            #center point of the character will be situated
            #left and right indices:
            il = np.where(rel_pos+w/2 >= l_fig)[0][-1]
            ir = np.where(rel_pos+w/2 <= l_fig)[0][0]

            #if we exactly hit a data point:
            if ir == il:
                ir += 1

            #how much of the letter width was needed to find il:
            used = l_fig[il]-rel_pos
            rel_pos = l_fig[il]

            #relative distance between il and ir where the center
            #of the character will be
            fraction = (w/2-used)/r_fig_dist[il]

            ##setting the character position in data coordinates:
            ##interpolate between the two points:
            x = self.__x[il]+fraction*(self.__x[ir]-self.__x[il])
            y = self.__y[il]+fraction*(self.__y[ir]-self.__y[il])

            #getting the offset when setting correct vertical alignment
            #in data coordinates
            t.set_va(self.get_va())
            bbox2  = t.get_window_extent(renderer=renderer)

            bbox1d = self.axes.transData.inverted().transform(bbox1)
            bbox2d = self.axes.transData.inverted().transform(bbox2)
            dr = np.array(bbox2d[0]-bbox1d[0])

            #the rotation/stretch matrix
            rad = rads[il]
            rot_mat = np.array([
                [math.cos(rad), math.sin(rad)*aspect],
                [-math.sin(rad)/aspect, math.cos(rad)]
            ])

            ##computing the offset vector of the rotated character
            drp = np.dot(dr,rot_mat)

            #setting final position and rotation:
            t.set_position(np.array([x,y])+drp)
            t.set_rotation(degs[il])

            t.set_va('center')
            t.set_ha('center')

            #updating rel_pos to right edge of character
            rel_pos += w-used

def FigSetup(xlab=r'Dark photon mass $m_{\gamma^\prime}$ [eV]',ylab='Kinetic mixing $\epsilon$',\
             chi_min = 1.0e-18,chi_max = 1.0e0,\
             m_min = 1.0e-15,m_max = 1e5,\
             lw=2.5,lfs=40,tfs=25,tickdir='out',\
             Grid=False,Shape='Rectangular',mathpazo=True,\
             TopAndRightTicks=False,FrequencyAxis=True,FrequencyLabels=True,UnitAxis=True,f_rescale=1,\
            tick_rotation = 20,width=20,height=10,upper_tickdir='out', upper_xlabel=None):

    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)
    

    
    
    if mathpazo:
        mpl.rcParams['text.latex.preamble'] = [r'\usepackage{mathpazo}']

    if Shape=='Wide':
        fig = plt.figure(figsize=(16.5,5))
    elif Shape=='Rectangular':
        fig = plt.figure(figsize=(16.5,11))
    elif Shape=='Square':
        fig = plt.figure(figsize=(14.2,14))
    elif Shape=='Custom':
        fig = plt.figure(figsize=(width,height))

    ax = fig.add_subplot(111)

    ax.set_xlabel(xlab,fontsize=lfs)
    ax.set_ylabel(ylab,fontsize=lfs)
    
    #ax.plot(m_dp,chi100,color="black",zorder=100)
    #ax.plot(m_dp,chi10,color="blue",zorder=100)
    #ax.plot(m_dp,chi1,color="orange",zorder=100)
    
    
    ax.tick_params(which='major',direction=tickdir,width=2.5,length=13,right=TopAndRightTicks,top=TopAndRightTicks,pad=7, labelsize=tfs)
    ax.tick_params(which='minor',direction=tickdir,width=1,length=10,right=TopAndRightTicks,top=TopAndRightTicks, labelsize=tfs)


    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim([m_min,m_max])
    ax.set_ylim([chi_min,chi_max])

    locmaj = mpl.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=50)
    locmin = mpl.ticker.LogLocator(base=10.0, subs=arange(2, 10)*.1,numticks=100)
    ax.xaxis.set_major_locator(locmaj)
    ax.xaxis.set_minor_locator(locmin)
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

    locmaj = mpl.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
    locmin = mpl.ticker.LogLocator(base=10.0, subs=arange(2, 10)*.1,numticks=100)
    ax.yaxis.set_major_locator(locmaj)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

    if Shape=='Rectangular':
        plt.xticks(rotation=tick_rotation)

    if Grid:
        ax.grid(zorder=0)

    if FrequencyAxis:
        ax2 = ax.twiny()



        ax2.set_xscale('log')
        ax2.tick_params(which='major',direction=upper_tickdir,width=2.5,length=13,pad=7, labelsize=tfs)
        ax2.tick_params(which='minor',direction=upper_tickdir,width=1,length=10, labelsize=tfs)
        locmaj = mpl.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=50)
        locmin = mpl.ticker.LogLocator(base=10.0, subs=arange(2, 10)*.1,numticks=100)
        ax2.xaxis.set_major_locator(locmaj)
        ax2.xaxis.set_minor_locator(locmin)
        ax2.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

        if FrequencyLabels:
            ax2.set_xticks([1e0,1e3,1e6,1e9,1e12,1*241.8*1e12,1000*241.8*1e12])
            ax2.set_xticklabels(['Hz','kHz','MHz','GHz','THz','eV','keV'])
        ax2.set_xlim([m_min*241.8*1e12/f_rescale,m_max*241.8*1e12/f_rescale])
        
        if (upper_xlabel is not None):
            ax2.set_xlabel(upper_xlabel, color='white', fontsize=lfs)

        plt.sca(ax)
    return fig,ax
    
    
import matplotlib.patheffects as pe
pek=[pe.Stroke(linewidth=7, foreground='k'), pe.Normal()]

    
def Haloscopes(ax,col=[0.75, 0.2, 0.2],fs=17,projection=True,text_on=True):
    y2 = ax.get_ylim()[1]
    zo = 0.3
    
    HAYSTAC_col = 'indianred'
    CAPP_col = 'crimson'
    QUAX_col = 'r'
    ADMX_col = 'firebrick'
    
    # ADMX
    costh = sqrt(0.0025)
    B = 7.6
    dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX.txt")
    dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=ADMX_col,zorder=0.1,lw=3)

    B = 6.8
    dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX2018.txt")
    dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=ADMX_col,zorder=0.1)

    B = 7.6
    dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX2019_1.txt")
    dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=ADMX_col,zorder=0.1)

    B = 7.6
    dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX2019_2.txt")
    dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=ADMX_col,zorder=0.1)

#     B = 3.11
#     dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX_Sidecar_AC.txt")
#     dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
#     plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=ADMX_col,facecolor=ADMX_col,zorder=0.1)

#     B = 5.0
#     dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/ADMX_SLIC.txt")
#     dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
#     plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=ADMX_col,facecolor=ADMX_col,zorder=100)

    
    
    B = 9
    dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/HAYSTAC_highres.txt")
    dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=HAYSTAC_col,zorder=0.1)
    dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/HAYSTAC_2020_highres.txt")
    dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=HAYSTAC_col,zorder=0.1)

    
    # CAPP
    B = 7.3
    dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/CAPP-1.txt")
    dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=CAPP_col,zorder=0.1)

    B = 7.8
    dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/CAPP-2.txt")
    dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=CAPP_col,zorder=0.1)

    B = 7.9
    dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/CAPP-3.txt")
    dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*costh*dat[:,0]))
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor='none',facecolor=CAPP_col,zorder=0.1)

    # CAPP-3 [KSVZ]
    dat_min = dat[argmin(dat[:,1]),:]
    dat_min[1] = dat_min[1]*costh/0.11
    plt.plot([dat_min[0],dat_min[0]],[1e-10,dat_min[1]],'-',color=CAPP_col,lw=1.5,zorder=0.1)


    B = 8.1
    dat = loadtxt(axionlimits_dir + "limit_data/AxionPhoton/QUAX.txt")
    dat[:,1] = 1e-9*dat[:,1]*(B/(1.444e-3*0.0046*dat[:,0]))
    plt.fill_between([dat[0,0],dat[0,0]],[y2,dat[0,1]],y2=y2,color=QUAX_col,zorder=0.1)


    if text_on: 
        plt.text(1.9e-6,0.5e-14,r'{\bf ADMX}',fontsize=fs,color=ADMX_col,rotation=90,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.8e-5,0.4e-13,r'{\bf CAPP}',fontsize=fs-2,color=CAPP_col,rotation=90,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.19e-4,8e-15,r'{\bf HAYSTAC}',fontsize=fs-5,color=HAYSTAC_col,rotation=90,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.55e-4,4e-11,r'{\bf QUAX}',fontsize=fs-5,color=QUAX_col,rotation=-90,rotation_mode='anchor',ha='center',va='center')

    return
    
    
def StellarBounds(ax, fs=20,text_on=True, Higgsed=True, e_D = 0.1):
    y2 = ax.get_ylim()[1]
    # Stellar physics constraints

    # Globular clusters 
    HB_col = [0.01, 0.75, 0.24]
    HB = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/RG.txt")

    if (Higgsed):
        HB = np.concatenate((np.atleast_2d([1e-3, 1.1e-13*(0.1/e_D)]), np.atleast_2d([HB[0,0], 1.1e-13*(0.1/e_D)]), HB))
                            
        #mask = (HB[:,0] < 100) & (HB[:,1] >  1.1e-13)
        #HB[mask,1] = HB[mask,1]*0 + 1.1e-13
    plt.plot(HB[:,0],HB[:,1],color='k',alpha=0.5,zorder=0.9,lw=2)
    plt.fill_between(HB[:,0],HB[:,1],y2=y2,edgecolor=None,facecolor=HB_col,zorder=0.9)
    
    # Globular clusters 
    HB_col = 'DarkGreen'
    HB = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/HB.txt")
    if (Higgsed):
        HB = np.concatenate((np.atleast_2d([1e-3, 2e-13*(0.1/e_D)]), np.atleast_2d([HB[0,0], 2e-13*(0.1/e_D)]), HB))
        #mask = (HB[:,0] < 100) & (HB[:,1] >  2e-13)
        #HB[mask,1] = HB[mask,1]*0 + 2e-13
    plt.plot(HB[:,0],HB[:,1],color='k',alpha=0.5,zorder=0.95,lw=2)
    plt.fill_between(HB[:,0],HB[:,1],y2=y2,edgecolor=None,facecolor=HB_col,zorder=0.95)

    # Solar bound
    Solar_col = 'ForestGreen'
    Solar = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Solar.txt")
    if (Higgsed):
        mask = (Solar[:,0] < 100) & (Solar[:,1] > 3.4e-13*(0.1/e_D))
        Solar[mask,1] = Solar[mask,1]*0 + 3.4e-13*(0.1/e_D)
    plt.plot(Solar[:,0],Solar[:,1],color='k',alpha=0.5,zorder=1,lw=2)
    plt.fill_between(Solar[:,0],Solar[:,1],y2=y2,edgecolor=None,facecolor=Solar_col,zorder=1)

    if not Higgsed:
        Solar = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Solar-Global.txt")
        plt.plot(Solar[:,0],Solar[:,1]/Solar[:,0],color='k',alpha=0.5,zorder=1.01,lw=2)
        plt.fill_between(Solar[:,0],Solar[:,1]/Solar[:,0],y2=y2,edgecolor=None,facecolor=Solar_col,zorder=1.01)

    if text_on:
        plt.text(1e2*(1-0.01),1.5e-14*(1+0.05),r'{\bf Solar}',fontsize=fs,color='k',rotation=-39,rotation_mode='anchor',ha='center',va='center')
        plt.text(1e3*(1-0.01),0.7e-14*(1+0.05),r'{\bf HB}',fontsize=fs,color='k',rotation=-36,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.8e4*(1-0.01),0.7e-14*(1+0.05),r'{\bf RG}',fontsize=fs,color='k',rotation=-35,rotation_mode='anchor',ha='center',va='center')
        plt.text(1e2,1.5e-14,r'{\bf Solar}',fontsize=fs,color='w',rotation=-39,rotation_mode='anchor',ha='center',va='center')
        plt.text(1e3,0.7e-14,r'{\bf HB}',fontsize=fs,color='w',rotation=-36,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.8e4,0.7e-14,r'{\bf RG}',fontsize=fs,color='w',rotation=-35,rotation_mode='anchor',ha='center',va='center')
    return

    
def Xenon(ax,col='crimson',fs=23,text_on=True, f_DP = 1.0):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Xenon1T.txt")
    dat[:,1] = dat[:,1]*sqrt(0.3/(0.45*f_DP))

    plt.plot(1e3*dat[:,0],dat[:,1],color='k',alpha=0.5,zorder=0.5,lw=2)
    plt.fill_between(1e3*dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5)
    
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Xenon1T_S1S2.txt")
    dat[:,1] = dat[:,1]*sqrt(0.3/(0.45*f_DP))

    plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.5,zorder=0.5,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5)
    
    
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Xenon100.txt")
    dat[:,1] = dat[:,1]*sqrt(0.3/(0.45*f_DP))

    plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.5,zorder=0.5,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.5)
    if text_on: 
        plt.text(8e2,3e-17,r'{\bf XENON}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
    
    return





def DAMIC(ax,col='salmon',fs=21,text_on=True, f_DP = 1.0):
    m1,y1 = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/DM_combined.txt",unpack=True)
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/DAMIC.txt")
    dat[:,1] = dat[:,1]*sqrt(0.3/(0.45*f_DP))

    y2 = interp(dat[:,0],m1,y1)
    dat[0,1] = y2[0]
    dat[-1,1] = y2[-1]
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.05,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.05)
    if text_on: 
        plt.text(6e-1,1.5e-14,r'{\bf DAMIC}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.plot([5e0,1e1],[3e-14,6e-14],'k-',lw=2.5,color=col)
    return


def FUNK(ax,col='red',fs=21,text_on=True, f_DP = 1.0):
    m1,y1 = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/DM_combined.txt",unpack=True)
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/FUNK.txt")
    dat[:,1] = dat[:,1]*sqrt(0.3/(0.45*f_DP))*sqrt(2/3/0.27)

    y2 = interp(dat[:,0],m1,y1)
    dat[0,1] = y2[0]/1.1
    dat[-1,1] = y2[-1]/1.1
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.075,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.075)
    if text_on: 
        plt.text(5e-1,1e-13,r'{\bf FUNK}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.plot([9e-1,3e0],[3e-13,1e-12],'k-',lw=2.5,color=col)
    return

def SENSEI(ax,col='firebrick',fs=21,text_on=True, f_DP = 1.0):
    #y2 = ax.get_ylim()[1]
    y2 = 1e-10
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/SENSEI.txt")
    dat[:,1] = dat[:,1]*sqrt(0.3/(0.45*f_DP))

    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.025,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.025)
    if text_on: 
        plt.text(1e0,3e-14,r'{\bf SENSEI}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.plot([2.5e0,4e0],[3e-14,5.3e-14],'k-',lw=2.5,color=col)
    return


def Nanowire(ax,col='pink',fs=22,text_on=True):
    m1,y1 = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/DM_combined.txt",unpack=True)
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/WSi_Nanowire.txt")
    dat[:,1] = dat[:,1]*sqrt(0.3/0.45)
    y2 = interp(dat[:,0],m1,y1)
    dat[0,1] = y2[0]/1.1
    dat[-1,1] = y2[-1]/1.1
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.3,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.3)
    if text_on: 
        plt.text(5e-4,1e-10,r'{\bf WSi Nanowire}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.plot([9e-3,3e-3],[3e-10,9e-10],'k-',lw=2.5,color=col)
    return



def Tokyo(ax,col='darkred',fs=15,text_on=True):
    m1,y1 = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/DM_combined.txt",unpack=True)
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Tokyo-Dish.txt")
    dat[:,1] = dat[:,1]*sqrt(2/3/0.5)
    y2 = interp(dat[:,0],m1,y1)
    dat[0,1] = y2[0]/1.1
    dat[-1,1] = y2[-1]/1.1
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.4,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.4)

    
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Tokyo-Knirck.txt")
    dat[:,1] = dat[:,1]*sqrt(1/3/0.17)
    plt.fill_between(dat[:,0],dat[:,1],y2=1e0,edgecolor='k',facecolor=col,zorder=1.09)
    
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Tokyo-Tomita.txt")
    plt.plot([dat[1,0],dat[1,0]],[dat[1,1],1e0],'-',color=col,lw=3,zorder=0.2)
    if text_on: 
        plt.text(2e-4,2.5e-10,r'{\bf Tokyo-3}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.5e-3,3e-8,r'{\bf Tokyo-2}',fontsize=fs-2,color='k',rotation=90,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.45e-1,4e-12,r'{\bf Tokyo-1}',fontsize=fs+4,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.plot([3e-1,4e0],[5e-12,8e-12],'-',lw=2.5,color=col)
    return
    

def Jupiter(ax,col='Green',fs=17,text_on=True):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Jupiter.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=2,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=2)
    if text_on: 
        plt.text(1e-14*(1-0.02),3e-1*(1+0.07),r'{\bf Jupiter}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.text(1e-14,3e-1,r'{\bf Jupiter}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')
    return

def Earth(ax,col='DarkGreen',fs=17,text_on=True):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Earth.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.9,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.9)
    if text_on: 
        plt.text(4e-13*(1-0.01),2e-1*(1+0.05),r'{\bf Earth}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.text(4e-13,2e-1,r'{\bf Earth}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')
    return


def Crab(ax,col=[0.1,0.4,0.1],fs=17,text_on=True):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Crab.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.09999,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.09999)
    
#     dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Crab_2.txt")
#     plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.9,lw=2)
#     plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.9)
    if text_on: 
        plt.text(0.5e-6*(1-0.02),3e-1*(1+0.07),r'{\bf Crab}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.5e-6,3e-1,r'{\bf Crab}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')
    
        plt.text(0.8e-6*(1-0.02),0.9e-1*(1+0.07),r'{\bf nebula}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.8e-6,0.9e-1,r'{\bf nebula}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')
    
    return


def SHUKET(ax,col='maroon',fs=13,text_on=True,edge_on=False,lw=0.8):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/SHUKET.txt")
    dat[:,1] = dat[:,1]*sqrt(0.3/0.45)*sqrt(1/3/0.0086)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.2)
    if edge_on:
        plt.plot(dat[:,0],dat[:,1],'k-',lw=lw,zorder=0.2)
    if text_on: 
        plt.text(3.5e-5,1.5e-12,r'{\bf SHUKET}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center')
    return

def DarkEfield(ax,col='darkred',fs=17,text_on=True,edge_on=False,lw=0.8):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/DarkEfield.txt")
    dat[:,1] = dat[:,1]*sqrt(0.3/0.45)*sqrt(1/3/0.079)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.2)
    if edge_on:
        plt.plot(dat[:,0],dat[:,1],'k-',lw=lw,zorder=0.2)
    if text_on: 
        plt.text(0.88e-7,0.2e-12,r'{\bf Dark}',fontsize=fs,color=col,rotation=90,rotation_mode='anchor',ha='center',va='center')
        plt.text(2e-7,0.2e-12,r'{\bf E-field}',fontsize=fs,color=col,rotation=90,rotation_mode='anchor',ha='center',va='center')
    return

def WISPDMX(ax,col='crimson',fs=12,text_on=True,edge_on=False,lw=0.8):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/WISPDMX.txt")
    dat[:,1] = dat[:,1]*sqrt(0.3/0.45)*sqrt(1/3/0.079)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.201)
    if edge_on:
        plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=0.201,lw=lw)

    if text_on: 
        plt.text(9e-7,4.1e-12,r'{\bf WISP}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.text(9e-7,1.8e-12,r'{\bf DMX}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')

    return

def SQuAD(ax,col=[0.7,0,0],fs=12,text_on=True,lw=0.5,point_on=False,ms=10):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/SQuAD.txt")
    dat[:,1] = dat[:,1]*sqrt(0.4/0.45)*sqrt(1/3/0.0025)
    plt.plot([dat[0,0],dat[0,0]],[y2,dat[0,1]],lw=lw,color=col,alpha=1,zorder=0.2)
    if point_on:
        plt.plot(dat[0,0],dat[0,1],'o',mfc=col,mec='k',mew=lw+1,zorder=0.2,markersize=ms)
    if text_on: 
        plt.text(36e-6,0.25e-13,r'{\bf SQuAD}',fontsize=fs,color=col,rotation=-90,rotation_mode='anchor',ha='center',va='center')
    return

def DMPathfinder(ax,col='pink',fs=13,text_on=True):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/DM-Pathfinder.txt")
    dat[:,1] = dat[:,1]*sqrt(1/0.028)
    plt.plot([dat[0,0],dat[0,0]],[y2,dat[0,1]],lw=2,color=col,alpha=1,zorder=0.6)
    if text_on: 
        plt.text(2.1e-9,0.5e-8,r'{\bf DM}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.text(2.1e-9,0.2e-8,r'{\bf Pathfinder}',fontsize=fs,color=col,rotation=0,rotation_mode='anchor',ha='center',va='center')

    return

def DarkMatter(ax,Witte_col='royalblue',Caputo_col='dodgerblue',Arias_col='navy',fs=20,projection=True,text_on=True):
    y2 = ax.get_ylim()[1]
    zo = 0.3
    
    # Combined limits
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/DM_combined.txt")
    plt.plot(dat[:,0],dat[:,1],'-',color='w',alpha=1,zorder=zo+0.1,lw=2.5,path_effects=pek)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor='lightgray',zorder=zo,alpha=1.0)
    plt.plot([1e-16,dat[0,0]],[dat[0,1],dat[0,1]],'--',color='w',alpha=1,zorder=zo+0.1,lw=2.5,path_effects=pek)
    plt.fill_between([1e-16,dat[0,0]],[dat[0,1],dat[0,1]],y2=y2,edgecolor=None,facecolor='lightgray',zorder=zo+0.1,alpha=1.0)
    plt.plot(dat[40:,0],dat[40:,1],'--',color='w',alpha=1,lw=2.5,zorder=1000,solid_capstyle='round')
    
    # Individual limits
    #dat2 = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Cosmology_Witte_inhomogeneous.txt")
    #dat4 = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Cosmology_Caputo_HeII.txt",delimiter=',')
    #dat5 = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Cosmology_Arias.txt")
    
    #plt.fill_between(dat2[:,0],dat2[:,1],y2=y2,edgecolor='k',facecolor=Witte_col,zorder=0.305,alpha=0.8)
    #plt.fill_between(dat4[:,0],dat4[:,1],y2=y2,edgecolor='k',facecolor=Caputo_col,zorder=0.305,alpha=0.8)
    #plt.fill_between(dat5[:,0],dat5[:,1],y2=y2,edgecolor='k',facecolor=Arias_col,zorder=0.306,alpha=1)

    """
    if text_on: 
        plt.gcf().text(0.21,0.42,r'{\bf DPDM} HeII',fontsize=15,color='w',ha='center')
        plt.gcf().text(0.21,0.4,r'Reionisation',fontsize=15,color='w',ha='center')
        plt.gcf().text(0.21,0.38,r'(Caputo et al.)',fontsize=13,color='w',ha='center')

        plt.gcf().text(0.297,0.37,r'{\bf DPDM}',fontsize=17,color='w',ha='center')
        plt.gcf().text(0.297,0.35,r'(Witte et al.)',fontsize=13,color='w',ha='center')

        plt.gcf().text(0.44,0.48,r'{\bf DPDM}',fontsize=18,color='w',ha='center')
        plt.gcf().text(0.44,0.46,r'(Arias et al.)',fontsize=16,color='w',ha='center')

        plt.arrow(0.24, 0.26, 0, -0.045, transform=fig.transFigure,figure=fig,
          length_includes_head=True,lw=2.5,
          head_width=0.012, head_length=0.028, overhang=0.13,
          edgecolor='k',facecolor='w',clip_on=False)
        plt.text(1e-12,0.4e-16,r'{\bf Dark}',fontsize=27,ha='center')
        plt.text(1e-12,0.7e-17,r'{\bf photon DM}',fontsize=27,ha='center')
    """
    return
    
def COBEFIRAS(ax,col=[0.1,0.2,0.5],text_on=True):
    y2 = ax.get_ylim()[1]   
    dat3 = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/COBEFIRAS.txt",delimiter=',')
    plt.fill_between(dat3[:,0],dat3[:,1],y2=y2,edgecolor='k',facecolor=col,zorder=0.5,alpha=1)
    if text_on: 
        plt.gcf().text(0.21,0.72,r'{\bf COBE/FIRAS}',fontsize=22,color='w',ha='center')
        plt.gcf().text(0.21,0.69,r'$\gamma \rightarrow X$',fontsize=22,color='w',ha='center')
    return


def IGM(ax,col='seagreen',fs=17,text_on=True):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/IGM.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.5,zorder=0.49,lw=2)

    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.49)
    
    if text_on: 
        plt.text(1.2e-12*(1-0.05),0.125e-7*(1+0.07),r'{\bf DPDM} (IGM)',fontsize=fs,color='k',rotation=-38,rotation_mode='anchor',ha='center',va='center')
        plt.text(1.2e-12,0.125e-7,r'{\bf DPDM} (IGM)',fontsize=fs,color='w',rotation=-38,rotation_mode='anchor',ha='center',va='center')
    return

def LeoT(ax,col='mediumseagreen',fs=19,text_on=True):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/LeoT.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=0.48,zorder=0.306,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=0.306)
    
    if text_on: 
        plt.text(1.1e-13*(1-0.07),1.29e-9*(1+0.1),r'{\bf DPDM} (Leo T)',fontsize=fs,color='k',rotation=-37,rotation_mode='anchor',ha='center',va='center')
        plt.text(1.1e-13,1.29e-9,r'{\bf DPDM} (Leo T)',fontsize=fs,color='w',rotation=-37,rotation_mode='anchor',ha='center',va='center')
    return


def LSW(ax,text_on=True):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/SPring-8.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.1001,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.45, 0.05, 0.1],zorder=1.1001)
    
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/ALPS.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.091,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.55, 0.0, 0.16],zorder=1.091)
    
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/LSW_UWA.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.09,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.6, 0.0, 0.2],zorder=1.09)
    
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/LSW_ADMX.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.089,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.65, 0.1, 0.24],zorder=1.089)
    
#     dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/LSW_CERN.txt")
#     plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.089,lw=2)
#     plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.65, 0.15, 0.2],zorder=1.089)
    
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/CROWS.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.08,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.7, 0.2, 0.2],zorder=1.08)
    
    if text_on: 
        plt.text(0.4e-6,0.15e-3,r'{\bf LSW-ADMX}',fontsize=17,color='w',rotation=-55,rotation_mode='anchor',ha='center',va='center')
        plt.text(1e-5,5e-5,r'{\bf LSW-UWA}',fontsize=14,color='w',rotation=-55,rotation_mode='anchor',ha='center',va='center')
        
        plt.text(0.6e0*(1-0.02),0.9e-4*(1+0.08),r'{\bf LSW-SPring-8}',fontsize=15,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.6e0,0.9e-4,r'{\bf LSW-SPring-8}',fontsize=15,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')

        
        plt.text(1.2e-4*(1-0.02),0.9e-5*(1+0.08),r'{\bf ALPS}',fontsize=25,color='k',rotation=-54,rotation_mode='anchor',ha='center',va='center')
        plt.text(1.2e-4,0.9e-5,r'{\bf ALPS}',fontsize=25,color='w',rotation=-54,rotation_mode='anchor',ha='center',va='center')

        plt.text(0.75e-7*(1-0.01),9.9e-5*(1+0.05),r'{\bf CROWS}',fontsize=24,color='k',rotation=-54,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.75e-7,9.9e-5,r'{\bf CROWS}',fontsize=24,color='w',rotation=-54,rotation_mode='anchor',ha='center',va='center')
    return

def Coulomb(ax,text_on=True):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Cavendish.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.07,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.7,0,0],zorder=1.07)
    
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/PlimptonLawton.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.071,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor='crimson',zorder=1.071)
    
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/Spectroscopy.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.11,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.4, 0.0, 0.13],zorder=1.11)
    
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/AFM.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.5,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=[0.4, 0.2, 0.2],zorder=1.5)
    if text_on: 
        plt.text(2e-10*(1-0.02),0.3e-1*(1+0.08),r'{\bf Plimpton-Lawton}',fontsize=18,color='k',rotation=-35,rotation_mode='anchor',ha='center',va='center')
        plt.text(2e-10,0.3e-1,r'{\bf Plimpton-Lawton}',fontsize=18,color='w',rotation=-35,rotation_mode='anchor',ha='center',va='center')
        
        plt.text(3e1*(1-0.02),3e-1*(1+0.08),r'{\bf AFM}',fontsize=20,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.text(3e1,3e-1,r'{\bf AFM}',fontsize=20,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')

        plt.text(0.7e-8*(1-0.02),3e-6*(1+0.08),r'{\bf Cavendish-Coulomb}',fontsize=23,color='k',rotation=-35,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.7e-8,3e-6,r'{\bf Cavendish-Coulomb}',fontsize=23,color='w',rotation=-35,rotation_mode='anchor',ha='center',va='center')
        
        plt.text(0.2e2*(1-0.01),1e-3*(1+0.08),r'{\bf Spectroscopy}',fontsize=23,color='k',rotation=-30,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.2e2,1e-3,r'{\bf Spectroscopy}',fontsize=23,color='w',rotation=-30,rotation_mode='anchor',ha='center',va='center')
    
    return


def CAST(ax,col='maroon',fs=27,text_on=True):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/CAST.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.1,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.1)
    if text_on: 
        plt.text(1e3*(1-0.01),0.8e-6*(1+0.08),r'{\bf CAST}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.text(1e3,0.8e-6,r'{\bf CAST}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')
    return

def SHIPS(ax,col='indianred',fs=20,text_on=True):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/SHIPS.txt")
    dat[:,1] = dat[:,1]/dat[:,0]
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.09,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.09)
    if text_on: 
        plt.text(0.6e-1*(1-0.05),0.08e-8*(1+0.1),r'{\bf SHIPS}',fontsize=fs,color='k',rotation=-30,rotation_mode='anchor',ha='center',va='center')
        plt.text(0.6e-1,0.08e-8,r'{\bf SHIPS}',fontsize=fs,color='w',rotation=-30,rotation_mode='anchor',ha='center',va='center')
    return

def TEXONO(ax,col=[0.5, 0.0, 0.13],fs=15,text_on=True):
    y2 = ax.get_ylim()[1]
    dat = loadtxt(axionlimits_dir + "limit_data/DarkPhoton/TEXONO.txt")
    plt.plot(dat[:,0],dat[:,1],color='k',alpha=1,zorder=1.101,lw=2)
    plt.fill_between(dat[:,0],dat[:,1],y2=y2,edgecolor=None,facecolor=col,zorder=1.101)
    if text_on: 
        plt.text(1.2e2*(1-0.01),0.1e-4*(1+0.08),r'{\bf TEXONO}',fontsize=fs,color='k',rotation=0,rotation_mode='anchor',ha='center',va='center')
        plt.text(1.2e2,0.1e-4,r'{\bf TEXONO}',fontsize=fs,color='w',rotation=0,rotation_mode='anchor',ha='center',va='center')
    return