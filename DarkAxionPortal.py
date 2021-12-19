import numpy as np
from Units import *

from scipy import special, interpolate, optimize

SUPPRESS_WARNINGS = True

if (SUPPRESS_WARNINGS):
    np.seterr(divide='ignore', invalid='ignore')
    print("> DarkAxionPortal.py: Suppressing divide-by-zero errors!")


#-------------------
#Specificiation of some useful constants
h           = 0.674   #https://arxiv.org/abs/1807.06209
s0          = 2889.2*cm**-3 #Entropy density today
rho_crit    = 1.05368e-5*h**2*GeV*cm**-3    #Critical density today
Omegah2_dm  = 0.12 #DM density Omega_DM*h^2
Omega_dm    = Omegah2_dm/h**2

Mpl_reduced = 2.435e18*GeV #Reduced Planck Mass
N_C         = 3 #Number of colors
alpha_em    = 0.0072973525693   #EM fine structure constant
alpha_s     = 0.12 #g_s^2 = 4*pi*alpha_s
g_s         = np.sqrt(4*np.pi*alpha_s) #Strong coupling constant

e           = np.sqrt(4*np.pi*alpha_em)   #Electron charge
m_e         = 510.99895*keV               #Electron mass

z_q         = 0.56 #z = m_u/m_d
f_pi        = 92*MeV #Pion decay constant
m_pi        = 135*MeV #Pion mass

#--------------------

#Anharmonic correction for axion energy density
#See Fig. 4 of https://journals.aps.org/prd/abstract/10.1103/PhysRevD.33.889
#Could be a correction of ~10
def load_anharmonic_correction():
    theta, F = np.loadtxt("data/AnharmonicCorrection.txt", unpack=True)
    F[F > 10] = 0.0*F[F > 10] + 10.0 #Truncate at 10 (because it tends to infinity as theta -> pi)
    return interpolate.interp1d(theta, F, bounds_error=False, fill_value=(F[0], F[-1]))
    
anharmonic_correction = load_anharmonic_correction()


class Model():
    """Class to represent a Dark Axion portal model with fixed couplings.
        
       The model is specified by the charges/couplings PQ_Phi, e_D and D_psi.
       This leaves the Dark Photon mass (m_dp), the axion decay constant (f_a,
       or equivalently axion mass, m_a), and kinetic mixing (epsilon) as free
       parameters.
    """
    
    def __init__(self, PQ_Phi = 1.0, e_D = 0.1, D_psi = 3):
        """Constructs a Dark Axion Portal Model with specified charges and couplings.
        
            Arguments:
            ---------
                PQ_Phi (float): PQ Charge of the Phi field (also the PQ color anomaly) (default: 1.0)
                e_D (float): Dark Coupling constant for U(1)_Dark (analogous to e for U(1)_EM )
                D_psi (float): Dark Charge of the psi field
   
        """
        
        #Specification of charges/couplings
        self.PQ_Phi  = PQ_Phi
        self.e_D     = e_D
        self.D_psi   = D_psi
        self.Q_psi   = 0
        self.N_ax    = 1.0*PQ_Phi

    #Couplings
    def calc_G_agg(self, f_a):
        """Calculate (dimensionful) axion-gluon-gluon coupling.
        
            Arguments:
            ---------
                f_a (float): Value of the axion decay constant
        
            Returns:
            ---------
                G_agg (float): axion-gluon-gluon coupling
        """
        return ((g_s**2)/(8*np.pi**2))*(self.PQ_Phi/f_a)

    def calc_G_app(self, f_a):
        """Calculate (dimensionful) axion-photon-photon coupling.
        
            Arguments:
            ---------
                f_a (float): Value of the axion decay constant
        
            Returns:
            ---------
                G_app (float): axion-photon-photon coupling.
        """
        return (e**2/(8*np.pi**2))*(self.PQ_Phi/f_a)*(2*N_C*Q_psi**2 - (2/3)*(4+z_q)/(1+z_q))

    def calc_G_apdp(self, f_a, epsilon):
        """Calculate (dimensionful) axion-photon-dark photon coupling.
        
            Arguments:
            ---------
                f_a (float): Value of the axion decay constant
        
            Returns:
            ---------
                G_app (float): axion-photon-dark photon coupling.
        """
        G_app = self.calc_G_app(f_a)
        return ((e*self.e_D)/(8*np.pi**2))*(self.PQ_Phi/f_a)*(2*N_C*self.D_psi*self.Q_psi) + epsilon*G_app
    
    def calc_G_adpdp(self, f_a, epsilon):
        """Calculate (dimensionful) axion-dark photon-dark photon coupling.
        
            Arguments:
            ---------
                f_a (float): Value of the axion decay constant
        
            Returns:
            ---------
                G_app (float): axion-dark photon-dark photon coupling.
        """
        G_apdp = self.calc_G_apdp(f_a, epsilon)
        return ((self.e_D**2)/(8*np.pi**2))*(self.PQ_Phi/f_a)*(2*N_C*self.D_psi**2) + 2*epsilon*G_apdp

    #Axion calculations
    #-------------------------------------
    def fa_to_ma(self, f_a):
        """Converts axion decay constant to axion mass.
        
            Arguments:
            ---------
                f_a (float): Axion decay constant
        
            Returns:
            ---------
                m_a (float): Axion mass
        """
        return self.N_ax*((np.sqrt(z_q)/(1+z_q))*f_pi*m_pi/f_a)
    
    def ma_to_fa(self, m_a):
        """Converts axion mass to axion decay constant
        
            Arguments:
            ---------
                m_a (float): Axion decay constant
        
            Returns:
            ---------
                f_a (float): Axion mass
        """
        return self.N_ax*((np.sqrt(z_q)/(1+z_q))*f_pi*m_pi/m_a)

    def calc_Omegah2_Axion(self, f_a, theta_i = np.pi/np.sqrt(3), gamma = 1, H_I = 0*GeV):
        """Calculate Axion density Omega_a*h^2 in the pre-inflationary scenario
        
            Arguments:
            ---------
                fa (float): Value of the axion decay constant
                theta_i (float): Initial misalignment angle (default: pi/sqrt(3))
                gamma (float): Dilution factor (default: 1.0)
                H_I (float): Energy scale of inflation (default: 0.0)
        
            Returns:
            ---------
                f_ax (float): Fraction of DM density in axions
        """
        f_corr = anharmonic_correction(theta_i*self.N_ax)
        return 0.43*(f_a/(self.N_ax*1e12*GeV))**(7/6)*(theta_i**2 + (H_I/(2*np.pi*f_a))**2)*f_corr*gamma

    def calc_fa_AxionDensity(self, f_ax, theta_i = np.pi/np.sqrt(3), gamma = 1):
        """Calculate axion decay constant f_a corresponding to a particular value of f_ax = Omega_axion/Omega_DM,
           in the pre-inflationary scenario, ignoring backreaction contribution.
        
            Arguments:
            ---------
                f_ax (float): Fraction of DM density in axions
                theta_i (float): Initial misalignment angle (default: pi/sqrt(3))
                gamma (float): Dilution factor (default: 1.0)
        
            Returns:
            ---------
                f_a (float): Value of axion decay constant.
        """
        #From Eq. (7) of https://arxiv.org/abs/hep-th/0409059
        Omegah2_a = f_ax*Omegah2_dm
        f_corr = anharmonic_correction(theta_i*self.N_ax)
        return 1e12*GeV*self.N_ax*(Omegah2_a/(0.43*theta_i**2*gamma*f_corr))**(6/7)

    def calc_fa_AxionDensity_full(self, f_ax, theta_i = np.pi/np.sqrt(3), gamma = 1, H_I = 0*GeV):
        """Calculate axion decay constant f_a corresponding to a particular value of f_ax = Omega_axion/Omega_DM,
           in the pre-inflationary scenario, including backreaction contribution.
        
            Arguments:
            ---------
                f_ax (float): Fraction of DM density in axions
                theta_i (float): Initial misalignment angle (default: pi/sqrt(3))
                gamma (float): Dilution factor (default: 1.0)
                H_I (float): Energy scale of inflation (default: 0.0)
        
            Returns:
            ---------
                f_a (float): Value of axion decay constant.
        """
    
        obj = lambda y: (self.calc_Omegah2_Axion(y*GeV, theta_i, gamma, H_I) - f_ax*Omegah2_dm)**2
        res = optimize.minimize(obj, 1e9, method='Nelder-Mead', options={'disp': True})
        return res.x[0]*GeV


    def calc_fa_AxionDensity_Buschmann(self, f_ax):
        """Calculate axion decay constant f_a corresponding to a particular value of f_ax = Omega_axion/Omega_DM,
           in the post-inflationary scenario, following Buschmann et al (https://arxiv.org/abs/1906.00967)
        
             Arguments:
             ---------
                 f_ax (float): Fraction of DM density in axions
        
             Returns:
             ---------
                 f_a (float): Value of axion decay constant.
        """
        f_mid = 2.25e11*GeV*(f_ax)**(1/1.187)
        err = (11.0/25.2)
        f_min = f_mid*(1-err)
        f_max = f_mid*(1+err)
        return f_min, f_max


    #Dark Photon calculations
    #-------------------------------------

    def calc_Omegah2_DP_GluonFusion(self, m_dp, f_a, T_RH):
        """Calculate Dark Photon relic abundance Omega_DP*h^2 due to gluon fusion.
    
            Arguments:
            ---------
                m_dp (float): Dark Photon mass.
                f_a (float): Value of axion decay constant.
                T_RH (float): Reheating temperature.
   
            Returns:
            ---------
                Omegah2_dp (float): Dark Photon abundance Omega_DP*h^2
        """
        return (self.e_D*self.D_psi/0.3)**4*(m_dp/(10*eV))*(1e9*GeV/f_a)**4*(T_RH/(1e9*GeV))**3

    def calc_fa_DPDensity_GluonFusion(self, m_dp, f_dp, T_RH):
        """Calculate axion decay constant f_a corresponding to a particular value of f_DP = Omega_DP/Omega_DM, 
           for Dark Photon (DP) production through gluon fusion.
    
            Arguments:
            ---------
                m_dp (float): Dark Photon mass.
                f_dp (float): Dark Photon fraction f_DP = Omega_DP/Omega_DM
                T_RH (float): Reheating temperature.
   
            Returns:
            ---------
                f_a (float): Value of axion decay constant.
        """
        Omega_dp = f_dp*Omega_dm
        g_star   = 106.75
    
        A = (s0/rho_crit)*(135*np.sqrt(10)/(8*np.pi**13))*alpha_s**2*N_C**2*Mpl_reduced
        A *= m_dp*(self.e_D*self.D_psi*self.PQ_Phi)**4/(g_star**(3/2)*Omega_dp)
        return (A*T_RH**3)**(1/4)

    def calc_m_dp_min(self, f_dp):
        """Calculate minimum Dark Photon mass to achieve Freeze-in through gluon fusion.
        
            Arguments:
            ---------
                f_dp (float): Desired Dark Photon fraction f_DP = Omega_DP/Omega_DM
   
            Returns:
            ---------
                m_dp_min (float): Minimum possible Dark Photon mass.
        
        """
        Omegah2_dp = f_dp*Omegah2_dm
        zeta = special.zeta(3)
        g_star   = 106.75
        return Omegah2_dp*(48*np.pi**4*(rho_crit/h**2)*g_star)/(1080*zeta*s0)


    #Dark Photon Decay widths
    #-------------------------------------
    def calc_width_dp_e_e(self, m_dp, epsilon):
        """Decay width of the Dark Photon into e+ e-.
        
            Arguments:
            ---------
                m_dp (float): Dark Photon Mass.
                epsilon (float): Dark Photon kinetic mixing.
   
            Returns:
            ---------
                Gamma_dp (float): Dark Photon decay width.
        """
        if (m_dp < 2*m_e):
            return 0
        else:
            return ((epsilon*e)**2/(12*np.pi))*m_dp*(1 - 4*m_e**2/m_dp**2)**(1/2)
    
    def calc_width_dp_p_a(self, m_dp, epsilon, f_a):
        """Decay width of the Dark Photon into photon + axion.
        
            Arguments:
            ---------
                m_dp (float): Dark Photon Mass.
                epsilon (float): Dark Photon kinetic mixing.
                f_a (float): Axion decay constant
   
            Returns:
            ---------
                Gamma_dp (float): Dark Photon decay width.
        """
        G_apdp = self.calc_G_apdp(f_a, epsilon)
        m_a = self.fa_to_ma(f_a)
    
        if (m_a > m_dp):
            return 0
        else:
            return ((G_apdp)**2/(96*np.pi))*m_dp**3*(1-m_a**2/m_dp**2)**3
    
    def calc_width_dp_3p(self, m_dp, epsilon):
        """Decay width of the Dark Photon into 3 photons
        
            Arguments:
            ---------
                m_dp (float): Dark Photon Mass.
                epsilon (float): Dark Photon kinetic mixing.
   
            Returns:
            ---------
                Gamma_dp (float): Dark Photon decay width.
        """
        return (5e-8)*epsilon**2*(e**2/(4*np.pi**2))**4*(m_dp**9/m_e**8)
    
    def calc_width_dp_total(self, m_dp, epsilon, f_a):
        """Total decay width of the Dark Photon (into e+ e-, photon + axion or 3 photons)
        
            Arguments:
            ---------
                m_dp (float): Dark Photon Mass.
                epsilon (float): Dark Photon kinetic mixing.
                f_a (float): Axion decay constant
   
            Returns:
            ---------
                Gamma_dp (float): Dark Photon decay width.
        """
        return self.calc_width_dp_e_e(m_dp, epsilon) + self.calc_width_dp_p_a(m_dp, epsilon, f_a) + self.calc_width_dp_3p(m_dp, epsilon)
    
    @np.vectorize
    def calc_lifetime_dp(self, m_dp, epsilon, f_a):
        """Dark Photon lifetime (1/decay width.)
        
            Arguments:
            ---------
                m_dp (float): Dark Photon Mass.
                epsilon (float): Dark Photon kinetic mixing.
                f_a (float): Axion decay constant
   
            Returns:
            ---------
                tau (float): Dark Photon decay lifetime
            
        """
        return 1/self.calc_width_dp_total(m_dp, epsilon, f_a)
    
    
    #Dark Photon Decay widths
    #-------------------------------------
    def calc_width_a_p_p(self, f_a):
        """Axion decay width into 2 photons.
        
            Arguments:
            ---------
                f_a (float): Axion decay constant
   
            Returns:
            ---------
                Gamma_a (float): axion decay width
        """
        G_app = self.calc_G_app(f_a)
        m_a = self.fa_to_ma(f_a)
        return ((G_app**2)/(64*np.pi))*m_a**3
    
    def calc_width_a_p_dp(self, m_dp, epsilon, f_a):
        """Axion decay width into photon + dark photon.
        
            Arguments:
            ---------
                m_dp (float): Dark Photon Mass.
                epsilon (float): Dark Photon kinetic mixing.
                f_a (float): Axion decay constant
   
            Returns:
            ---------
                Gamma_a (float): axion decay width
        """  
        G_apdp = self.calc_G_apdp(f_a, epsilon)
        m_a = self.fa_to_ma(f_a)
        if (m_a < m_dp):
            return 0
        else:    
            return ((G_apdp**2)/(32*np.pi))*m_a**3*(1 - m_dp**2/m_a**2)**3

    def calc_width_a_dp_dp(self, m_dp, epsilon, f_a):
        """Axion decay width into two dark photons.
        
            Arguments:
            ---------
                m_dp (float): Dark Photon Mass.
                epsilon (float): Dark Photon kinetic mixing.
                f_a (float): Axion decay constant
   
            Returns:
            ---------
                Gamma_a (float): axion decay width
        """
        G_adpdp = self.calc_G_adpdp(f_a, epsilon)
        m_a = self.fa_to_ma(f_a)
    
        if (m_a < 2*m_dp):
            return 0
        else:
            return ((G_adpdp**2)/(64*np.pi))*m_a**3*(1 - 4*m_dp**2/m_a**2)**(3/2)
    
    def calc_width_a_total(self, m_dp, epsilon, f_a):
        """Total decay width of the axion (into 2 photons, photon + dark photon, or 2 dark photons)
        
            Arguments:
            ---------
                m_dp (float): Dark Photon Mass.
                epsilon (float): Dark Photon kinetic mixing.
                f_a (float): Axion decay constant
   
            Returns:
            ---------
                Gamma_a (float): Total axion decay width
        """
        return self.calc_width_a_p_p(f_a) + self.calc_width_a_p_dp(m_dp, epsilon, f_a) + self.calc_width_a_dp_dp(m_dp, epsilon, f_a)

    @np.vectorize
    def calc_lifetime_a(self, m_dp, epsilon, f_a):
        """Axion lifetime (1/decay width.)
        
            Arguments:
            ---------
                m_dp (float): Dark Photon Mass.
                epsilon (float): Dark Photon kinetic mixing.
                f_a (float): Axion decay constant
   
            Returns:
            ---------
                tau (float): Axion decay lifetime
            
        """
        return 1/self.calc_width_a_total(m_dp, epsilon, f_a)
    
    
    
    
    
    
