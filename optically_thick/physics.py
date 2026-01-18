import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator

G = 6.674 *(10**(-11))
c = 299792458 # m/s
h_bar = 1.054571817*(10**(-34)) # J s
h = 6.62607015e-34
m_hydrogen = 1.6735575*(10**(-27)) #kg
m_e = 9.1093837139e-31
m_p = 1.67262192e-27

k_b = 1.38 * (10**(-23)) # J/K
sigma = 5.670374419 * (10**(-8)) # Stefan-Boltzmann constant 

rho_c_sun = 150000 # kg/m^3
L_sun = 3.846 * (10**26) # J/s
M_sun = 1.9891e30 # Kg
eV = 1.602176634e-19 # electron volt in joules

a_rad = 4 * sigma / c


# Import cloudy abundances to use
cloudy_data = np.genfromtxt('data/Cloudy_ISM_Abundances.txt', comments='#', skip_header=3,skip_footer=24, usecols=(1,))
number_fractions = cloudy_data/sum(cloudy_data)

# Weights from https://iupac.qmul.ac.uk/AtWt/
elemental_weights = m_p * np.array([1, 4, 6.94, 10.81, 12, 14, 16, 19, 20.18, 23, 24.3, 26.98,28.09, 30.97, 32.06, 35.45, 39.95, 39.1, 40.08,47.87, 50.94, 52, 54.94,55.85,58.93,58.69,63.55,65.38])
mass_fractions = (number_fractions * elemental_weights)/sum(number_fractions*elemental_weights)
X = mass_fractions[0]
Y = mass_fractions[1]
Z = sum(mass_fractions[2:])
# print(f'X,Y,Z = {X:.4f}, {Y:.4f}, {Z:.4f}')

mmm_unionized = sum(number_fractions*elemental_weights)
mmm = mmm_unionized/2
# print(f'mu/m_p = {mmm/m_p:.4f}')


def g_MS(rho_c_MS, r):
    return (4/3) * np.pi * G * rho_c_MS * r

def g_self(m, r):
    return G * m * (r**(-2))

def hminusopacity(rho, T):
    return 0.1 * np.sqrt(rho/10) * np.power(T/3250, 9) #in m^2/kg

def freefreeopacity(rho,T):
    return (X + Y) * (1 + X) * (rho/1000) * np.power(2.87 * (10**6) / T, 3.5) *0.1 #in m^2/kg

def boundfreeopacity(rho,T):
    return (Z/0.02) * (1 + (X/7)) * (rho/1000) * np.power(6.22 * (10**6) / T, 3.5) *0.1 #in m^2/kg

def electron_scattering_opacity():
    return (1+X)*0.02 # 1cm^2/g=0.1m^2/kg 

saha_const = 2/(h_bar**3/(2*np.pi*m_e*k_b)**1.5) # 2/electron wavelength^3 without the T^3/2
def H_ionization_fraction(T, rho):
    exponential = np.exp(-13.59844*eV/(k_b*T))

    # ionized has e- & p+ which both have 2 states so g_1=4, degeneracy in H atom is 2n^2 times 2 for proton spin and n=1 so g_0=4 
    degeneracy_fraction = 1

    saha_eq = saha_const*(T**1.5)*degeneracy_fraction*exponential

    # This is equivalent to below method, but avoids division by zero/overflow
    big_X = np.divide((X*rho/m_hydrogen), saha_eq, where= saha_eq >1e-50, out = np.zeros_like(saha_eq))
    return 2/(1+np.sqrt(1+4*big_X))

    # Following https://www.physics.unlv.edu/~jeffery/astro/educational_notes/100_saha.pdf somewhat
    # saha_eq=n_e*n_ionized/n_unionized, but here n_ionized = n_e, and the total is n=n_e+n_unionized
    # so 1/x=: saha_eq/n_e = n_e/n-n_e = x/X-x for x=n_e/saha, X=n/saha, so that I=x/X is the ionization fraction
    # thus x = X-x/x => x^2+x-X = 0, and can use quadratic equation
    # x = -1 + sqrt(1+4X)/2 = 2X/1+sqrt(1+4X) so I = 2/1+sqrt(1+4X)
    # Notice X = rho_H/rho = n*m_h/rho, so n = X*rho*m_h
    n = X*rho/m_hydrogen 
    big_X = n/saha_eq

    I = 2/(1+np.sqrt(1+4*big_X))
    return I

def kappa_old(T, rho):
    hminus = hminusopacity(rho, T)
    free = freefreeopacity(rho, T) + boundfreeopacity(rho, T)
    scattering = H_ionization_fraction(T,rho) * electron_scattering_opacity()

    result = np.nan_to_num(np.where(free < hminus, free, hminus)) + scattering# + molecular 
    return result

# Import Freedman Opacity table and make an interpolator from it
T_freedman, rho_freedman, kap_freedman = np.loadtxt('data/freedman_opacities_2014_[M_H]=0.txt', skiprows=38, unpack=True, usecols = (0,2,3))
kappa_interpolator_freedman = CloughTocher2DInterpolator(list(zip(np.log10(rho_freedman*1000), np.log10(T_freedman))), np.log10(kap_freedman/10))
del T_freedman, rho_freedman, kap_freedman
def kappa_freedman(T, rho):
    points = np.stack([np.log10(rho), np.log10(T)], axis = -1)
    logkap = kappa_interpolator_freedman(points)
    if np.shape(T) == (): return logkap[0]
    return logkap


def kappa(T, rho):
    logkap_freedman = kappa_freedman(T, rho)
    logkap_import = logkap_freedman
    kap_old = kappa_old(T, rho)
    return np.where(np.isfinite(logkap_import),10**logkap_import, kap_old) # if the imported exists, use it, otherwise default to analytic

def T__pressure(p_rad):
    return np.power(3 * p_rad / a_rad, 0.25)

def P_rad(T, rho):
    return (1/3) * a_rad * (T**4)

def P_gas(T, rho):
    return (k_b/mmm) * rho * T

def P(T, rho):
    return P_rad(T, rho) + P_gas(T, rho)

def rho__pressure(p, p_rad):
    return (p-p_rad) * mmm / (k_b * T__pressure(p_rad))

def rho_convective(rho_c, R, r):
    return rho_c * ((1 - ((r/R)**2))**(3/2))

def T_convective(T_c, R, r):
    return T_c * (1 - ((r/R)**2))

def L_heating(rho_c_MS, xi, m):
    alpha = (X/m_hydrogen) * (rho_c_MS/rho_c_sun) * xi
    return m * alpha

def L_blackbody(R, T):
    return 4 * np.pi * (R**2) * (T**4) * sigma


def tau_capture(M, R, rho_c_MS, M_MS = M_sun):
    '''Return value is in Myr, 
    M & M_MS in kg, 
    R in m, 
    rho_c_MS in kg/m^3 '''
    return 65.0 * (M/1e17) * (1e7/R)**2 * (M_sun/M_MS)**(2/3) * (rho_c_sun/rho_c_MS)**(1/3)

def tau_capture_wd(M, R, rho_c_MS, r_H = 1):
    '''Return value is in Myr, 
    M & M_MS in kg, 
    R in m, 
    rho_c_MS in kg/m^3 '''
    return 1e2 * (M/1e17) * (1e7/R)**2 * (rho_c_sun/rho_c_MS)**(1/3) * r_H**(4/3)


# Kippenhhan & weigert pg 124 eq 13.12/13.15
def nabla_ad(beta): # = (dlnT/dlnP)_ad
    term = (1-beta) * (4+beta)
    num = 1 + (term/(beta**2))
    denom = 2.5 + 4*(term/(beta**2))
    result = num/denom
    return result

def gamma_ad(beta):
    alpha = 1 / beta
    delta = (4/beta) - 3
    denom = alpha - (delta * nabla_ad(beta))
    return 1/denom


# Gradients
def grad_P_rad_rad(rho_c_MS, xi, p_rad, rho, m, r):
    T = T__pressure(p_rad)
    L = L_heating(rho_c_MS, xi, m)
    F_rad = L / (4*np.pi*(r**2))
    return -1 * kappa(T, rho) * rho * F_rad / c

def grad_P_rad_ad(rho_c_MS, p, p_rad, rho, r):
    beta = 1 - (p_rad/p)
    return -4* g_MS(rho_c_MS, r) * rho * nabla_ad(beta) * p_rad/p

def dmdr(rho, r):
    return 4 * np.pi * (r**2) * rho

# Pols 5.18
def nabla_rad(rho_c_MS, xi, T, rho, m, r):
    constants = 3/ (16 * np.pi * a_rad * c* G)
    L_r = L_heating(rho_c_MS, xi, m)
    m_tot = (4/3) * np.pi * (r**3) * rho_c_MS + m
    return constants * kappa(T, rho) * L_r * P(T, rho) / (m_tot * (T**4)) 

def is_convectively_unstable(rho_c_MS, xi, T, rho, m, r):
    # if not np.isfinite(rho):
    #     return True
    beta = 1 - (P_rad(T, rho)/P(T, rho))
    return nabla_rad(rho_c_MS, xi, T, rho, m, r) > nabla_ad(beta)