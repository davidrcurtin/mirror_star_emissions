import numpy as np
import scipy.integrate
import scipy.interpolate
from physics import *
import pickle

def get_profiles_combined_pressures_adaptive(rho_c_MS, rho_c, T_c, xi, unpack = False, convective_tagging = False, tol = 1e-4): 
    '''Returns rows of [[r, P, P_rad, rho, m], ...]. 
    If Unpack = False and convective_tagging = True, the tags are converted to floats as 1 or 0, 1 if convective 0 if not.
    Usually written as r, p, p_rad, rho, m, tag = get_profiles_combined_pressures_adaptive(...)
    '''
    p_c = P(T_c, rho_c)
    p_rad_c = P_rad(T_c, rho_c)

    reduction_factor = 1e-30
    p_stop = p_c * reduction_factor
    p_rad_stop = p_rad_c * reduction_factor
    rho_stop = rho_c * reduction_factor

    def sim_complete(p, p_rad, rho, count):
        return p <= p_stop or p_rad <= p_rad_stop or rho <= rho_stop or count >= 1e9

    def take_step(p, p_rad, rho, m, r, dr):
        # Hydrostatic
        dp = -1 * (g_MS(rho_c_MS, r) + g_self(m,r)) * rho * dr
        new_m = m + dmdr(rho, r) * dr
        new_p = p + dp
        new_r = r + dr
        
        # Consider convective step
        dp_rad_convective = grad_P_rad_ad(rho_c_MS, p, p_rad, rho, r) * dr
        p_rad_convective = p_rad + dp_rad_convective 
        rho_convective = rho__pressure(new_p, p_rad_convective)

        # Check if the convective step would be convectively unstable
        if is_convectively_unstable(rho_c_MS, xi, T__pressure(p_rad_convective), rho_convective, new_m, new_r):
            # If so, complete the convective step
            convective_used = True
            new_p_rad = p_rad_convective 
            new_rho = rho_convective

        else: # Otherwise take a radiative step
            convective_used = False
            dp_rad_radiative = grad_P_rad_rad(rho_c_MS, xi, p_rad, rho, m, r) * dr
            new_p_rad = p_rad + dp_rad_radiative
            new_rho = rho__pressure(new_p, p_rad)
        
        return np.array([new_p, new_p_rad, new_rho, new_m]), new_r, convective_used
    
    # Actually start the integration
    r, p, p_rad, rho, m, dr = 1e-6, p_c, p_rad_c, rho_c, 1e-6, 1e-8
    vals = [[r, p, p_rad, rho, m]]

    # Check if the first zone will be convective
    p_rad_convective_test = p_rad + grad_P_rad_ad(rho_c_MS, p, p_rad, rho, r) * dr/2
    p_test = p - g_MS(rho_c_MS, r) * rho * dr/2
    m_test = m + dmdr(rho, r) * dr/2
    r_test = r + dr/2
    convective_tag = [is_convectively_unstable(rho_c_MS, xi, T__pressure(p_rad_convective_test), rho__pressure(p_test, p_rad_convective_test), m_test, r_test)]

    count = 0
    while not sim_complete(p, p_rad, rho, count):
        count += 1

        y_1a, r_1a, convective_1a = take_step(p, p_rad, rho, m, r, dr/2)
        if np.any(y_1a < 0):
            break
        y_1b, r_1b, convective_1b = take_step(*y_1a, r_1a, dr/2)
        if np.any(y_1b < 0):
            break

        y_2, r_2, convective_2 = take_step(p, p_rad, rho, m, r, dr)

        solution_diff = y_1b - y_2
        error = solution_diff[0] # Take error as error in p
        rtol = tol*max(y_1b[0], y_2[0])
        thing = (rtol/(2*abs(error)))**0.5 if error != 0 else 1.5 
        dr = 0.9 * dr * min(max(thing, 0.2), 1.5)

        vals.append([r_1a, *y_1a])
        vals.append([r_1b, *y_1b])
        convective_tag.append(True if convective_1a else False)
        convective_tag.append(True if convective_1b else False)

        p, p_rad, rho, m = y_1b
        r = r_1b

    if convective_tagging:
        if unpack:
            return [*np.transpose(vals), np.array(convective_tag)]
        return np.append(np.array(vals), np.array(convective_tag).reshape((len(vals), 1)), axis = 1)
    
    if unpack:
        return np.transpose(vals)
    return np.array(vals)


def get_photosphere_quantities(r, T, rho, m, desired_tau = 2/3):
    '''Returns (photo_index, R_photo, T_photo, rho_photo, m_photo)'''
    if len(r) < 5:
        return -1, -1, -1, -1, -1
    # Interpolate data
    
    integrand_vals = kappa(T, rho) * rho

    # Interpolate the integrand: idk if its better to interpolate the individual functions first
    integrand_interp = scipy.interpolate.make_interp_spline(r, integrand_vals)

    def integrand(s):
        r_value = r.max() - s
        return integrand_interp(r_value)
    
    s_vals = np.geomspace(1, r.max(), 10000)
    
    # Integrate up to every point tau
    taus = scipy.integrate.cumulative_simpson(integrand(s_vals), x=s_vals, initial =0)

    # If tau doesnt reach desired value by core skip the rest and return error values
    if taus[-1] < desired_tau:
        return -1, -1, -1, -1, -1

    # Find first point from the exterior with high enough tau and convert to r
    photo_s_val = s_vals[taus > desired_tau][0]
    photo_r_val = r.max() - photo_s_val

    # Find the nearest index just for returning (not usually used)
    photo_index = np.where(r > photo_r_val)[0][0]

    # interpolate physical quantities to use the interpolated r value
    T_interp = scipy.interpolate.make_interp_spline(r, T)
    rho_interp = scipy.interpolate.make_interp_spline(r, rho)
    m_interp = scipy.interpolate.make_interp_spline(r, m)

    return photo_index, photo_r_val, T_interp(photo_r_val), rho_interp(photo_r_val), m_interp(photo_r_val)

class Nugget:
    rho_c_MS: float
    rho_c_MS_frac: float
    xi: float
    T_c: float
    rho_c: float

    r: np.ndarray
    m: np.ndarray
    rho: np.ndarray
    p: np.ndarray
    p_rad: np.ndarray
    convective_tag: np.ndarray

    has_photosphere: bool
    photo_index: int
    R_photo: float
    rho_photo: float
    T_photo: float
    M_photo: float

    def __init__(self, rho_c_MS_frac, xi, rho_c, T_c, solver_tol = 1e-5, rho_c_cgs = False):
        self.rho_c_MS_frac = rho_c_MS_frac
        self.xi = xi
        self.T_c = T_c
        if rho_c_cgs: self.rho_c = rho_c * 1e3
        else:       self.rho_c = rho_c

        profile_result = get_profiles_combined_pressures_adaptive(self.rho_c_MS(), self.rho_c, self.T_c, self.xi, unpack=True, convective_tagging=True, tol = solver_tol)
        self.r, self.p, self.p_rad, self.rho, self.m, self.convective_tag = profile_result

        photosphere_result = get_photosphere_quantities(self.r, self.T(), self.rho, self.m)
        self.has_photosphere = photosphere_result[0] != -1
        self.photo_index, self.R_photo, self.T_photo, self.rho_photo, self.M_photo = photosphere_result

    def rho_c_MS(self) -> float:
        return rho_c_sun * self.rho_c_MS_frac

    def T(self) -> np.ndarray:
        return T__pressure(self.p_rad)
    
    def L_photo(self) -> float:
        return L_blackbody(self.R_photo, self.T_photo)
        
    def L_heat(self) -> float:
        return L_heating(self.rho_c_MS(), self.xi, self.M_photo)
    
    def L_ratio(self) -> float:
        '''L_photo/L_heating'''
        return self.L_photo()/self.L_heat()
        
    def log_L_ratio(self) -> float:
        '''log10(L_photo/L_heating)'''
        return np.log10(self.L_ratio())
    
    def rho_cgs(self) -> np.ndarray:
        return self.rho / 1e3
    
    def g_photo(self, cgs = False) -> float:
        g_photo = g_MS(self.rho_c_MS(), self.R_photo) + g_self(self.M_photo, self.R_photo)
        if cgs: g_photo *= 100
        return g_photo
    
    def log_g_photo(self, cgs = False) -> float:
        return np.log10(self.g_photo(cgs))


def save_nugget_file(filepath: str, nugget: Nugget):
    '''Writes as a pickle file'''
    with open(f'{filepath}', 'wb') as file:
        pickle.dump(nugget, file)

def load_nugget_file(filepath:str) -> Nugget:
    '''Load from a pickle file at filepath'''
    with open(f'{filepath}', 'rb') as file:
        nugget = pickle.load(file)
    return nugget
