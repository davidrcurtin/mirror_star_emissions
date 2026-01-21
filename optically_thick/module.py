import json
import numpy as np
import scipy.integrate
import scipy.interpolate
from physics import *
import pickle
from math import isclose

def nugget_profiles(rho_c_MS, rho_c, T_c, xi, unpack = False, convective_tagging = False, tol = 1e-4): 
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


def photosphere(r, T, rho, m, desired_tau = 2/3):
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
    # Check taus[-2] to avoid having r_photo = 0
    if taus[-2] < desired_tau:
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
    _rho_c_MS: float
    _rho_c_MS_frac: float
    _xi: float
    _T_c: float
    _rho_c: float

    _r: np.ndarray
    _m: np.ndarray
    _rho: np.ndarray
    _p: np.ndarray
    _p_rad: np.ndarray
    _convective_tag: np.ndarray

    _has_photosphere: bool
    _photo_index: int
    _R_photo: float
    _rho_photo: float
    _T_photo: float
    _M_photo: float

    def __init__(self, rho_c_MS_frac, xi, rho_c, T_c, solver_tol = 1e-5):
        '''        
        :param rho_c_MS_frac: Mirror Star central density, in solar central densities
        :param xi: Heating rate xi
        :param rho_c: Central nugget density, in kg/m^3
        :param T_c: Central nugget temperature, in K
        :param solver_tol: Tolerance to use in the stellar structure solver's step size estimation, smaller tolerance increases number of steps taken
        '''
        self._rho_c_MS_frac = rho_c_MS_frac
        self._xi = xi
        self._T_c = T_c
        self._rho_c = rho_c

        self._r, self._p, self._p_rad, self._rho, self._m, self._convective_tag =\
            nugget_profiles(rho_c_MS_frac*rho_c_sun, rho_c, T_c, xi, unpack=True, convective_tagging=True, tol = solver_tol)

        self._photo_index, self._R_photo, self._T_photo, self._rho_photo, self._M_photo =\
            photosphere(self._r, T__pressure(self._p_rad), self._rho, self._m)
        self._has_photosphere = self._photo_index != -1

    def rho_c_MS(self, cgs=False) -> float:
        '''
        The mirror star central density in kg/m^3

        :param cgs: If True, returns the result in cgs units (g/cm^3)
        '''
        rho_c_MS = self._rho_c_MS_frac * rho_c_sun
        if cgs: return rho_c_MS / 1e3 
        return rho_c_MS 
    
    def rho_c_MS_frac(self) -> float:
        '''The mirror star central density, in units of solar central density'''
        return self._rho_c_MS_frac
    
    def xi(self) -> float:
        '''The heating rate xi for the nugget'''
        return self._xi
    
    def T_c(self) -> float:
        ''' The nugget central density, in K'''
        return self._T_c
    
    def rho_c(self, cgs=False) -> float:
        '''Returns the nugget central density in kg/m^3

        :param cgs: If True, returns the result in cgs units (g/cm^3)
        '''
        if cgs: return self._rho_c / 1e3 
        return self._rho_c 

    def r(self, km=False) -> np.ndarray:
        '''Returns the radius values that index the other nugget profiles, in m 
        
        :param km: If true, returns the result in km'''
        if km: return self._r/1000
        return self._r

    def m(self) -> np.ndarray:
        '''Returns the mass values indexed by self.r(), in kg'''
        return self._m

    def rho(self, cgs=False) -> np.ndarray:
        '''Returns the density values indexed by self.r(), in kg/m^3\n
        self.rho()[0] == self.rho_c()
        
        :param cgs: If True, returns the result in cgs units (g/cm^3)'''
        if cgs: return self._rho / 1e3 
        return self._rho
    
    def T(self) -> np.ndarray:
        '''Returns the temperature values indexed by self.r(), in K'''
        return T__pressure(self.P_rad())
    
    def P(self) -> np.ndarray:
        '''Returns the total pressure values indexed by self.r(), in Pa = N/m^2'''
        return self._p
    
    def P_rad(self) -> np.ndarray:
        '''Returns the radiation pressure values indexed by self.r(), in Pa = N/m^2'''
        return self._p_rad

    def convective_tag(self) -> np.ndarray:
        '''Returns a np array of either 1 or 0, indexed by self.r(), where 1 means the nugget is convective at that radius and 0 means it is radiative.'''
        return self._convective_tag

    def has_photosphere(self) -> bool:
        '''Returns true if the nugget has a valid photosphere (total optical depth tau(0) > 2/3)'''
        return self._has_photosphere
    
    def photo_index(self) -> int:
        '''Returns the index of the nugget profiles closest to the photosphere location\n
        If the nugget does not have a photosphere, returns None'''
        if self.has_photosphere():
            return self._photo_index

    def r_photo(self, km = False) -> float:
        '''Returns the radius of the photosphere, in m \n
        If the nugget does not have a photosphere, returns None

        :param km: If true, returns the result in km'''
        if self.has_photosphere():
            if km: return self._R_photo/1000
            return self._R_photo
    
    def rho_photo(self) -> float:
        '''Returns the density at the photosphere, in kg/m^3 \n
        If the nugget does not have a photosphere, returns None'''
        if self.has_photosphere():
            return self._rho_photo
    
    def T_photo(self) -> float:
        '''Returns the temperature at the photosphere, in K \n
        If the nugget does not have a photosphere, returns None'''
        if self.has_photosphere():
            return self._T_photo
    
    def M_photo(self) -> float:
        '''Returns the nugget mass within the photosphere, in kg\n
        If the nugget does not have a photosphere, returns None'''
        if self.has_photosphere():
            return self._M_photo
    
    def L_photo(self) -> float:
        '''Returns the blackbody luminosity taken at the photosphere: 4pi*sigma*r_photo*T_photo^4\n
        If the nugget does not have a photosphere, returns None'''
        if self.has_photosphere():
            return L_blackbody(self.r_photo(), self.T_photo())
        
    def L_heat(self) -> float:
        '''Returns the heating luminosity at the photosphere: (X/m_hydrogen) * (rho_c_MS/rho_c_sun) * xi * M_photo\n
        If the nugget does not have a photosphere, returns None'''
        if self.has_photosphere():
            return L_heating(self.rho_c_MS(), self.xi(), self.M_photo())
    
    def L_ratio(self) -> float:
        '''Returns L_photo/L_heating\n
        If the nugget does not have a photosphere, returns None'''
        if self.has_photosphere():
            return self.L_photo()/self.L_heat()
        
    def log_L_ratio(self) -> float:
        '''Returns log10(L_photo/L_heating)\n
        If the nugget does not have a photosphere, returns None'''
        if self.has_photosphere():
            return np.log10(self.L_ratio())
    
    def g_photo(self, cgs = False) -> float:
        '''Returns the gravitational acceleration at the photosphere, with contributions from both the nugget and mirror star, in m/s^2\n
        If the nugget does not have a photosphere, returns None
        
        :param cgs: If True, returns the result in cgs units (cm/s^2)'''
        if self.has_photosphere():
            g_photo = g_MS(self.rho_c_MS(), self.r_photo()) + g_self(self.M_photo(), self.r_photo())
            if cgs: g_photo *= 100
            return g_photo
    
    def log_g_photo(self, cgs = False) -> float:
        '''Returns log10(g_photo)\n
        If the nugget does not have a photosphere, returns None
        
        :param cgs: If True, takes the log in cgs units (cm/s^2)'''
        if self.has_photosphere():
            return np.log10(self.g_photo(cgs))
        
    def total_tau(self) -> float:
        '''Returns the total optical depth of the nugget, integrating kappa*rho out to infinity'''
        if len(self.r()) < 5:
            return -1
        dtau_vals = kappa(self.T(), self.rho()) * self.rho()
        dtau = scipy.interpolate.make_interp_spline(self.r(), dtau_vals)

        return scipy.integrate.quad(dtau, 0, max(self.r()), limit=200, full_output=1)[0]

    
def save_nugget_file(filepath: str, nugget: Nugget):
    '''Writes to a pickle file at filepath'''
    with open(f'{filepath}', 'wb') as file:
        pickle.dump(nugget, file)

def load_nugget_file(filepath:str) -> Nugget:
    '''Load from a pickle file at filepath'''
    with open(f'{filepath}', 'rb') as file:
        nugget = pickle.load(file)
    return nugget

def bisection(rho_c_MS_frac, xi, T_c, L_ratio_error, log_rho_c_l=-10, log_rho_c_r=5, max_steps = 50, solver_tol = 1e-5) -> tuple[Nugget, list[Nugget], bool]:
    '''Performs a binary search in log_rho_c to find a nugget with given rho_c_MS, xi, and T_c with L_photo/L_heat within L_ratio_error of 1.\n 
    Exits if the max_steps count is surpassed or if the endpoints come within machine epsilon of eachother.
    
    :return: the final nugget, a list of all nuggets created, True if a satisfactory nugget was found '''

    def luminosity_function(log_rho_c):
        rho_c = 10**log_rho_c
        nugget = Nugget(rho_c_MS_frac, xi, rho_c, T_c, solver_tol=solver_tol)
        if nugget.has_photosphere():
            return nugget.log_L_ratio(), nugget
        else:
            return 1, nugget # Treat no photosphere as rho_c being too low
    
    log_L_ratio_error = np.log10(L_ratio_error + 1) #Ensures that the bisection in log(L_ratio) will give a result accurate to at least L_ratio_error
    max_rho_diff = np.finfo(float).eps * (abs(log_rho_c_l) + abs(log_rho_c_r))/2 # Numerical Recipies 9.1.1

    log_rho_c_mid = (log_rho_c_l + log_rho_c_r)/2
    log_L_ratio_mid, nugget = luminosity_function(log_rho_c_mid)
    nuggets = [nugget]
    num_steps = 1
    failed = False
    while abs(log_L_ratio_mid) > log_L_ratio_error and not failed: 
        if log_L_ratio_mid < 0:
            log_rho_c_r = log_rho_c_mid
        else: log_rho_c_l = log_rho_c_mid
        
        log_rho_c_mid = (log_rho_c_r + log_rho_c_l)/2
        log_L_ratio_mid, nugget = luminosity_function(log_rho_c_mid)
        nuggets.append(nugget)
        num_steps += 1
        failed = num_steps >= max_steps or isclose(log_rho_c_l,log_rho_c_r, abs_tol=max_rho_diff, rel_tol=0)
        
    return nugget, nuggets, not failed

def loop_bisection(rho_c_MS_fracs, xis, T_cs, L_ratio_error, log_rho_c_l=-10, log_rho_c_r=5, max_steps = 50, solver_tol = 1e-5) -> list[Nugget]:
    '''Perform bisection function over inputs with given parameters and return a list of the good nuggets'''
    good_nuggets = []
    for i, rho_c_MS_frac in enumerate(rho_c_MS_fracs):
        for j, xi in enumerate(xis):
            for k, T_c in enumerate(T_cs):
                print(f'{i+1}/{len(rho_c_MS_fracs)} rho_c_MS = {rho_c_MS_frac} rho_c_sun; {j+1}/{len(xis)} xi = {xi:.1e}; {k+1}/{len(T_cs)} T_c = {T_c:.2e}', end='\r', flush=True)
                good_nugget, nuggets, succesful = bisection(rho_c_MS_frac, xi, T_c, L_ratio_error, log_rho_c_l, log_rho_c_r, max_steps, solver_tol)
                if succesful: good_nuggets.append(good_nugget)
    print(f'Out of {len(rho_c_MS_fracs) * len(xis) * len(T_cs)} desired nuggets, {len(good_nuggets)} were found.                                                                   ', flush=True)
    return good_nuggets

def convert_dictionary(nuggets: list[Nugget]):
    '''Converts a list of Nugget objects into a dictionary of the form used by some other functions.\n
    It has keys of rho_c_MS_frac, with values being a list of 'contours' for that rho_c_MS with a single xi. 
    Here a contour is a list of properties (xis, log_rho_c_cgs, log_T_c, log_L_heating, log_R_photo, log_T_photo, log_rho_photo, log_M_photo, log_rad_radius, log_L_ratio, log_R_max, log_M_max, log_tau_total, log_g_photo) of the nuggets in the contour.
    '''

    def contour_from_nuggets(xi, nugget_list: list[Nugget]):
        properties = [np.log10(
            [
            n.rho_c(True), n.T_c(), n.L_heat(), n.r_photo(), n.T_photo(), n.rho_photo(), n.M_photo(),
            n.r()[np.nonzero(n.convective_tag()-1)[0][0] if len(np.nonzero(n.convective_tag()-1)[0]) > 0 else -1], # First index where its radiative, for non-on-liner see condense profiles
            n.L_ratio(), n.r()[-1], n.m()[-1], 
            n.total_tau(),
            n.g_photo()
            ]
            ) for n in nugget_list]
        return xi * np.ones(len(nugget_list)), *np.transpose(properties)

    nugget_dict = {}
    for nugget in nuggets:
        if nugget.rho_c_MS_frac() not in nugget_dict:
            nugget_dict[nugget.rho_c_MS_frac()] = {}

        if nugget.xi() not in nugget_dict[nugget.rho_c_MS_frac()]:
            nugget_dict[nugget.rho_c_MS_frac()][nugget.xi()] = []

        nugget_dict[nugget.rho_c_MS_frac()][nugget.xi()].append(nugget)
    
    contour_dict = {}
    for rho_c_MS_frac, contours in nugget_dict.items():
        if rho_c_MS_frac not in contour_dict:
            contour_dict[rho_c_MS_frac] = []

        for xi, contour in contours.items():
            contour_dict[rho_c_MS_frac].append(contour_from_nuggets(xi, contour))

    return contour_dict

def save_dictionary(filepath: str, contour_dict: dict):
    '''Save a dictionary (of the form given by convert_dictionary) as a json file at filepath'''
    with open(filepath, 'w') as file:
        json.dump({key: [[list(prop) for prop in contour] for contour in value] for key, value in contour_dict.items()}, file)

def load_dictionary(filepath: str):
    '''Load a dictionary from filepath, given it is of the for given by convert_dictionary'''
    with open(filepath, 'r') as file:
        contour_dict = {float(key): [[np.array(prop) for prop in contour] for contour in value] for key, value in json.load(file).items()}
    return contour_dict

def load_thin_data() -> tuple[dict,dict]:
    '''returns two dictionaries of the form {rho_c_MS_frac: {keys = 'xi', 'M_nugget', 'T_c', 'R', 'L', 'G_bp-G_rp'}}\n
    The first is the actual cloudy solutions, the second is the interpolated points.'''
    dictionary_cut, dictionary_interpolated = {}, {}
    for rho_c in [0.01, 0.1, 1, 10, 100]:
        data_strings_cut = np.genfromtxt(f'../optically_thin/{str(rho_c)}RhoDM/{str(rho_c)}PlottedNuggetProperties.csv',
                                      skip_header=1, delimiter=',', unpack = True, dtype=str)
        data_strings_interpolated = np.genfromtxt(f'../optically_thin/{str(rho_c)}RhoDM/{str(rho_c)}InterpolatedNuggets.csv',
                                      skip_header=1, delimiter=',', unpack = True, dtype=str)
        data = []
        for col in data_strings_cut:
            data.append(np.array([np.float64(value.strip('"')) for value in col]))
        
        cut_dict = {
            'xi': np.power(10, data[0]), 
            'rho_c': np.power(10, data[1]), # g/cm^3
            'M_nugget': np.power(10, data[2]), # grams
            'T_c': np.power(10, data[3]), #K
            'R': np.power(10, data[4])/100/1000, # cm --> km
            'L': np.power(10, data[5]), # L_sol
        }
        cut_dict['g'] = g_MS(rho_c*rho_c_sun, cut_dict['R']/1000) + g_self(cut_dict['M_nugget']/1000,cut_dict['R']/1000)
        dictionary_cut[rho_c] = cut_dict


        data = []
        for col in data_strings_interpolated:
            data.append(np.array([np.float64(value.strip('"')) for value in col]))

        interp_dict = {
            'xi': np.power(10, data[0]), 
            'rho_c': data[1], # g/cm^3 #interpolated data has rho_c saved as not log
            'M_nugget': np.power(10, data[2]), # grams
            'T_c': np.power(10, data[3]), #K
            'R': np.power(10, data[4])/100/1000, # cm --> km
            'L': np.power(10, data[5]), # L_sol
        }
        interp_dict['g'] = g_MS(rho_c*rho_c_sun, interp_dict['R']/1000) + g_self(interp_dict['M_nugget']/1000,interp_dict['R']/1000)*100 # m/s^2 --> cm/s^2
        dictionary_interpolated[rho_c] = interp_dict
        
    return dictionary_cut, dictionary_interpolated
