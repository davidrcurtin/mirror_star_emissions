import numpy as np
from module import bisection, save_nugget_file, save_dictionary, loop_bisection, convert_dictionary, trim_dictionary, interpolate_paramspace
from plotting import plot_paramspace_tau
import os
import timeit

# Path of folder containing script outputs (relative to this file's location) 
output_path = 'paramspace_test'

# Set parameters
rho_c_MS_fracs  = [1e-5, 0.01, 0.1, 1.0, 10.0, 100.0, 1e5]
xis             = np.logspace(-26, -16, 21)
T_cs            = np.logspace(2,8, 60)
L_ratio_error   = 0.001
solver_tol      = 1e-5 

# Approximate runtime
def test():
    bisection(rho_c_MS_fracs[len(rho_c_MS_fracs)//2], xis[len(xis)//2], T_cs[len(T_cs)//2], L_ratio_error, solver_tol=solver_tol)
approx_sec_p_nug = timeit.timeit(test, number=1, setup='gc.enable()')
print(f'Trying to find {len(rho_c_MS_fracs) * len(xis) * len(T_cs)} nuggets. At ~{approx_sec_p_nug:.1f}s per nugget, this will take order {len(rho_c_MS_fracs) * len(xis) * len(T_cs) * approx_sec_p_nug /60/60:.1f} hours.', flush=True)

# Manipulate paths to be usable
output_path = os.path.join(os.path.dirname(__file__), output_path)
os.makedirs(output_path, exist_ok=True)
nugget_list_path = os.path.join(output_path, 'all_nuggets.nuggets')
dictionary_path  = os.path.join(output_path, 'contour_dict.json')
plot_path        = os.path.join(output_path, 'full_paramspace.png')
 
# Make the nuggets (Will take a long time, order 100 hours)
start = timeit.default_timer()
nuggets = loop_bisection(rho_c_MS_fracs, xis, T_cs, L_ratio_error, solver_tol=solver_tol)
print(f'Took {(timeit.default_timer()-start)/60/60:.2f} hours to generate data\n', flush=True)

# Save all the nuggets to one (large) file
save_nugget_file(nugget_list_path, nuggets)
print(f'NOTE: A file has been created containing all nugget data at: \n    {nugget_list_path}\n    With file size {os.path.getsize(nugget_list_path)/1024/1024:.3f}MB\n')

# Save the nuggets to a contour dictionary
contour_dict = convert_dictionary(nuggets)
save_dictionary(dictionary_path, contour_dict)
print(f'NOTE: A file has been created containing a dictionary of nugget data at: \n    {dictionary_path}\n    With file size {os.path.getsize(dictionary_path)/1024/1024:.3f}MB\n')

# Create the parameter space plot as shown in paper, first trimming it and interpolating in the gaps
contour_dict_trimmed = trim_dictionary(contour_dict)
contour_dict_interpolated = interpolate_paramspace(contour_dict_trimmed)

# Plot it
plot_paramspace_tau(contour_dict_trimmed, contour_dict_interpolated, plotpath=plot_path)
print(f'NOTE: A plot of the parameter space has been saved at the path: \n    {plot_path}')