import numpy as np
from module import save_nugget_file, save_dictionary, loop_bisection, convert_dictionary, trim_dictionary
from plotting import plot_paramspace
from physics import L_sun
import os
from timeit import default_timer as timer

# Script inputs
nugget_list_path = 'paramspace/all_nuggets.nuggets'
dictionary_path  = 'paramspace/contour_dict.json'
plot_path        = 'paramspace/full_paramspace.png'
approx_sec_p_nug = 12

# Set parameters
rho_c_MS_fracs  = [1e-5, 0.01, 0.1, 1.0, 10.0, 100.0, 1e5]
xis             = np.logspace(-26, -16, 21)
T_cs            = np.logspace(2,8, 60)
L_ratio_error   = 0.001
solver_tol      = 1e-5 
print(f'Trying to find {len(rho_c_MS_fracs) * len(xis) * len(T_cs)} nuggets. At {approx_sec_p_nug}s per nugget, this will take {len(rho_c_MS_fracs) * len(xis) * len(T_cs) * approx_sec_p_nug /60:.2f} minutes.', flush=True)

# Make the nuggets (Will take a long time, order 100 hours)
start = timer()
nuggets = loop_bisection(rho_c_MS_fracs, xis, T_cs, L_ratio_error, solver_tol=solver_tol)
print(f'Took {(timer()-start)/60:.2f} minutes to generate data\n', flush=True)

# Save all the nuggets to one (large) file
save_nugget_file(nugget_list_path, nuggets)
print(f'NOTE: A file has been created containing all nugget data at: \n    {nugget_list_path}\n    With file size {os.path.getsize(nugget_list_path)/1024/1024:.3f}MB\n')

# Save the nuggets to a contour dictionary
contour_dict = convert_dictionary(nuggets)
save_dictionary(dictionary_path, contour_dict)
print(f'NOTE: A file has been created containing a dictionary of nugget data at: \n    {dictionary_path}\n    With file size {os.path.getsize(dictionary_path)/1024/1024:.3f}MB\n')

# Create the parameter space plot as shown in paper, first trimming it
contour_dict_trimmed = trim_dictionary(contour_dict)

# Plot it
plot_paramspace(contour_dict_trimmed, plotpath=plot_path)
print(f'NOTE: A plot of the parameter space saved at the path: \n    {plot_path}')