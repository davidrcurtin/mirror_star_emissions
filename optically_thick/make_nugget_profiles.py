import os
from plotting import plot_nugget_paper
from module import Nugget, save_nugget_file


# Define paths to where the output should be stored, relative to this file 
convective_path = os.path.join(os.path.dirname(__file__), 'profiles', 'convective_nugget')
radiative_path = os.path.join(os.path.dirname(__file__), 'profiles', 'radiative_nugget')


# Values for rho_c and T_c were chosen by hand from a selection of nuggets found using the bisection method
convective_ex = Nugget(1, 1e-21, 0.0003453333529746241, 27080.0)
radiative_ex = Nugget(1, 1e-21, 0.0020057765691810954, 100000.0)

save_nugget_file(convective_path +'.nugget', convective_ex)
save_nugget_file(radiative_path +'.nugget', radiative_ex)

plot_nugget_paper(convective_ex, convective_path + '.png')
plot_nugget_paper(radiative_ex, radiative_path + '.png')