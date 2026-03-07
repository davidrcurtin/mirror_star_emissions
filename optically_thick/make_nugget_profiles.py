import os
from plotting import plot_nugget_paper
from module import Nugget, save_nugget_file


# Define paths to where the output should be stored, relative to this file 
convective_path = os.path.join(os.path.dirname(__file__), 'profiles', 'convective_nugget_2')
radiative_path = os.path.join(os.path.dirname(__file__), 'profiles', 'radiative_nugget_2')
plot_Pgas_Prad = False


# Values for rho_c and T_c were chosen by hand from a selection of nuggets found using the bisection method
convective_ex = Nugget(1, 1e-21, 0.00034261281441637663, 27074.199332437096)
radiative_ex = Nugget(1, 1e-21, 0.0021121561632486117, 93260.3346851357)

save_nugget_file(convective_path +'.nugget', convective_ex)
save_nugget_file(radiative_path +'.nugget', radiative_ex)

figsize = (4.5, 6.3) if plot_Pgas_Prad else (4.5, 5)
plot_nugget_paper(convective_ex, convective_path + '.png', figsize=figsize, plot_Pgas_Prad = plot_Pgas_Prad)
plot_nugget_paper(radiative_ex, radiative_path + '.png',  figsize=figsize, plot_Pgas_Prad = plot_Pgas_Prad)