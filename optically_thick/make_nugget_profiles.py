from plotting import plot_nugget_paper
from module import Nugget, save_nugget_file

convective_path = 'profiles/convective'
radiative_path = 'profiles/radiative'

# Values for rho_c and T_c were chosen by hand from a selection of nuggets made using the bisection method
convective_ex = Nugget(1, 1e-21, 0.00034261281441637663, 27074.199332437096)
radiative_ex = Nugget(1, 1e-21, 0.0021121561632486117, 93260.3346851357)

save_nugget_file(convective_path +'.nugget', convective_ex)
save_nugget_file(radiative_path +'.nugget', radiative_ex)

plot_nugget_paper(convective_ex, convective_path + '.png')
plot_nugget_paper(radiative_ex, radiative_path + '.png')