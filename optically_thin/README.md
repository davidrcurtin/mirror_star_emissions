This is the public repository for "Electromagnetic Signatures of Mirror Stars" by Isabella Armstrong, Berkin Gurbuz, David Curtin, Christopher D. Matzner.

For further questions regarding the following, please contact either of first authors by email:\
Berkin Gurbuz (berkin.gurbuz@mail.utoronto.ca, berkingurbuz@gmail.com)\
Isabella Armstrong (isabella.armstrong@mail.utoronto.ca)

Here we provide all files for you to replicate our results, including data files directly obtained from Cloudy, as well as further figures that were not present on the paper, such as spectra and profile (Rho(r), T(r)) of each nugget. Within each DM density file, you'll also find a "Properties" .csv file. This is a master file that contains any and all physical properties of each nugget, ranging from xi and mass to line ratios and Gaia magnitudes. Finally, we provide a tutorial on how to exactly replicate our results using Cloudy yourselves, including necessary conversions between Cloudy's conventions and those utilized in our paper.

#Properties#\
The properties file gives the following data on each nugget per column\
A) log(Xi) - Please see paper for Xi's definition as the heating rate\
B) Log(Rho) - This is the core density (g/cm^3) of the nugget\
C) Log(Lumi) - Luminosity (Erg/s) of the nugget\
D) Log(Mass) - Mass (Solar Masses)\
E) Log(Temp) - Core temperature (K) of the nugget\
F) Gaia magnitude G\
G) Gaia G_bp - Grp color\
H) Log(Radius) - (cutoff) Radius of the nugget\
I Through S) Log(Line Emission) - Energy of line emissions\
T through V) Log(Line Ratio) - Key line ratios\

#Spectra and Profile#\
The figures from spectra and profile of each nugget following the naming convention:\
Nugget(X)_Xi(Y)_Mass(Z)_T.pdf and Nugget(X)_Xi(Y)_Mass(Z)_Rho.pdf\
Where (X) is the 'nugget number', aligned with rows of the properties master file, and (Y), (Z) are the heating rate and mass (g) of the nugget, respectively. Note the difference in units for mass between properties file and the naming convention, for better identification of nuggets via phase space plots.

Do note that the spectra figures given here are not complete, and are instead focused on peak spectra. Please make use of the provided Cloudy output continuum file (.con), as described below.

#Cloudy Files#\
We provide our .con continuum, .out output, .ovr overview, line.txt line ratios and .prs pressure files directly as obtained by Cloudy. Additionally, we provide our .in input cards for anybody wishing to replicate calculations of our exact nuggets using a newer version of Cloudy, or otherwise, themselves. We absolutely recommend anybody intending to work directly with Cloudy files to read Cloudy's guide document "Hazy1", which can be found on Cloudy's official page: https://gitlab.nublado.org/cloudy/cloudy/-/wikis/home. Otherwise, you can see the section below, where we define each column of each file that we make use of, and how they can be converted for my practical calculations.

#Replicating Our Results#\
Here we give an in-depth overview on recreating our results. We begin by introducing Cloudy's input cards, how we set them up, and how we can 'hack' them to replicate the physical conditions of a nugget in hydrostatic equilibrium that is heated by the dark mirror star, and how we iteratively run them. Then, we go over defining any necessary outputs from Cloudy, their units and 'cases' based on our input cards, and how we convert them to be used in our calculations of nugget properties. By the end of this section, you should be able to create your Cloudy input cards, run Cloudy over the phase space to replicate our data, and follow our conversion methods to recreate the figures of our paper.

. . . coming very soon . . .

Last updated: March 12th 2024
