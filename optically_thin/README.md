This is the public repository for "Electromagnetic Signatures of Mirror Stars" by Isabella Armstrong, Berkin Gurbuz, David Curtin, Christopher D. Matzner. (https://arxiv.org/abs/2311.18086)

For further questions regarding the following, please contact either of the first authors by email:\
Berkin Gurbuz (berkin.gurbuz@mail.utoronto.ca, berkingurbuz@gmail.com)\
Isabella Armstrong (isabella.armstrong@mail.utoronto.ca)

Here we provide all files for you to replicate our results, including data files directly obtained from Cloudy, as well as further figures that were not present on the paper, such as spectra and profile (Rho(r), T(r)) of each nugget. Within each DM density file, you'll also find a "Properties" .csv file. This is a master file that contains any and all physical properties of each nugget, ranging from xi and mass to line ratios and Gaia magnitudes. Finally, we provide a tutorial on how to exactly replicate our results using Cloudy yourselves, including necessary conversions between Cloudy's conventions and those utilized in our paper.

#Properties#\
The properties file gives the following data on each nugget per column\
A) log(Xi) - Please see the paper for Xi's definition as the heating rate\
B) Log(Rho) - This is the core density (g/cm^3) of the nugget\
C) Log(Lumi) - Luminosity (L_Sun) of the nugget\
D) Log(Mass) - Mass (grams)\
E) Log(Temp) - Core temperature (K) of the nugget\
F) Gaia magnitude G\
G) Gaia G_bp - Grp color\
H) Log(Radius) - (cm, cutoff) Radius of the nugget\
I Through S) Log(Line Emission) - Energy of line emissions\
T through V) Log(Line Ratio) - Key line ratios

Making use of the "Plotted" and "Interpolated" propety lines, which utilize our physical cutoffs as described in the paper, you can reproduce figure 3.

#Spectra and Profile#\
The figures from spectra and profile of each nugget following the naming convention:\
Nugget(X)_Xi(Y)_Mass(Z)_T.pdf and Nugget(X)_Xi(Y)_Mass(Z)_Rho.pdf\
Where (X) is the 'nugget number', aligned with rows of the properties master file, and (Y), and (Z) are the heating rate and mass (g) of the nugget, respectively. Note the difference in units for mass between properties file and the naming convention, for better identification of nuggets via phase space plots.

Do note that the spectra figures given here are not complete, and are instead focused on visible frequency. If you need the full spectrum of the nugget, please make use of the provided Cloudy output continuum file (.con), as described below.

#Cloudy Files#\
We provide our .con continuum, .out output, .ovr overview, line.txt line ratios, and .prs pressure files directly as obtained by Cloudy. Additionally, we provide our .in input cards for anybody wishing to replicate calculations of our exact nuggets using a newer version of Cloudy, or otherwise, themselves. We recommend anybody intending to work directly with Cloudy files to read Cloudy's guide document "Hazy1", which can be found on Cloudy's official page: https://gitlab.nublado.org/cloudy/cloudy/-/wikis/home. Otherwise, you can see the section below, where we define each column of each file that we make use of, and how they can be converted for more practical calculations.

#Replicating Our Results#\
Here we give an in-depth overview of recreating our results. We begin by introducing Cloudy's input cards, how we set them up, how we can 'hack' them to replicate the physical conditions of a nugget in hydrostatic equilibrium that is heated by the dark mirror star, and how we iteratively run them. Then, we go over defining any necessary outputs from Cloudy, their units and 'cases' based on our input cards, and how we convert them to be used in our calculations of nugget properties. By the end of this section, you should be able to create your Cloudy input cards, run Cloudy over the phase space to replicate our data and follow our conversion methods to recreate the figures of our paper.

We use Cloudy to find the hydrostatic equilibrium solution of the nuggets and obtain their spectra and profiles. Thus, the main 'gas' given to Cloudy is just the nugget, and all the following outputs from Cloudy are for the nugget only. However, we keep the mirror star in the picture by introducing a 'background' heating rate on the gas and placing the gas in a gravitational potential equivalent to that of the mirror star's core.

Here is an example input card we give to Cloudy for our nuggets:

##############################################\
hextra -11.5 density\
CMB\
#--------------------\
sphere\
radius 0 10\
set dr 6\
hden 11.5\
#--------------------\
abundances "ISM.abn" no grains\
#--------------------\
constant pressure\
set pressure convergence 0.01\
set pressure ionize 20000\
gravity spherical\
gravity external 9.9e19 0.1 3\
#--------------------\
stop temperature 4.9 exceeds\
iterate 3\
#--------------------\
stop zone 10000\
print last iteration\
save continuum "m.con" units micron last\
save overview ".ovr" last\
save pressure ".prs" last\
save line list "line.txt" "LineList.dat"\
#--------------------\
luminosity 0.0001 linear \
blackbody t = 0.0001 K linear\
##############################################

Let's go over these lines, and explain why we chose to make use of them:
- Hextra X density: This is how we introduce the heating of the nugget due to the mirror star (photon dark-photon kinetic mixing). The density keyword makes this heating scale with the local hydrogen density. Here, X is the log of the volume heating rate (erg cm^-3 s^-1). We compute the conversion to be xi = 10^(X - 7) (rho_solar/rho_DMcore). For example, for DM core density = 0.1 rho_solar, we have that xi = 10^(X - 7) * 1/0.1 = 10^(X - 6).
- CMB: This introduces background CMB radiation. We choose to include this as otherwise Cloudy returns quite a few errors.
- Sphere: This tells Cloudy that we are interested in a spherically symmetric geometry. The other option is a plane-like solution, but we are not interested in that.
- Radius X Y: This creates a 'vacuum box' in which the solution will the contained within. X and Y are then the log of the inner and outer radius in cm of this box, respectively.
- Set dr X: This breaks the space into radius increments (called zones) of length 10^X cm. Although Cloudy will shorten/lengthen zone sizes to minimize error, we chose to keep it constant as it helped eliminate a bunch of errors.
- hden X: This determines the core density 10^X g/cm^3 of the nugget. Hydrostatic equilibrium is solved outwards starting from this core density.
- abundances "ISM.abn" no grains: This tells Cloudy to set the elemental composition of the nugget to the interstellar medium's abundances. This takes an average of the warm and cold phases of the ISM.
- Constant pressure: This command tells Cloudy to solve for hydrostatic equilibrium.
- Set pressure convergence X: This is our allowed relative error (X = 0.01 = 1%) between the zones of our nugget.
- Set pressure ionize X: This is the number of calls to the ionization table that are made. Keeping this number very large resolves a particular error code, without adding any significant run time.
- Gravity spherical: Creates a spherically symmetric gravitational potential.
- Gravity external X Y Z: Determines the strength of the gravitational potential due to the mirror star. X is the mass of the source, Y is the extent in parsecs, and Z is the power law distribution of the mass. The huge mass alongside the 0.1 parsec extent imitates the constant dark matter density, where the mass is picked accordingly so that this constant is as desired. For example, 9.9e19g at 0.1 parsecs corresponds to DM core density = 0.01 rho_solar.
- Stop temperature X exceeds: Cloudy requires some stop condition. Our primary stop condition is based on the corona temperature of the Corona. At the end of its stable radius, the nugget's temperature grows several orders of magnitude. We stop this from crashing our code by setting a stop condition at 10^X Kelvin. Cloudy completes an iteration of the zones when a zone exceeds 10^4.9 Kelvin.
- Iterate X: This is how many times Cloudy will run through the zones of the nugget, reducing errors and equilibrating the zones each time. This is greatly at the cost of run time, so we set X = 3 iterations, as recommended by Hazy1.
- Stop Zone X: This is a safety stop condition. If a nugget never reaches the 10^4.9 Kelvin criteria and keeps growing larger and larger, this condition will stop it. We have X = 10000, but nuggets usually don't get this far.
- Print last iteration: This command gives us back an output file of the data of the last and least erroneous iteration.
- Save continuum "m.con" units micron last: This tells Cloudy to save an output file of the continuum/spectra of the nugget. With 'units micron' set, this file has columns 0 and 6 corresponding to the wavelength (in microns) and the SED, respectively.
- Save overview ".ovr" last: This tells Cloudy to save an output file of, amongst other things, the density and temperature profile of the nugget. Columns 0, 1, and 3 correspond to the radius, temperature, and density, respectively.
- Save pressure ".prs" last: Similar to the above, this tells Cloudy to save the pressure within the nugget versus the depth into the nugget.
- save line list "line.txt" "LineList.dat": Cloudy's output does not keep track of each line intensity, and certain lines are dropped depending on their relative energy to H_beta. To force Cloudy to keep track of these lines, we give them to Cloudy in "LineList.dat", and Cloudy returns their intensity in "Nugget(X)_Xi(Y)_Mass(Z)_line.txt". LineList.dat can be found within each dark matter density subdirectory and should be saved in a common directory with .in files when iteratively running Cloudy.
- Luminosity X linear & Blackbody t = Y K linear: This is how we trick Cloudy into switching to the "Luminosity Case". By introducing a Blackbody of negligible temperature/emissions at the center of our nugget, Cloudy detects input radiation given in Luminosity units. As a result, Cloudy switches to the case where it outputs the continuum in luminosity (nu L_nu) units. By default, Cloudy returns the continuum in the "Intensity Case", and since hextra works as background heating instead of incident radiation, it does not lead Cloudy to change to the "Luminosity Case". Therefore, it is crucial to put this dummy luminosity source into the code to get out the intensity units.

A phase space of constant DM density nuggets corresponds to a grid .in files with (hextra, hden) changing across the grid. Typically, Hextra runs from -20 to -6 and Hden runs from 6 to 19 across this phase space, with each nugget keeping any properties relating to the geometry, gravitational potential, and stop condition unchanged. Between DM density phase spaces, the mass (and thus gravitational potential) is scaled by factors of 10.

We make use of the Digital Research Alliance of Canada's Niagara cluster supercomputer. Within Niagara, we run shell code to iterative run the (hextra, hden) grid of the phase space. These shell codes can be found in our repository under mirror_star_emissions/optically_thin/Iteration. This code is specific to our installation of Cloudy on the cluster, and you'll have to rewrite these to fit your setup. The general idea revolves around 4 shell files, MakeNuggets.sh holds within it the input card with hextra and hden as variables. GenerateNuggets.sh calls to MakeNuggets.sh and submits the variable grid of (hextra, hden) and creates all .in files. SubmitCloudy.sh runs Cloudy on a script of arbitrary file name. SubmitJob.sh calls on SubmitCloudy.sh for each .in file in the directory. Thus, MakeNuggets.sh and GenerateNuggets.sh together create all input cards to be run, whereas SubmitCloudy.sh and SubmitJob.sh together iteratively run all these input cards.

Extracting physical properties from Cloudy output files is quite straightforward, and Hazy1 clearly defines each column. However, the columns relevant to our analysis are
- Radius: Column 0 of .ovr in cm
- Temperature: Column 1 of .ovr in Kelvin
- Density: Column 3 of .ovr in g/cm^3
- Wavelength: Column 0 of .con in microns. Should be turned to frequency to use with the continuum.
- Continuum: Column 6 of .con in Erg (nu L_nu). Divide out the frequency to obtain the luminosity per frequency.
- Pressure: Column 2 of .prs in Dyne cm^-2. Note that this column starts at the outer shell and moves inward.

All physical properties discussed in the paper (Mass, Luminosity, ...) can be integrated for using these columns.

We make use of GAIA DR3 passbands, available here: https://www.cosmos.esa.int/web/gaia/edr3-passbands to get G and Gbp-Grp using the spectra available in .con files.

Last updated: March 21st 2024
