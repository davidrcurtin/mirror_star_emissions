# Electromagnetic Signatures of Optically Thick Mirror Nuggets
This directory contains the code and cached data used to analyze the electromagnetic signatures of the optically thick nuggets generated in `mirror_star_emissions/optically_thick`. In particular, it contains the plotting workflow used to produce figures 3, 4, 5, and 6 of "Generalized Predictions for the Electromagnetic Signatures of Mirror Stars" by Franco Cabral, Stuart Williamson, David Curtin, Christopher D. Matzner (https://arxiv.org/abs/2604.00106).

For further questions regarding the following, please Franco Cabral by email (franco.cabral@mail.utoronto.ca).

This subdirectory is designed to be used together with the main `optically_thick` workflow. The plotting script here reads nugget-derived summary data from `optically_thick/paramspace`, example nugget files from `optically_thick/profiles`, and cached stellar/spectroscopic data from `optically_thick/em_signatures/data`.

## Description of files and directories
The main plotting entrypoint is <code>make_plots.py</code>. This script generates the EM-signature plots in <code>optically_thick/em_signatures/plots</code>, including:
- example interpolated spectra for the convective and radiative profile nuggets,
- photospheric temperature-surface gravity comparison plots,
- 2D Gaia Hertzsprung-Russell overlays, and
- 3D Gaia Hertzsprung-Russell overlays.

The <code>model_spectra_data.py</code> file contains utilities for reading MPS-Atlas atmosphere and spectra files, including model atmospheres, center-to-limb spectra, disk-integrated spectra, and the tabulated angular grid.

The <code>stellar_parameters.py</code> file contains helper routines for matching user-supplied metallicity, effective temperature, and surface gravity values to the nearest entries on the MPS-Atlas grid stored in <code>data/grids</code>.

The <code>data</code> directory contains local inputs used by the plotting scripts:
- <code>grids</code>: tabulated metallicity, effective temperature, and surface gravity grids used for MPS-Atlas interpolation,
- <code>gaia_cache</code>: cached Gaia-derived HR-diagram data used in the 2D and 3D comparison plots,
- <code>Cloudy ISM Abundances.txt</code>: ISM abundances used by the helper physics routines.

The <code>plots</code> directory is the default output location for <code>make_plots.py</code>. By default it contains the following subdirectories:
- <code>spectra_examples</code>: interpolated disk-integrated spectra for the example convective and radiative nuggets,
- <code>logg_T</code>: temperature-surface gravity comparison plots,
- <code>HR</code>: 2D and 3D Hertzsprung-Russell comparison plots,
- <code>data_exports</code>: machine-readable CSV files corresponding to the plotted model outputs.

The <code>plots/data_exports</code> directory contains:
- <code>spectra_examples/*.csv</code>: wavelength-flux tables for the example spectra plots, with columns <code>wavelength_nm</code> and <code>flux_erg_cm^-2_s^-1_nm^-1</code>,
- <code>logg_T/*_points.csv</code>: the plotted model points for the temperature-surface gravity figures, including core-density slice, photospheric temperature, surface gravity, and luminosity,
- <code>logg_T/*_boundaries.csv</code>: the fitted linear boundary parameters used to draw the lower- and upper-temperature envelopes in the temperature-surface gravity figures,
- <code>HR/stacked_HR_threecol_model_points.csv</code>: the model points used in the 2D HR plots, including color, absolute magnitude, surface gravity, nugget mass, heating rate, and photospheric temperature,
- <code>HR/stacked_HR_3D_model_points.csv</code>: the model points used in the 3D HR plots, including color, absolute magnitude, surface gravity, and photospheric temperature.

## Use of python scripts

To regenerate all EM-signature plots from the `optically_thick` directory, run:

<code>python em_signatures/make_plots.py --data-root em_signatures/data --plots-root em_signatures/plots</code>

The recommended working directory is <code>mirror_star_emissions/optically_thick</code>.

Running this script requires the following python packages:
- numpy
- scipy
- matplotlib
- pandas

The <code>make_plots.py</code> script expects the following repository inputs:
- <code>optically_thick/paramspace/contour_dict.json</code>: nugget summary contours used in the temperature-gravity and HR plots,
- <code>optically_thick/profiles/convective_nugget.nugget</code> and <code>optically_thick/profiles/radiative_nugget.nugget</code>: example nuggets used for the spectra plots,
- <code>optically_thick/em_signatures/data/grids/*.txt</code> and <code>optically_thick/em_signatures/data/gaia_cache/*.csv</code>: cached interpolation and Gaia data.

The spectra-example portion of the script also requires the external MPS-Atlas spectra grid. The full <code>set2</code> dataset is too large to store in this repository. It should be extracted so that the script can find a directory structure of the form:

<code>&lt;spectra-root&gt;/set2/MH.../teff.../logg.../mpsa_flux_spectra.dat</code>

If the spectra files are stored outside this repository, pass their parent directory with:

<code>python em_signatures/make_plots.py --data-root em_signatures/data --plots-root em_signatures/plots --spectra-root "C:\path\to\parent\of\set2"</code>

The script provides the following optional flags:
- <code>--set-type {set1,set2}</code>: choose which MPS-Atlas set to use for the spectra interpolation,
- <code>--spectra-root</code>: override the default spectra location of <code>&lt;data-root&gt;/spectra</code>,
- <code>--skip-spectra-examples</code>: skip the example spectra plots,
- <code>--skip-loggT</code>: skip the temperature-surface gravity plots,
- <code>--skip-hr</code>: skip the 2D HR plots,
- <code>--skip-hr-3d</code>: skip the 3D HR plots.

Because the HR and logg-T figures are built from <code>contour_dict.json</code>, updating the nugget set for this workflow typically means regenerating or replacing that file first, then rerunning <code>make_plots.py</code>.

When <code>make_plots.py</code> is rerun, it updates both the figure files and the CSV exports in <code>plots/data_exports</code>.


Last updated: July 9th 2026
