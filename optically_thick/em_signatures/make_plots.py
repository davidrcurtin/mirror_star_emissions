#!/usr/bin/env python3

"""
Reproducible plot generator for EM signatures/mirror-star contours.

Recommended usage:
  python em_signatures/make_plots.py --data-root em_signatures/data --plots-root em_signatures/plots

Inputs expected under --data-root:
  grids/grid_mh.txt
  grids/grid_teff.txt
  grids/grid_logg.txt
  contours/contours_new__rho_c_MS=...json
  nuggets/radiative_ex.nugget
  nuggets/convective_ex.nugget
  gaia_cache/gaia_observed_hr.csv
  gaia_cache/gaia_observed_hr_logg.csv
  gaia_cache/gaia_white_dwarfs_hr_logg.csv

Spectra data (IMPORTANT!)
- The full set2 spectra grid is too large to store in this repo.
- Download `set2.zip` from the MPS-Atlas dataset page:
  https://edmond.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.NJ56TR
- Extract it so that you have a subfolder structure in data/ like:
  spectra/set2/MH.../teff.../logg.../mpsa_flux_spectra.dat
- Note: this download is large (multi-GB).

White dwarf Gaia cache (this is just a note)
- The provided `gaia_cache/gaia_white_dwarfs_hr_logg.csv` was derived from the FITS catalog:
  https://warwick.ac.uk/fac/sci/physics/research/astro/research/catalogues/gaiaedr3_wd_main.fits.gz
- Note: the FITS file is also large.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.interpolate import interpn
from scipy.optimize import curve_fit, fsolve

# Local modules
from stellar_parameters import GridOfStellarParameters
from model_spectra_data import ReadData

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


import importlib.util
import sys
from pathlib import Path as _Path
_module_path = _Path(__file__).resolve().parents[1] / "module.py"
spec = importlib.util.spec_from_file_location("module", str(_module_path))
_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_module)
sys.modules['module'] = _module
load_nugget_file = _module.load_nugget_file

_plot_path = _Path(__file__).resolve().parents[1] / "plotting.py"
spec_plot = importlib.util.spec_from_file_location("repo_plotting", str(_plot_path))
_plot_mod = importlib.util.module_from_spec(spec_plot)
spec_plot.loader.exec_module(_plot_mod)
custom_cmap = _plot_mod.custom_cmap


# Plot style
plt.rcParams.update({
    "font.size": 24,
    "legend.fontsize": 15,
    "legend.title_fontsize": 16
})


# Physical constants
G = 6.674e-8  # cm^3 g^-1 s^-2
L_sun_W = 3.828e26  # W
L_sun_cgs = 3.828e33  # erg/s
sigma_cgs = 5.670374419e-5  # erg cm^-2 s^-1 K^-4

rho_c_sun = 150  # g/cm^3
logg_sun = np.log10(27400)  # Sun surface gravity in cgs


# Simple fit models
def linear_fit(x, a, b):
    return a * x + b


# Helpers
def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def format_rho(rho: float) -> str:
    # Return a mathtext-ready string so legend renders exponents correctly.
    # Use isclose to avoid float equality pitfalls.
    if np.isclose(rho, 1e-5):
        return r"10^{-5}"
    if np.isclose(rho, 1e5):
        return r"10^{5}"
    # For round/simple values prefer integer-like formatting, otherwise use
    # a compact decimal representation.
    if np.isclose(rho, round(rho)):
        return fr"{int(round(rho))}"
    return fr"{rho:g}"

def invert_teff_to_color(teff_array: np.ndarray) -> np.ndarray:
    """
    Invert:
      log10(Teff) = 3.999 - 0.654*C + 0.709*C^2 - 0.316*C^3
    for C = (G_BP - G_RP).

    Uses fsolve per element.
    """
    teff_array = np.asarray(teff_array, dtype=float)
    log_teff = np.log10(teff_array)

    def inversion_func(C, logT):
        return 3.999 - 0.654 * C + 0.709 * C**2 - 0.316 * C**3 - logT

    colors = np.empty_like(log_teff)
    for i, lt in enumerate(log_teff):
        colors[i] = fsolve(inversion_func, x0=1.0, args=(lt,))[0]
    return colors


# Star region boxes (T vs logg)
def build_star_regions() -> Tuple[Dict[str, Tuple[Tuple[float, float], Tuple[float, float]]], Dict[str, str]]:
    spectral_data = {
        "O-type": {"T": (33000, 50000), "M": [16, 100], "R": [6.6, 20]},
        "B-type": {"T": (10000, 33000), "M": [2.1, 16], "R": [1.8, 6.6]},
        "A-type": {"T": (7300, 10000), "M": [1.4, 2.1], "R": [1.4, 1.8]},
        "F-type": {"T": (6000, 7300), "M": [1.04, 1.4], "R": [1.15, 1.4]},
        "G-type": {"T": (5300, 6000), "M": [0.8, 1.04], "R": [0.96, 1.15]},
        "K-type": {"T": (3900, 5300), "M": [0.45, 0.8], "R": [0.7, 0.96]},
        "M-type": {"T": (2300, 3900), "M": [0.08, 0.45], "R": [0.1, 0.7]},
    }

    star_regions = {
        sp: (
            props["T"],
            (
                round(logg_sun + np.log10(props["M"][0]) - 2 * np.log10(props["R"][1]), 2),
                round(logg_sun + np.log10(props["M"][1]) - 2 * np.log10(props["R"][0]), 2),
            ),
        )
        for sp, props in spectral_data.items()
    }

    star_regions.update({
        "Brown Dwarfs": ((250, 3000), (4.5, 5.5)),
        "White Dwarfs": ((5500, 40000), (7.5, 9.0)),
        "Red Giants": ((2500, 5000), (0.0, 2.5)),
    })

    region_colors = {
        "O-type": "darkviolet", "B-type": "royalblue", "A-type": "dodgerblue",
        "F-type": "deepskyblue", "G-type": "gold", "K-type": "orange", "M-type": "red",
        "Brown Dwarfs": "sienna", "White Dwarfs": "lightgray", "Red Giants": "peachpuff",
    }
    return star_regions, region_colors


# Contours loading/trimming
def load_contours(contour_path: Path) -> Dict[float, List[List[np.ndarray]]]:
    """
    JSON structure: {rho_frac_str: [contour, contour, ...]}
    contour = [xis, rho_cs, T_cs, L_photos, R_photos, T_photos, rho_photos, m_photos, r_rads,
               L_ratios, R_max, M_max, tau_center]
    Each element is list-like -> convert to numpy arrays.
    """
    with contour_path.open("r", encoding="utf-8") as f:
        raw = json.load(f)

    contour_dict: Dict[float, List[List[np.ndarray]]] = {}
    for k, v in raw.items():
        rho_frac = float(k)
        contour_dict[rho_frac] = [[np.array(prop) for prop in contour] for contour in v]
    return contour_dict


def trim_contours_tau(contour_dict: Dict[float, List[List[np.ndarray]]], tau_min: float = 10.0) -> Dict[float, List[List[np.ndarray]]]:
    """
    Keep only indices where (10**tau_center) > tau_min
    """
    out: Dict[float, List[List[np.ndarray]]] = {}
    for rho_frac, contours in contour_dict.items():
        for contour in contours:
            tau_center_log10 = contour[12]
            condition = (10 ** tau_center_log10) > tau_min
            trimmed = [arr[condition] for arr in contour]
            out.setdefault(rho_frac, []).append(trimmed)
    return out


# Spectra grid interpolation
class Grid:
    def __init__(self, mh_input: float, teff_input: float, logg_input: float, set_type: str,
                 grids_dir: Path, spectra_root: Path):
        self.mh_input = float(mh_input)
        self.teff_input = float(teff_input)
        self.logg_input = float(logg_input)
        self.set_type = str(set_type)
        self.grids_dir = grids_dir
        self.spectra_root = spectra_root

    def get_surrounding_models(self, n: int = 6):
        mh_grid = np.loadtxt(self.grids_dir / "grid_mh.txt")
        teff_grid = np.loadtxt(self.grids_dir / "grid_teff.txt")
        logg_grid = np.loadtxt(self.grids_dir / "grid_logg.txt")

        mh_near = np.sort(mh_grid[np.argsort(np.abs(mh_grid - self.mh_input))[:n]])
        teff_near = np.sort(teff_grid[np.argsort(np.abs(teff_grid - self.teff_input))[:n]])
        logg_near = np.sort(logg_grid[np.argsort(np.abs(logg_grid - self.logg_input))[:n]])
        return mh_near, teff_near, logg_near

    def read_disk_integrated_spectra(self, mh: float, teff: float, logg: float):
        teff_i = int(round(teff))
        file_name = self.spectra_root / self.set_type / f"MH{mh}" / f"teff{teff_i}" / f"logg{logg}" / "mpsa_flux_spectra.dat"
        if not file_name.exists():
            raise FileNotFoundError(f"Missing spectra file: {file_name}")
        return np.loadtxt(file_name, skiprows=1, unpack=True)

    def get_interpolated_disk_integrated_spectrum(self):
        mh_near, teff_near, logg_near = self.get_surrounding_models()

        flux_grid = {}  # wavelength -> 3D grid (teff, logg, mh)

        for mh in mh_near:
            for teff in teff_near:
                for logg in logg_near:
                    wln, flux = self.read_disk_integrated_spectra(mh, teff, logg)

                    teff_index = int(np.where(teff_near == teff)[0][0])
                    logg_index = int(np.where(logg_near == logg)[0][0])
                    mh_index = int(np.where(mh_near == mh)[0][0])

                    for i, wavelength in enumerate(wln):
                        if wavelength not in flux_grid:
                            flux_grid[wavelength] = np.zeros((len(teff_near), len(logg_near), len(mh_near)), dtype=float)
                        flux_grid[wavelength][teff_index, logg_index, mh_index] = flux[i]

        # Ensure deterministic ordering
        wln_sorted = np.array(sorted(flux_grid.keys()), dtype=float)

        flux_interp_all = np.empty_like(wln_sorted, dtype=float)
        for j, wavelength in enumerate(wln_sorted):
            vals = flux_grid[wavelength]
            interp_val = interpn(
                (teff_near, logg_near, mh_near),
                vals,
                (self.teff_input, self.logg_input, self.mh_input),
                method="linear",
                bounds_error=False,
                fill_value=None,
            )
            flux_interp_all[j] = interp_val[0]

        return wln_sorted, flux_interp_all


# Plot parts
def process_spectra_examples(
    data_root: Path,
    plots_root: Path,
    set_type: str = "set2",
    mh_input: float = -0.135,
    spectra_root: Path | None = None,
) -> None:
    """
    Reproduces the "interpolated vs nearest model" plot(s) for the example nuggets.
    """
    # Use nuggets from the repository `profiles/` directory under optically_thick
    nuggets_dir = Path(__file__).resolve().parents[1] / "profiles"
    grids_dir = data_root / "grids"
    # Spectra root can be supplied explicitly (CLI); otherwise default to <data-root>/spectra.
    spectra_root = spectra_root or (data_root / "spectra")

    # Sanity check that the expected set folder exists.
    set_dir = spectra_root / set_type
    if not set_dir.exists():
        raise FileNotFoundError(
            "Could not find spectra set folder: "
            f"{set_dir}\n"
            "If you didn't extract set2 under <data-root>/spectra/, pass the correct location via --spectra-root.\n"
            "Expected structure: <spectra-root>/<set-type>/MH.../teff.../logg.../mpsa_flux_spectra.dat"
        )


    # Expect convective.nugget and radiative.nugget in the `profiles/` folder
    nugget_files = [nuggets_dir / "convective.nugget", nuggets_dir / "radiative.nugget"]
    for f in nugget_files:
        if not f.exists():
            raise FileNotFoundError(f"Missing nugget file in profiles/: {f}")

    nuggets = [load_nugget_file(str(f)) for f in nugget_files]

    out_folder = plots_root / "spectra_examples"
    ensure_dir(out_folder)

    for nug in nuggets:
        if not (hasattr(nug, "has_photosphere") and nug.has_photosphere()):
            continue

        teff_input = float(nug.T_photo())
        R_input_m = float(nug.r_photo())   # assume meters
        xi_input = float(nug.xi())
        rho_c_input = float(nug.rho_c())
        T_c_input = float(nug.T_c())

        # Surface gravity: use precomputed gravity from the nugget (log10 in cgs)
        logg_input = float(nug.log_g_photo(cgs=True))

        # Interpolate spectra
        grid = Grid(mh_input, teff_input, logg_input, set_type, grids_dir=grids_dir, spectra_root=spectra_root)
        wavelength_interp, flux_interp = grid.get_interpolated_disk_integrated_spectrum()

        # Nearest model via helper classes
        nearest = GridOfStellarParameters(mh_input, teff_input, logg_input)
        mh, teff, logg = nearest.get_closest_mh_teff_logg()
        print(f"NEAREST: Teff={teff}, logg={logg}, M/H={mh}")

        data = ReadData(spectra_root=spectra_root)
        wln, spectra = data.read_disk_integrated_spectra(mh, teff, logg, set_type)

        wln = np.array(wln).flatten()
        spectra = np.array(spectra).flatten()

        # Map interpolated flux onto the nearest-model wavelength grid
        flux_on_wln = np.interp(wln, wavelength_interp, flux_interp)

        plt.figure(figsize=(10, 8))
        plt.plot(wavelength_interp, flux_interp, label="Interpolated", linewidth=1.2)
        #plt.plot(wln, spectra, label="Nearest model", linewidth=1.0, alpha=0.8)

        plt.ylabel(r"Flux (erg/cm$^2$/s/nm)")
        plt.xlabel("Wavelength (nm)")
        plt.title(rf"$T_{{\mathrm{{eff}}}}={teff_input:.0f}\,$K, $\log g={logg_input:.2f}$, $[M/H]={mh_input:.3f}$")
        plt.grid(True, alpha=0.4)
        plt.xlim(0, 10000)
        #plt.legend()

        fname = out_folder / f"xi={xi_input:.0e}_rho={rho_c_input:.3e}_Tc={T_c_input:.3e}.png"
        plt.tight_layout()
        plt.savefig(fname, dpi=300, bbox_inches="tight")
        plt.close()


def combined_plot(
    contour_dict_trimmed: Dict[float, List[List[np.ndarray]]],
    rho_arr: List[float],
    plots_root: Path,
    luminosity_band: Tuple[float, float],
    filename: str,
) -> None:
    """
    Combined T_photo vs logg plot, with analytic luminosity contours, star-region boxes,
    and per-rho scatter points.

    luminosity_band is in log10(L/Lsun): e.g. (-9, -3) or (-3, 0)
    """
    star_regions, region_colors = build_star_regions()

    ensure_dir(plots_root / "logg_T")

    fig, ax = plt.subplots(figsize=(10, 8))
    colors_list = plt.cm.viridis(np.linspace(0, 1, len(rho_arr)))

    Tmin_points = []
    Tmax_points = []
    density_handles = []

    lo, hi = luminosity_band

    for color, rho_frac in zip(colors_list, rho_arr):
        if rho_frac not in contour_dict_trimmed:
            continue
        data = contour_dict_trimmed[rho_frac]
        if len(data) == 0:
            continue

        # concatenate contours within luminosity band
        properties = [[] for _ in data[0]]
        for contour in data:
            L_rel_raw = 10 ** contour[3]
            lumi = np.log10(L_rel_raw / L_sun_W)
            inds = (lumi > lo) & (lumi < hi)

            for i, prop in enumerate(contour):
                arr = np.array(prop)
                properties[i] += list(arr[inds])

        arrays = [np.array(p) for p in properties]
        if arrays[0].size == 0:
            continue

        T_photo = 10 ** arrays[5]  # K
        R_m = 10 ** arrays[4]      # m
        L_W = 10 ** arrays[3]      # W
        logg = arrays[13] + 2  # g_photo from contour_dict is log10(m/s^2), convert to cm/s^2 by adding log10(100)

        ax.scatter(T_photo, logg, color=color, s=20)

        min_idx = int(np.argmin(T_photo))
        max_idx = int(np.argmax(T_photo))
        Tmin_points.append((np.log10(T_photo[min_idx]), logg[min_idx]))
        Tmax_points.append((np.log10(T_photo[max_idx]), logg[max_idx]))

        # Legend entry: show only the numeric core-density label next to the
        # colored marker; the full legend title is added later.
        density_handles.append(
            Line2D([0], [0], marker="o", color="w",
                   label=rf"${format_rho(rho_frac)}$",
                   markerfacecolor=color, markersize=10)
        )

    if len(Tmin_points) < 2 or len(Tmax_points) < 2:
        raise RuntimeError("Not enough points to fit Tmin/Tmax boundary lines.")

    params_min = curve_fit(linear_fit, *map(np.array, zip(*Tmin_points)))[0]
    params_max = curve_fit(linear_fit, *map(np.array, zip(*Tmax_points)))[0]

    # boundary lines
    logT_space = np.linspace(0, 7, 100)
    ax.plot(10 ** logT_space, linear_fit(logT_space, *params_min), color="black")
    ax.plot(10 ** logT_space, linear_fit(logT_space, *params_max), color="black")

    # star region boxes
    for label, ((Tmin, Tmax), (gmin, gmax)) in star_regions.items():
        rect = Rectangle((Tmin, gmin), Tmax - Tmin, gmax - gmin,
                         facecolor=region_colors[label], edgecolor="k",
                         linewidth=0.4, alpha=0.5, label=label)
        ax.add_patch(rect)

    ax.set_xscale("log")
    ax.set_xlim(5e1, 5e5)
    ax.set_ylim(-0.5, 9.5)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_xlabel(r"$T_{\mathrm{photo}}$ (K)")
    ax.set_ylabel(r"log$_{10}$(g [cm/s$^2$])")

    def pow10_tex(exp: int) -> str:
        return "1" if exp == 0 else rf"10^{{{exp}}}"

    lo_exp = int(round(lo))
    hi_exp = int(round(hi))

    lumi_text = rf"$\log_{{10}}(L/L_\odot) \in [{lo_exp},\,{hi_exp}]$"
    ax.text(0.99, 0.35, lumi_text, transform=ax.transAxes, ha="right", va="top",
        fontsize=13, bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))

    rho_title = r"Mirror star core density $\rho_{\mathrm{core}}/\rho_{\mathrm{core},\odot}$"
    rho_legend = ax.legend(handles=density_handles, title=rho_title, loc="upper left",
                frameon=True, bbox_to_anchor=(0.02, 0.98), ncol=3, fontsize=11)
    rho_legend.get_frame().set_alpha(0.6)
    ax.add_artist(rho_legend)

    # Analytic luminosity contours (restricted between boundary lines)
    T_vals = np.logspace(1.5, 5.5, 300)
    g_vals = np.logspace(-0.5, 9.5, 300)
    TT, GG = np.meshgrid(T_vals, g_vals)

    logL_analytic = (
        np.log10((TT**4) * (GG**2))
        - 2 * np.log10(rho_c_sun)
        - np.log10(L_sun_cgs)
        + np.log10(9 * sigma_cgs / (4 * np.pi * G**2))
    )

    logT = np.log10(TT)
    logG = np.log10(GG)

    g_min_line = linear_fit(logT, *params_min)
    g_max_line = linear_fit(logT, *params_max)

    mask = (logG < g_max_line) | (logG > g_min_line)
    logL_masked = np.copy(logL_analytic)
    logL_masked[mask] = np.nan

    contour_levels = np.arange(np.floor(np.nanmin(logL_masked)), np.ceil(np.nanmax(logL_masked)), 1.0)
    cs = ax.contour(TT, logG, logL_masked, levels=contour_levels,
                    colors="black", linestyles="dashed", linewidths=0.9)
    label_texts = ax.clabel(cs, fmt="%d", fontsize=15, colors="black")
    for txt in label_texts:
        txt.set_zorder(30)

    contour_proxy = mlines.Line2D([], [], color="black", linestyle="dashed", linewidth=0.9)

    handles, labels = ax.get_legend_handles_labels()
    handles.append(contour_proxy)
    labels.append(r"$\log_{10}(L/L_\odot)$ contours" + "\n" + r"for $\rho_c/\rho_\odot = 1$")

    star_legend = ax.legend(handles, labels, title="Legend", loc="lower right", ncol=2, frameon=True, fontsize=11)
    star_legend.get_frame().set_alpha(0.6)
    try:
        star_legend.get_title().set_fontsize(10)
    except Exception:
        pass
    ax.add_artist(star_legend)

    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout(rect=[0.08, 0, 1, 0.90])
    outpath = plots_root / "logg_T" / filename
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


# Gaia cache
def load_gaia_hr_cache(gaia_cache_dir: Path) -> Tuple[np.ndarray, np.ndarray]:
    fp = gaia_cache_dir / "gaia_observed_hr.csv"
    if not fp.exists():
        raise FileNotFoundError(
            f"Missing Gaia cache {fp}."
        )
    df = pd.read_csv(fp)
    return df["M_G"].values, df["BP_RP"].values


def load_gaia_hr_logg_cache(gaia_cache_dir: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    fp = gaia_cache_dir / "gaia_observed_hr_logg.csv"
    if not fp.exists():
        raise FileNotFoundError(
            f"Missing Gaia cache {fp}."
        )
    df = pd.read_csv(fp)
    return df["M_G"].values, df["BP_RP"].values, df["logg"].values


def load_wd_hr_logg_cache(gaia_cache_dir: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    fp = gaia_cache_dir / "gaia_white_dwarfs_hr_logg.csv"
    if not fp.exists():
        raise FileNotFoundError(
            f"Missing WD cache {fp}."
        )
    df = pd.read_csv(fp)
    return df["M_G"].values, df["BP_RP"].values, df["logg_H"].values

# HR plots
def stacked_hr_plots(
    contour_dict_trimmed: Dict[float, List[List[np.ndarray]]],
    rho_arr: List[float],
    plots_root: Path,
    gaia_cache_dir: Path,
) -> None:
    ensure_dir(plots_root / "HR")

    m_g, bp_rp = load_gaia_hr_cache(gaia_cache_dir)

    fig, axes = plt.subplots(
        nrows=len(rho_arr), ncols=3,
        figsize=(20, 6 * len(rho_arr)),
        sharex=True, sharey=True
    )

    for (ax1, ax2, ax3), rho_frac in zip(axes, rho_arr):
        data = contour_dict_trimmed[rho_frac]
        properties = [[] for _ in data[0]]

        for contour in data:
            L_rel_raw = 10 ** contour[3]
            lumi = np.log10(L_rel_raw / L_sun_W)
            inds = (lumi > -9) & (lumi < 0)
            for i, prop in enumerate(contour):
                properties[i] += list(np.array(prop)[inds])

        arrays = [np.array(p) for p in properties]
        T_photo = 10 ** arrays[5]
        R_m = 10 ** arrays[4]
        m_photo_kg = 10 ** arrays[7]
        L_W = 10 ** arrays[3]
        xi = arrays[0]
        logg = arrays[13] + 2  # g_photo from contour_dict is log10(m/s^2), convert to cm/s^2 by adding log10(100)

        M_G = 4.83 - 2.5 * np.log10(L_W / L_sun_W)

        colors_hr = invert_teff_to_color(T_photo)
        logmass = np.log10(m_photo_kg * 1000.0)  # kg -> g
        logxi = np.log10(xi)

        # Col 1: color=logg (use discrete custom_cmap)
        Ncol = custom_cmap.N
        #norm1 = mpl.colors.BoundaryNorm(np.linspace(0, 9, Ncol + 1), Ncol)
        norm1 = mpl.colors.Normalize(0, 7+1/3)
        sc1 = ax1.scatter(colors_hr, M_G, c=logg, cmap=custom_cmap, norm=norm1, s=10, zorder=3)
        ax1.hexbin(bp_rp, m_g, gridsize=150, cmap="Greys", bins="log", zorder=1, alpha=1)
        ax1.annotate(rf"$\rho_c={format_rho(rho_frac)}\,\rho_\odot$", xy=(-0.1, 0.5),
                     xycoords="axes fraction", ha="right", va="center", rotation=90)

        # Col 2: color=logmass (discrete custom_cmap)
        #norm2 = mpl.colors.BoundaryNorm(np.linspace(11, 32, Ncol + 1), Ncol)
        norm2 = mpl.colors.Normalize(11, 33)
        sc2 = ax2.scatter(colors_hr, M_G, c=logmass, cmap=custom_cmap, norm=norm2, s=10, zorder=3)
        ax2.hexbin(bp_rp, m_g, gridsize=150, cmap="Greys", bins="log", zorder=1, alpha=1)

        # Col 3: color=logxi (discrete custom_cmap)
        #norm3 = mpl.colors.BoundaryNorm(np.linspace(-26, -16, Ncol + 1), Ncol)
        norm3 = mpl.colors.Normalize(-26, -15)
        sc3 = ax3.scatter(colors_hr, M_G, c=logxi, cmap=custom_cmap, norm=norm3, s=10, zorder=3)
        ax3.hexbin(bp_rp, m_g, gridsize=150, cmap="Greys", bins="log", zorder=1, alpha=1)

        for ax in (ax1, ax2, ax3):
            ax.invert_yaxis()
            ax.grid(True, alpha=0.4)
            ax.set_xlim(-1, 5.5)
            ax.set_ylim(28.5, -5)

    # Shared discrete colorbars using the repo's custom_cmap with nice integer ticks
    cbar_ax1 = fig.add_axes([0.125, 0.93, 0.22, 0.01])
    sm1 = mpl.cm.ScalarMappable(norm=norm1, cmap=custom_cmap)
    #bins1 = np.linspace(0, 9, Ncol + 1)
    cb1 = fig.colorbar(sm1, cax=cbar_ax1, orientation="horizontal", ticks=[0, 2, 4, 6], extend="max")
    cb1.set_label(r"log$_{10}(g$ [cm/s$^2$])")

    cbar_ax2 = fig.add_axes([0.425, 0.93, 0.22, 0.01])
    sm2 = mpl.cm.ScalarMappable(norm=norm2, cmap=custom_cmap)
    #bins2 = np.linspace(11, 32, Ncol + 1)
    cb2 = fig.colorbar(sm2, cax=cbar_ax2, orientation="horizontal", ticks=[13, 17, 21, 25, 29, 33], extend="max")
    cb2.set_label(r"log$_{10}(M$ [g])")

    cbar_ax3 = fig.add_axes([0.7, 0.93, 0.22, 0.01])
    sm3 = mpl.cm.ScalarMappable(norm=norm3, cmap=custom_cmap)
    #bins3 = np.linspace(-26, -16, Ncol + 1)
    cb3 = fig.colorbar(sm3, cax=cbar_ax3, orientation="horizontal", ticks=[-26, -22, -18, -16], extend="max")
    cb3.set_label(r"log$_{10}(\xi)$")

    fig.text(0.5, 0.04, r"$G_{\rm BP}-G_{\rm RP}$", ha="center")
    fig.text(0.04, 0.5, r"$M_G$", va="center", rotation="vertical")

    plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.91])
    outpath = plots_root / "HR" / "stacked_HR_threecol.png"
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()

def stacked_hr_plots_3d(
    contour_dict_trimmed: dict,
    rho_fracs: list,
    plots_root: Path,
    gaia_cache_dir: Path,
) -> None:
    """
    Single-column stacked 3D HR plot (BP-RP, M_G, logg) combining:
      - model points (logg computed from rho_c_sun * R_photo)
      - Gaia observed stars with logg_gspphot (from cache)
      - Gaia WD catalog with logg_H (from cache)

    Writes: plots_root/HR/stacked_HR_3D.png
    """
    # load caches (must exist)
    gaia_fp = gaia_cache_dir / "gaia_observed_hr_logg.csv"
    wd_fp = gaia_cache_dir / "gaia_white_dwarfs_hr_logg.csv"

    if not gaia_fp.exists():
        raise FileNotFoundError(f"Missing Gaia cache: {gaia_fp}")
    if not wd_fp.exists():
        raise FileNotFoundError(f"Missing WD cache: {wd_fp}")

    gaia_data = pd.read_csv(gaia_fp)
    wd_data = pd.read_csv(wd_fp)

    m_g = gaia_data["M_G"].values
    bp_rp = gaia_data["BP_RP"].values
    logg_gaia = gaia_data["logg"].values

    bp_rp_wd = wd_data["BP_RP"].values
    M_G_wd = wd_data["M_G"].values
    logg_wd = wd_data["logg_H"].values

    plot_folder = plots_root / "HR"
    ensure_dir(plot_folder)

    n = len(rho_fracs)
    nrows, ncols = 3, 3
    fig = plt.figure(figsize=(25, 20), dpi=300)
    
    plot_positions = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 1)]  # 7 plots
    axes_list = []
    for row, col in plot_positions:
        ax = fig.add_subplot(nrows, ncols, row * ncols + col + 1, projection="3d")
        axes_list.append(ax)
        
    for ax_idx, (ax, rho_frac) in enumerate(zip(axes_list, rho_fracs)):
        data = contour_dict_trimmed[rho_frac]
        properties = [[] for _ in data[0]]

        # Collect model properties in luminosity window
        for contour in data:
            L_rel_raw = 10 ** contour[3]
            lumi = np.log10(L_rel_raw / L_sun_W)
            inds = (lumi > -9) & (lumi < 0)
            for i, prop in enumerate(contour):
                properties[i] += list(np.array(prop)[inds])

        arrays = [np.array(p) for p in properties]
        T_photo = 10 ** arrays[5]   # K
        R_m = 10 ** arrays[4]       # m
        L_W = 10 ** arrays[3]       # W
        logg_model = arrays[13] + 2  # g_photo from contour_dict is log10(m/s^2), convert to cm/s^2 by adding log10(100)

        M_G_model = 4.83 - 2.5 * np.log10(L_W / L_sun_W)
        colors_hr = invert_teff_to_color(T_photo)

        # Combine all points for the 3D scatter
        x_all = np.concatenate([colors_hr, bp_rp, bp_rp_wd])
        y_all = np.concatenate([M_G_model, m_g, M_G_wd])
        z_all = np.concatenate([logg_model, logg_gaia, logg_wd])

        # Discrete 3D scatter using the repository custom_cmap
        Ncol = custom_cmap.N
        #norm_main = mpl.colors.BoundaryNorm(np.linspace(0, 9, Ncol + 1), Ncol)
        norm_main = mpl.colors.Normalize(0, 7+1/3)
        sc = ax.scatter(
            x_all, y_all, z_all,
            c=z_all, cmap=custom_cmap, norm=norm_main,
            s=6
        )

        # "Legs" + shadows
        step_model = 50
        step_gaia = 5000
        z_shadow_plane = 0.0
        # Gaia legs
        for i in range(0, len(bp_rp), step_gaia):
            x, y, z = bp_rp[i], m_g[i], logg_gaia[i]
            ax.plot([x, x], [y, y], [z_shadow_plane, z],
                    color="gray", alpha=0.25, linewidth=0.4, zorder=1)
        ax.scatter(bp_rp, m_g, z_shadow_plane, color="gray", s=2, alpha=0.15, zorder=1)

        # Model legs
        for i in range(0, len(colors_hr), step_model):
            x, y, z = colors_hr[i], M_G_model[i], logg_model[i]
            ax.plot([x, x], [y, y], [z_shadow_plane, z],
                    color="black", alpha=0.4, linewidth=0.5, zorder=2)
        ax.scatter(colors_hr, M_G_model, z_shadow_plane, color="black", s=3, zorder=2)

        # WD legs
        for i in range(0, len(bp_rp_wd), step_gaia):
            x, y, z = bp_rp_wd[i], M_G_wd[i], logg_wd[i]
            ax.plot([x, x], [y, y], [z_shadow_plane, z],
                    color="gray", alpha=0.25, linewidth=0.4, zorder=1)
        ax.scatter(bp_rp_wd, M_G_wd, z_shadow_plane, color="gray", s=2, alpha=0.15, zorder=1)
        # Labels/orientation
        ax.set_xlabel(r"$G_{\rm BP} - G_{\rm RP}$", labelpad=10)
        ax.set_ylabel(r"$M_G$", labelpad=10)
        ax.set_zlabel(r"log$_{10}(g$ [cm/s$^2$])", labelpad=10)

        ax.invert_yaxis()
        ax.view_init(elev=25, azim=-30)
        ax.set_xlim(-1, 5.5)
        ax.set_ylim(28.5, -5)
        ax.set_zticks([0, 5, 10])

        ax.text2D(
            -0.06, 0.5, rf"$\rho_c = {format_rho(rho_frac)}\,\rho_\odot$",
            transform=ax.transAxes,
            fontweight="bold",
            ha="right", va="center",
            rotation=90,
        )
    

    cbar_slot = fig.add_subplot(nrows, ncols, 2 * ncols + 0 + 1)
    cbar_slot.set_axis_off() 

    cbar_ax = inset_axes(
        cbar_slot,
        width="95%",  
        height="12%",  
        loc="lower left",
        borderpad=0.6,
    )

    sm_main = mpl.cm.ScalarMappable(norm=norm_main, cmap=custom_cmap)
    cb = fig.colorbar(
        sm_main,
        cax=cbar_ax,
        orientation="horizontal",
        ticks=[0, 2, 4, 6],
        extend="max"
    )

    cb.set_label(r"log$_{10}(g$ [cm/s$^2$])", labelpad=2)
    cb.ax.tick_params(labelsize=20, pad=1)
    cb.ax.xaxis.set_label_position("bottom")

    plt.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.1, hspace=0.2, wspace=0.1)
    plt.savefig(plot_folder / "stacked_HR_3D.png", dpi=300, bbox_inches="tight", pad_inches=0.45)
    plt.close()


# Main
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root",
        type=str,
        default="em_signatures/data",
        help="Root directory containing grids/, contours/, nuggets/, spectra/, gaia_cache/, etc."
    )
    parser.add_argument(
        "--plots-root",
        type=str,
        default="em_signatures/plots",
        help="Where to write plots (will be created)."
    )
    parser.add_argument("--set-type", type=str, default="set2", choices=["set1", "set2"])
    parser.add_argument(
        "--spectra-root",
        type=str,
        default=None,
        help=(
            "Optional path to the spectra root folder (expects set1/ or set2/ beneath it). "
            "If omitted, defaults to <data-root>/spectra."
        ),
    )
    parser.add_argument("--skip-spectra-examples", action="store_true")
    parser.add_argument("--skip-loggT", action="store_true")
    parser.add_argument("--skip-hr", action="store_true")
    parser.add_argument("--skip-hr-3d", action="store_true")
    args = parser.parse_args()

    data_root = Path(args.data_root).expanduser().resolve()
    plots_root = Path(args.plots_root).expanduser().resolve()
    ensure_dir(plots_root)

    # Optional spectra root override
    spectra_root = None
    if args.spectra_root:
        spectra_root = Path(args.spectra_root).expanduser().resolve()

    # Data subdirs
    gaia_cache_dir = data_root / "gaia_cache"

    # Keep as defaults but read from data_root
    rho_arr = [1e-5, 0.01, 0.1, 1, 10, 100, 1e5]

    # Load contour file from paramspace
    contour_path = Path(__file__).parent.parent / "paramspace" / "contour_dict.json"

    contour_dict = load_contours(contour_path)
    contour_dict_trimmed = trim_contours_tau(contour_dict, tau_min=10.0)

    if not args.skip_spectra_examples:
        process_spectra_examples(
            data_root=data_root,
            plots_root=plots_root,
            set_type=args.set_type,
            spectra_root=spectra_root,
        )

    if not args.skip_loggT:
        combined_plot(
            contour_dict_trimmed=contour_dict_trimmed,
            rho_arr=rho_arr,
            plots_root=plots_root,
            luminosity_band=(-9, -3),
            filename="combined_lowL.png",
        )
        combined_plot(
            contour_dict_trimmed=contour_dict_trimmed,
            rho_arr=rho_arr,
            plots_root=plots_root,
            luminosity_band=(-3, 0),
            filename="combined.png",
        )

    if not args.skip_hr:
        stacked_hr_plots(
            contour_dict_trimmed=contour_dict_trimmed,
            rho_arr=rho_arr,
            plots_root=plots_root,
            gaia_cache_dir=gaia_cache_dir,
        )
    
    if not args.skip_hr_3d:
        stacked_hr_plots_3d(contour_dict_trimmed, rho_arr, plots_root, gaia_cache_dir)

    print(f"Done. Plots written under: {plots_root}")


if __name__ == "__main__":
    main()
