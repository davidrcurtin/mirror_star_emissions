import numpy as np
from pathlib import Path

class ReadData:
    def __init__(self, spectra_root=None):
        # spectra_root should be the folder that contains set1/ and set2/
        self.spectra_root = Path(spectra_root) if spectra_root is not None else Path(".")

    def _p(self, rel_path: str) -> str:
        return str(self.spectra_root / rel_path)

    def read_model_atmosphere(self, mh, teff, logg, set_type):
        file_name = self._p(f"{set_type}/MH{mh}/teff{teff}/logg{logg}/mpsa_model_atmosphere.dat")
        return np.genfromtxt(file_name, skip_header=2, skip_footer=23)

    def read_clv_spectra(self, mh, teff, logg, set_type):
        file_name = self._p(f"{set_type}/MH{mh}/teff{teff}/logg{logg}/mpsa_intensity_spectra.dat")
        data = np.loadtxt(file_name, skiprows=2)
        return data[:, 0], data[:, 1:]

    def read_disk_integrated_spectra(self, mh, teff, logg, set_type):
        file_name = self._p(f"{set_type}/MH{mh}/teff{teff}/logg{logg}/mpsa_flux_spectra.dat")
        return np.loadtxt(file_name, skiprows=1, unpack=True)

    def read_mu_positions(self, mh, teff, logg, set_type):
        file_name = self._p(f"{set_type}/MH{mh}/teff{teff}/logg{logg}/mpsa_intensity_spectra.dat")
        f = open(file_name, "r")
        data = f.readlines()
        muval = data[1].split()[2:]
        return np.array(muval).astype(float)
