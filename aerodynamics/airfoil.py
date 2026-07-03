import urllib.request
import os
from pathlib import Path
import numpy as np
from typing import Callable

class AeroFoil:
    """Class to download, parse, and fit camber line and aerodynamic performance data of airfoils."""

    def __init__(self, airfoil_name: str, data_folder: str):
        self.airfoil_name = airfoil_name.lower()
        self.data_folder = Path(data_folder)
        self.data_folder.mkdir(parents=True, exist_ok=True)
        
        self.airfoil_file = self.data_folder / f"{self.airfoil_name}_coordinates.dat"
        self.airfoil_url = f"https://m-selig.ae.illinois.edu/ads/coord_updates/{airfoil_name.upper()}.DAT"

        if not self.airfoil_file.is_file():
            print(f"Downloading airfoil coordinates for {self.airfoil_name}...")
            try:
                urllib.request.urlretrieve(self.airfoil_url, self.airfoil_file)
            except Exception:
                fallback_url = f"https://m-selig.ae.illinois.edu/ads/coord/{self.airfoil_name}.dat"
                urllib.request.urlretrieve(fallback_url, self.airfoil_file)
                
        self.x_base = np.array([])
        self.y_base = np.array([])
        self._read_coords()

    def _read_coords(self) -> None:
        x_coords = []
        y_coords = []
        with open(self.airfoil_file, "r") as f:
            for line in f:
                values = line.split()
                if len(values) < 2:
                    continue
                try:
                    # Validate they are numeric
                    float(values[0])
                    float(values[1])
                except ValueError:
                    continue
                
                val_x = float(values[0])
                val_y = float(values[1])
                # Skip the Selig header (e.g., "66 66" or "80 79") which specifies number of points
                if val_x > 1.05 and val_y > 1.05 and abs(val_x - val_y) < 2.0:
                    continue
                x_coords.append(val_x)
                y_coords.append(val_y)
                     
        x_arr = np.array(x_coords)
        y_arr = np.array(y_coords)
        
        # Scale coordinates if they are in percentage format (0 to 100)
        if len(x_arr) > 0 and np.max(x_arr) > 1.05:
            x_arr = x_arr / 100.0
            y_arr = y_arr / 100.0
            
        self.x_base = x_arr
        self.y_base = y_arr

    def fitCamberLine(self, degree: int = 4) -> Callable[[float], float]:
        """Fits the mean camber line of the airfoil using a Cubic Spline."""
        from scipy.interpolate import interp1d, CubicSpline
        
        n_pts = len(self.x_base)
        
        # Auto-detect format and split upper/lower surfaces
        if self.x_base[0] > 0.5:
            # Selig format: TE upper -> LE -> TE lower
            le_idx = np.argmin(self.x_base)
            x_upper = self.x_base[:le_idx+1]
            y_upper = self.y_base[:le_idx+1]
            x_lower = self.x_base[le_idx:]
            y_lower = self.y_base[le_idx:]
        else:
            # Lednicer format: LE -> TE upper, then LE -> TE lower
            mid_idx = n_pts // 2
            x_upper = self.x_base[:mid_idx]
            y_upper = self.y_base[:mid_idx]
            x_lower = self.x_base[mid_idx:]
            y_lower = self.y_base[mid_idx:]
        
        # Sort and make unique for upper
        u_sort = np.argsort(x_upper)
        xu_u, u_uniq = np.unique(x_upper[u_sort], return_index=True)
        yu_u = y_upper[u_sort][u_uniq]
        
        # Sort and make unique for lower
        l_sort = np.argsort(x_lower)
        xl_u, l_uniq = np.unique(x_lower[l_sort], return_index=True)
        yl_l = y_lower[l_sort][l_uniq]
        
        # Define common grid from 0 to 1
        x_common = np.linspace(0.0, 1.0, 100)
        
        # Interpolate both surfaces onto the common grid
        yu_interp = interp1d(xu_u, yu_u, kind="linear", fill_value="extrapolate")(x_common)
        yl_interp = interp1d(xl_u, yl_l, kind="linear", fill_value="extrapolate")(x_common)
        
        # Compute mean camber coordinates
        y_camber = 0.5 * (yu_interp + yl_interp)
        
        # Fit Cubic Spline on the smooth camber line coordinates
        spline = CubicSpline(x_common, y_camber)
        return lambda x: float(spline(x))

    def downloadPerformance(self, re_target: int = 1000000) -> None:
        """Downloads the airfoil polar performance data from Airfoil Tools for the target Reynolds number."""
        self.url_string_perf = f"http://airfoiltools.com/polar/text?polar=xf-{self.airfoil_name}-il-{re_target}"
        self.perf_file_name = self.data_folder / f"{self.airfoil_name}_performance_{re_target}.txt"

        if not self.perf_file_name.is_file():
            print(f"Downloading airfoil performance polar for {self.airfoil_name} at Re={re_target}...")
            urllib.request.urlretrieve(self.url_string_perf, self.perf_file_name)
        self._read_perf()

    def _read_perf(self) -> None:
        self.alfa = []
        self.CL = []
        self.CD = []
        with open(self.perf_file_name, "r") as f:
            for line in f:
                values = line.split()
                numbers = [el for el in values if len(el) > 3]
                if len(numbers) < 3:
                    continue
                try:
                    float(numbers[0])
                    float(numbers[1])
                    float(numbers[2])
                except ValueError:
                    continue
                self.alfa.append(float(numbers[0]))
                self.CL.append(float(numbers[1]))
                self.CD.append(float(numbers[2]))
        
        self.alfa = np.array(self.alfa)
        self.CL = np.array(self.CL)
        self.CD = np.array(self.CD)

    def fitCL(self, degree: int = 5, convert_radians: bool = True) -> Callable[[float], float]:
        """Fits the CL vs alpha curve of the airfoil using a Cubic Spline."""
        from scipy.interpolate import CubicSpline
        alfa_fit = np.radians(self.alfa) if convert_radians else self.alfa
        
        # Sort and unique
        sort_idx = np.argsort(alfa_fit)
        alfa_sorted = alfa_fit[sort_idx]
        cl_sorted = self.CL[sort_idx]
        
        alfa_unique, unique_idx = np.unique(alfa_sorted, return_index=True)
        cl_unique = cl_sorted[unique_idx]
        
        spline = CubicSpline(alfa_unique, cl_unique, extrapolate=True)
        return lambda x: float(spline(x))

    def fitCD(self, degree: int = 5, convert_radians: bool = True) -> Callable[[float], float]:
        """Fits the CD vs alpha curve of the airfoil using a Cubic Spline."""
        from scipy.interpolate import CubicSpline
        alfa_fit = np.radians(self.alfa) if convert_radians else self.alfa
        
        # Sort and unique
        sort_idx = np.argsort(alfa_fit)
        alfa_sorted = alfa_fit[sort_idx]
        cd_sorted = self.CD[sort_idx]
        
        alfa_unique, unique_idx = np.unique(alfa_sorted, return_index=True)
        cd_unique = cd_sorted[unique_idx]
        
        spline = CubicSpline(alfa_unique, cd_unique, extrapolate=True)
        return lambda x: max(0.0, float(spline(x)))
