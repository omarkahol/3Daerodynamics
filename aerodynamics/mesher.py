import numpy as np
from math import radians, cos, sin
from typing import Dict, List, Any, Tuple
from .airfoil import AeroFoil
from .panel import VortexPanel

class WingMesher:
    """Generates a 3D panel mesh grid for a wing based on station-based multi-section configuration."""

    def __init__(self, config: Dict[str, Any], data_folder: str = "airfoilDATA"):
        self.config = config
        self.wing_cfg = config["wing"]
        self.mesh_cfg = config["meshing"]
        self.data_folder = data_folder

        # Parse stations sorted by spanwise coordinate y (root to tip)
        self.stations = sorted(self.wing_cfg["stations"], key=lambda st: float(st["y"]))
        
        # Total wingspan is double the y-coordinate of the outboard tip station
        self.span = 2.0 * float(self.stations[-1]["y"])
        
        self.n_span = int(self.mesh_cfg["n_span"])
        self.n_chord = int(self.mesh_cfg["n_chord"])

        # Initialize and download required airfoils
        self.airfoils = {}
        self.camber_fits = {}
        for st in self.stations:
            name = st.get("airfoil_name")
            if name and name.lower() != "none":
                name = name.lower()
                if name not in self.airfoils:
                    foil = AeroFoil(name, self.data_folder)
                    self.airfoils[name] = foil
                    self.camber_fits[name] = foil.fitCamberLine(4)

    def get_chord(self, y: float) -> float:
        """Returns the interpolated local chord at spanwise location y."""
        y_abs = abs(y)
        st1, st2 = self.stations[0], self.stations[-1]
        for idx in range(len(self.stations) - 1):
            if self.stations[idx]["y"] <= y_abs <= self.stations[idx + 1]["y"]:
                st1 = self.stations[idx]
                st2 = self.stations[idx + 1]
                break

        dy = float(st2["y"]) - float(st1["y"])
        eta = (y_abs - float(st1["y"])) / dy if dy > 1e-12 else 0.0
        return float((1.0 - eta) * float(st1["chord"]) + eta * float(st2["chord"]))

    def _interpolate_station_properties(self, y: float, s: float) -> Tuple[float, float, float, float, float]:
        """Interpolates geometry properties at a given spanwise y and chordwise fraction s.
        
        Returns:
            chord, x_le, z_le, twist_rad, z_camber
        """
        y_abs = abs(y)
        
        # Find the bounding stations
        st1, st2 = self.stations[0], self.stations[-1]
        for idx in range(len(self.stations) - 1):
            if self.stations[idx]["y"] <= y_abs <= self.stations[idx + 1]["y"]:
                st1 = self.stations[idx]
                st2 = self.stations[idx + 1]
                break

        dy = float(st2["y"]) - float(st1["y"])
        eta = (y_abs - float(st1["y"])) / dy if dy > 1e-12 else 0.0

        # Linear interpolation of properties
        chord = (1.0 - eta) * float(st1["chord"]) + eta * float(st2["chord"])
        x_le = (1.0 - eta) * float(st1["x_le"]) + eta * float(st2["x_le"])
        z_le = (1.0 - eta) * float(st1.get("z_le", 0.0)) + eta * float(st2.get("z_le", 0.0))
        twist_deg = (1.0 - eta) * float(st1.get("twist_deg", 0.0)) + eta * float(st2.get("twist_deg", 0.0))
        twist_rad = radians(twist_deg)

        # Linear interpolation of camber lines
        name1 = st1.get("airfoil_name")
        name2 = st2.get("airfoil_name")
        camber1 = self.camber_fits[name1.lower()](s) if (name1 and name1.lower() in self.camber_fits) else 0.0
        camber2 = self.camber_fits[name2.lower()](s) if (name2 and name2.lower() in self.camber_fits) else 0.0
        z_camber = chord * ((1.0 - eta) * camber1 + eta * camber2)

        return chord, x_le, z_le, twist_rad, z_camber

    def generate_mesh(self) -> List[VortexPanel]:
        """Generates the mesh grid and returns a list of VortexPanel objects."""
        y_nodes = np.linspace(-0.5 * self.span, 0.5 * self.span, self.n_span)
        s_nodes = np.linspace(0.0, 1.0, self.n_chord)

        grid = np.zeros((self.n_span, self.n_chord, 3))

        for i, y in enumerate(y_nodes):
            for j, s in enumerate(s_nodes):
                chord, x_le, z_le, twist_rad, z_camber = self._interpolate_station_properties(y, s)
                
                dx = s * chord
                dz = z_camber

                # Apply twist rotation (pitching chordwise grid nodes about the leading edge)
                x_rot = dx * cos(twist_rad) - dz * sin(twist_rad)
                z_rot = dx * sin(twist_rad) + dz * cos(twist_rad)

                grid[i, j, 0] = x_le + x_rot
                grid[i, j, 1] = y
                grid[i, j, 2] = z_le + z_rot

        # Generate panels from grid
        panels = []
        for i in range(self.n_span - 1):
            for j in range(self.n_chord - 1):
                p1 = tuple(grid[i, j])
                p2 = tuple(grid[i + 1, j])
                p3 = tuple(grid[i, j + 1])
                p4 = tuple(grid[i + 1, j + 1])
                
                panel = VortexPanel(p1, p2, p3, p4)
                panels.append(panel)

        return panels
