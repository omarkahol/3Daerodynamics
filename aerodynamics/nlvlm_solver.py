import numpy as np
from math import radians, degrees, pi
from typing import List, Dict, Any
from .panel import VortexPanel
from .vlm_solver import VortexLatticeSolver
from .gamma_tool import VortexLineTool
from .airfoil import AeroFoil

# ANSI Escape Codes for CLI styling
RESET = "\033[0m"
BLUE = "\033[34m"

class NonLinearVortexLatticeSolver:
    """Coupled Non-linear Vortex Lattice Method (NL-VLM) solver using decambering method."""

    def __init__(self, panels: List[VortexPanel], mesher: Any):
        self.panels = panels
        self.mesher = mesher
        self.n_span_sections = mesher.n_span - 1
        self.n_chord_sections = mesher.n_chord - 1
        self.nPoints = len(panels)
        
        # Instantiate a standard linear VLM solver
        self.vlm = VortexLatticeSolver(panels)

        # Initialize decambering angles for each panel (default to 0)
        self.decamber_angles = np.zeros(self.nPoints)

    def _get_airfoil_lift(self, y: float, alfa_eff: float) -> float:
        """Finds target 2D Cl for spanwise position y by interpolating between station airfoils."""
        y_abs = abs(y)
        st1_idx, st2_idx = 0, len(self.mesher.stations) - 1
        for idx in range(len(self.mesher.stations) - 1):
            if self.mesher.stations[idx]["y"] <= y_abs <= self.mesher.stations[idx + 1]["y"]:
                st1_idx = idx
                st2_idx = idx + 1
                break

        st1 = self.mesher.stations[st1_idx]
        st2 = self.mesher.stations[st2_idx]
        
        dy = float(st2["y"]) - float(st1["y"])
        eta = (y_abs - float(st1["y"])) / dy if dy > 1e-12 else 0.0

        cl1 = 0.0
        af1 = self.station_airfoils[st1_idx]
        if af1 is not None:
            cl1 = af1.fitCL(5)(alfa_eff)
            
        cl2 = 0.0
        af2 = self.station_airfoils[st2_idx]
        if af2 is not None:
            cl2 = af2.fitCL(5)(alfa_eff)

        return (1.0 - eta) * cl1 + eta * cl2

    def _get_airfoil_drag(self, y: float, alfa_eff: float) -> float:
        """Finds target 2D Cd for spanwise position y by interpolating between station airfoils."""
        y_abs = abs(y)
        st1_idx, st2_idx = 0, len(self.mesher.stations) - 1
        for idx in range(len(self.mesher.stations) - 1):
            if self.mesher.stations[idx]["y"] <= y_abs <= self.mesher.stations[idx + 1]["y"]:
                st1_idx = idx
                st2_idx = idx + 1
                break

        st1 = self.mesher.stations[st1_idx]
        st2 = self.mesher.stations[st2_idx]

        dy = float(st2["y"]) - float(st1["y"])
        eta = (y_abs - float(st1["y"])) / dy if dy > 1e-12 else 0.0

        cd1 = 0.0
        af1 = self.station_airfoils[st1_idx]
        if af1 is not None:
            cd1 = af1.fitCD(5)(alfa_eff)
            
        cd2 = 0.0
        af2 = self.station_airfoils[st2_idx]
        if af2 is not None:
            cd2 = af2.fitCD(5)(alfa_eff)

        return (1.0 - eta) * cd1 + eta * cd2

    def solve(self, alfa_wing_deg: float, v_inf: float, damping: float = 0.05, 
              tol: float = 1e-3, max_iter: int = 150) -> tuple[np.ndarray, np.ndarray, int, float, float, float, float, float, float]:
        """Runs the iterative decambering NL-VLM loop.
        
        Returns:
            control_points: spanwise control point coordinates
            leading_edge_gamma: spanwise circulation distribution
            iterations: count of iterations run
            error: final residual error
            CL: Lift coefficient
            CDi: Induced drag coefficient
            CDp: Profile drag coefficient
            CD: Total drag coefficient
            efficiency: Oswald efficiency factor
        """
        alfa_wing = radians(alfa_wing_deg)
        self.vlm.setDynamicData(alfa_wing_deg, v_inf)
        
        # Dynamically load and cache station airfoils at their correct Reynolds numbers
        self.station_airfoils = []
        print(f"\n{BLUE}[NL-VLM]{RESET} Resolving Station Reynolds Numbers (Vinf = {v_inf:.2f} m/s):")
        for idx, st in enumerate(self.mesher.stations):
            af_name = st.get("airfoil_name")
            if af_name and af_name.lower() != "none":
                name_lower = af_name.lower()
                chord_val = float(st["chord"])
                
                # Compute Reynolds number: Re = Vinf * c / nu
                re = (v_inf * chord_val) / 1.46e-5
                re_target = min([50000, 100000, 200000, 500000, 1000000], key=lambda x: abs(x - re))
                
                print(f"  Station {idx} ({af_name}): Chord = {chord_val:.2f}m | Calculated Re = {re:,.0f} -> Matched Re = {re_target:,}")
                
                key = (name_lower, re_target)
                if key not in self.mesher.airfoils:
                    af = AeroFoil(af_name, self.mesher.data_folder)
                    af.downloadPerformance(re_target)
                    self.mesher.airfoils[key] = af
                self.station_airfoils.append(self.mesher.airfoils[key])
            else:
                self.station_airfoils.append(None)

        # Initialize decambering angles
        self.decamber_angles = np.zeros(self.nPoints)
        error = tol + 1.0
        iteration = 0
        
        # Cache spanwise coordinate y_i and local chord c_i for each section
        y_sec = []
        chord_sec = []
        twist_sec = []
        
        for i in range(self.n_span_sections):
            indices = list(range(i * self.n_chord_sections, (i + 1) * self.n_chord_sections))
            y_i = np.mean([self.panels[idx].Ycontrol for idx in indices])
            y_sec.append(y_i)
            chord_sec.append(self.mesher.get_chord(y_i))
            
            # Interpolate twist at this section
            _, _, _, twist_rad_i, _ = self.mesher._interpolate_station_properties(y_i, 0.5)
            twist_sec.append(twist_rad_i)
            
        y_sec = np.array(y_sec)
        chord_sec = np.array(chord_sec)
        twist_sec = np.array(twist_sec)

        while (error > tol) and (iteration < max_iter):
            iteration += 1
            
            # 1. Update decambering angles in VLM solver and solve VLM
            self.vlm.decamber_angles = np.copy(self.decamber_angles)
            self.vlm.resetPanels(self.panels) # Rebuilds RHS with new decambering angles
            self.vlm.solveSystem()

            # 2. Get spanwise circulation distribution
            gamma_sec = []
            for i in range(self.n_span_sections):
                indices = list(range(i * self.n_chord_sections, (i + 1) * self.n_chord_sections))
                gamma_i = sum(self.panels[idx].gamma for idx in indices)
                gamma_sec.append(gamma_i)
            gamma_sec = np.array(gamma_sec)

            # 3. Compute local induced angle of attack via VortexLineTool
            vortex_tool = VortexLineTool(gamma_sec, y_sec, v_inf)
            alfa_i = vortex_tool.getInducedAngle() # in radians

            # 4. Calculate effective angle of attack and target Cl
            alfa_eff = alfa_wing + twist_sec - alfa_i
            cl_target = np.array([self._get_airfoil_lift(y_sec[idx], alfa_eff[idx]) for idx in range(self.n_span_sections)])
            cl_vlm = 2.0 * gamma_sec / (v_inf * chord_sec)

            # 5. Compute decambering angle correction (using 2*pi linear slope)
            # delta_decamber_sec = alfa_eff - cl_target / (2*pi)
            delta_decamber_target = (cl_vlm - cl_target) / (2.0 * pi)

            # Distribute section decambering correction back to panels in the section
            decamber_new = np.zeros(self.nPoints)
            for i in range(self.n_span_sections):
                indices = list(range(i * self.n_chord_sections, (i + 1) * self.n_chord_sections))
                for idx in indices:
                    decamber_new[idx] = self.decamber_angles[idx] + damping * delta_decamber_target[i]

            error = float(np.linalg.norm(decamber_new - self.decamber_angles))
            self.decamber_angles = decamber_new

        # Final solve and coefficients calculation
        self.vlm.decamber_angles = np.copy(self.decamber_angles)
        self.vlm.resetPanels(self.panels)
        self.vlm.solveSystem()

        # Compute final local effective angle of attack for profile drag lookup
        gamma_sec = []
        for i in range(self.n_span_sections):
            indices = list(range(i * self.n_chord_sections, (i + 1) * self.n_chord_sections))
            gamma_i = sum(self.panels[idx].gamma for idx in indices)
            gamma_sec.append(gamma_i)
        gamma_sec = np.array(gamma_sec)

        vortex_tool = VortexLineTool(gamma_sec, y_sec, v_inf)
        alfa_i = vortex_tool.getInducedAngle()
        alfa_eff = alfa_wing + twist_sec - alfa_i

        # Compute profile drag for each section from polars
        self.cd_profile_sec = np.array([self._get_airfoil_drag(y_sec[idx], alfa_eff[idx]) for idx in range(self.n_span_sections)])
        
        # Integrate profile drag: C_Dp = sum(cd_profile * c * dy) / S
        dy_sec = []
        for i in range(self.n_span_sections):
            idx = i * self.n_chord_sections
            dy_sec.append(self.panels[idx].DY)
        dy_sec = np.array(dy_sec)
        
        # Compute coefficients using VLM (with local force calculations)
        CDi, CS, CL = self.vlm.getAerodynamicCoefficients()
        
        # Normalize VLM coefficients by Vinf^2 (since getAerodynamicCoefficients does not normalize by Vinf^2)
        CDi = CDi / (v_inf**2)
        CL = CL / (v_inf**2)
        
        S = self.vlm.surfaceArea
        CDp = float(sum(self.cd_profile_sec * chord_sec * dy_sec) / S)
        CD = CDi + CDp
        
        AR = (self.mesher.span**2) / S
        efficiency = (CL**2) / (pi * AR * CDi) if CDi > 1e-12 else 0.0

        return self.vlm.controlPoints, self.vlm.leadingEdgeGamma, iteration, error, CL, CDi, CDp, CD, efficiency
