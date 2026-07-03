import numpy as np
from math import sin, cos, radians, pi
from typing import List
from .panel import VortexPanel

class VortexLatticeSolver:
    """Non-viscous solver for the Vortex Lattice Method (VLM)."""

    def __init__(self, panels: List[VortexPanel]):
        self.panels = panels
        self.surfaceArea = 0.0
        self.alfa = 0.0
        self.Vinf = 0.0
        self.nPoints = len(self.panels)
        self.velocityMatrix = np.ones((self.nPoints, self.nPoints))
        self.RHS = np.ones(self.nPoints)
        
        self.iVector = np.array([1.0, 0.0, 0.0])
        self.jVector = np.array([0.0, 1.0, 0.0])
        self.kVector = np.array([0.0, 0.0, 1.0])
        
        self.controlPoints = []
        self.leadingEdgeGamma = []

    def setDynamicData(self, alfa: float, Vinf: float, isDegrees: bool = True) -> bool:
        """Sets the angle of attack and freestream velocity, then builds the system of equations."""
        self.alfa = radians(alfa) if isDegrees else alfa
        self.Vinf = abs(Vinf)
        self.rotationMatrix = np.array([
            [cos(self.alfa), 0.0, sin(self.alfa)],
            [0.0, 1.0, 0.0],
            [-sin(self.alfa), 0.0, cos(self.alfa)]
        ])
        self._buildSystem()
        return True

    def resetPanels(self, panels: List[VortexPanel]) -> bool:
        """Resets the panel distribution and rebuilds the system matrix."""
        self.panels = panels
        self.nPoints = len(self.panels)
        self.velocityMatrix = np.ones((self.nPoints, self.nPoints))
        self.RHS = np.ones(self.nPoints)
        self._buildSystem()
        return True

    def _buildSystem(self) -> None:
        self.controlPoints = []
        self.leadingEdgeGamma = []
        
        # Pre-extract panel geometric data for vectorization / optimization
        xm = np.array([p.Xcontrol for p in self.panels])
        ym = np.array([p.Ycontrol for p in self.panels])
        zm = np.array([p.Zcontrol for p in self.panels])

        x1n = np.array([p.Xapply1 for p in self.panels])
        x2n = np.array([p.Xapply2 for p in self.panels])
        y1n = np.array([p.Yapply1 for p in self.panels])
        y2n = np.array([p.Yapply2 for p in self.panels])
        z1n = np.array([p.Zapply1 for p in self.panels])
        z2n = np.array([p.Zapply2 for p in self.panels])

        phi_n = np.array([p.phi for p in self.panels])
        delta_n = np.array([p.delta for p in self.panels])

        # Precompute normal vectors
        nVectors = np.stack([
            -np.sin(delta_n) * np.cos(phi_n),
            -np.cos(delta_n) * np.sin(phi_n),
            np.cos(phi_n) * np.cos(delta_n)
        ], axis=1) # shape: (nPoints, 3)

        # Build RHS: self.RHS[m] = -self.Vinf * sin(alfa - (delta_m + decamber_m)) * cos(phi_m)
        # Note: delta and phi here are of panel m
        phi_m = np.array([p.phi for p in self.panels])
        delta_m = np.array([p.delta for p in self.panels])
        if not hasattr(self, 'decamber_angles') or len(self.decamber_angles) != self.nPoints:
            self.decamber_angles = np.zeros(self.nPoints)
        self.RHS = -self.Vinf * np.sin(self.alfa - (delta_m + self.decamber_angles)) * np.cos(phi_m)

        # We can optimize the matrix building using numpy broadcasting
        # Let's compute influence coefficients
        # For each control point m (0 to nPoints-1) and panel n (0 to nPoints-1)
        for m in range(self.nPoints):
            # Control point coordinates
            x_m, y_m, z_m = xm[m], ym[m], zm[m]

            # Vector from apply1 to apply2 for all n panels
            r0 = np.stack([x2n - x1n, y2n - y1n, z2n - z1n], axis=1) # shape: (nPoints, 3)
            # Vector from apply1 to control point m for all n panels
            r1 = np.stack([x_m - x1n, y_m - y1n, z_m - z1n], axis=1) # shape: (nPoints, 3)
            # Vector from apply2 to control point m for all n panels
            r2 = np.stack([x_m - x2n, y_m - y2n, z_m - z2n], axis=1) # shape: (nPoints, 3)

            # Norms of r1 and r2
            norm_r1 = np.linalg.norm(r1, axis=1)
            norm_r2 = np.linalg.norm(r2, axis=1)

            # Cross product of r1 and r2
            cross_r1_r2 = np.cross(r1, r2)
            norm_cross_sq = np.sum(cross_r1_r2**2, axis=1)
            # Avoid division by zero
            norm_cross_sq[norm_cross_sq < 1e-12] = 1e-12

            # F10 and F20 for segment AB
            F10 = cross_r1_r2 / norm_cross_sq[:, np.newaxis]
            
            # Prevent division by zero in r1 / norm_r1
            norm_r1_safe = np.where(norm_r1 < 1e-12, 1e-12, norm_r1)
            norm_r2_safe = np.where(norm_r2 < 1e-12, 1e-12, norm_r2)
            r1_unit = r1 / norm_r1_safe[:, np.newaxis]
            r2_unit = r2 / norm_r2_safe[:, np.newaxis]
            F20 = np.sum(r0 * (r1_unit - r2_unit), axis=1)

            VAB = (0.25 / pi) * F10 * F20[:, np.newaxis]

            # Segment A-infinity
            denom_A = (z_m - z1n)**2 + (y1n - y_m)**2
            denom_A[denom_A < 1e-12] = 1e-12
            
            # F11 shape: (nPoints, 3)
            F11 = (0.25 / pi) * np.stack([
                np.zeros_like(z1n),
                (z_m - z1n) / denom_A,
                (y1n - y_m) / denom_A
            ], axis=1)
            F21 = 1.0 + (x_m - x1n) / norm_r1_safe
            VAI = F11 * F21[:, np.newaxis]

            # Segment B-infinity
            denom_B = (z_m - z2n)**2 + (y2n - y_m)**2
            denom_B[denom_B < 1e-12] = 1e-12
            F12 = (0.25 / pi) * np.stack([
                np.zeros_like(z2n),
                (z_m - z2n) / denom_B,
                (y2n - y_m) / denom_B
            ], axis=1)
            F22 = 1.0 + (x_m - x2n) / norm_r2_safe
            VBI = -F12 * F22[:, np.newaxis]

            VTOT = VAB + VAI + VBI # shape: (nPoints, 3)

            # Dot product with normal vectors of all n panels
            self.velocityMatrix[m, :] = np.sum(VTOT * nVectors, axis=1)

    def solveSystem(self) -> bool:
        """Solves the system equations for gamma strength on each panel."""
        self.solution = np.linalg.solve(self.velocityMatrix, self.RHS)
        self.surfaceArea = 0.0
        
        gammaDict = {}
        for i, gamma in enumerate(self.solution):
            self.panels[i].setGamma(gamma)
            self.surfaceArea += self.panels[i].surfaceArea
            
            # Use key based on spanwise control point coordinate
            key = str(round(self.panels[i].Ycontrol, 5))
            if key in gammaDict:
                gammaDict[key].append(self.panels[i].gamma)
            else:
                gammaDict[key] = [self.panels[i].gamma]
                
        for key in sorted(gammaDict.keys(), key=float):
            self.leadingEdgeGamma.append(sum(gammaDict[key]))
            self.controlPoints.append(float(key))
            
        self.leadingEdgeGamma = np.array(self.leadingEdgeGamma)
        self.controlPoints = np.array(self.controlPoints)
        return True

    def _computeVelocityField(self) -> None:
        # Pre-extract panel coordinates (using bound vortex midpoint instead of control point for force calculation)
        x_mid = np.array([0.5 * (p.Xapply1 + p.Xapply2) for p in self.panels])
        y_mid = np.array([0.5 * (p.Yapply1 + p.Yapply2) for p in self.panels])
        z_mid = np.array([0.5 * (p.Zapply1 + p.Zapply2) for p in self.panels])

        x1n = np.array([p.Xapply1 for p in self.panels])
        x2n = np.array([p.Xapply2 for p in self.panels])
        y1n = np.array([p.Yapply1 for p in self.panels])
        y2n = np.array([p.Yapply2 for p in self.panels])
        z1n = np.array([p.Zapply1 for p in self.panels])
        z2n = np.array([p.Zapply2 for p in self.panels])
        gammas = np.array([p.gamma for p in self.panels])

        for m in range(self.nPoints):
            # Start with Freestream velocity
            VTOT_m = self.Vinf * np.array([cos(self.alfa), 0.0, sin(self.alfa)])
            
            x_m, y_m, z_m = x_mid[m], y_mid[m], z_mid[m]

            r0 = np.stack([x2n - x1n, y2n - y1n, z2n - z1n], axis=1)
            r1 = np.stack([x_m - x1n, y_m - y1n, z_m - z1n], axis=1)
            r2 = np.stack([x_m - x2n, y_m - y2n, z_m - z2n], axis=1)

            norm_r1 = np.linalg.norm(r1, axis=1)
            norm_r2 = np.linalg.norm(r2, axis=1)

            cross_r1_r2 = np.cross(r1, r2)
            norm_cross_sq = np.sum(cross_r1_r2**2, axis=1)
            norm_cross_sq[norm_cross_sq < 1e-12] = 1e-12

            F10 = cross_r1_r2 / norm_cross_sq[:, np.newaxis]
            norm_r1_safe = np.where(norm_r1 < 1e-12, 1e-12, norm_r1)
            norm_r2_safe = np.where(norm_r2 < 1e-12, 1e-12, norm_r2)
            F20 = np.sum(r0 * (r1 / norm_r1_safe[:, np.newaxis] - r2 / norm_r2_safe[:, np.newaxis]), axis=1)

            VAB = (0.25 / pi) * F10 * F20[:, np.newaxis]
            # Exclude self-induction of bound segment
            VAB[m] = 0.0

            denom_A = (z_m - z1n)**2 + (y1n - y_m)**2
            denom_A[denom_A < 1e-12] = 1e-12
            F11 = (0.25 / pi) * np.stack([
                np.zeros_like(z1n),
                (z_m - z1n) / denom_A,
                (y1n - y_m) / denom_A
            ], axis=1)
            F21 = 1.0 + (x_m - x1n) / norm_r1_safe
            VAI = F11 * F21[:, np.newaxis]

            denom_B = (z_m - z2n)**2 + (y2n - y_m)**2
            denom_B[denom_B < 1e-12] = 1e-12
            F12 = (0.25 / pi) * np.stack([
                np.zeros_like(z2n),
                (z_m - z2n) / denom_B,
                (y2n - y_m) / denom_B
            ], axis=1)
            F22 = 1.0 + (x_m - x2n) / norm_r2_safe
            VBI = -F12 * F22[:, np.newaxis]

            V_induced = (VAB + VAI + VBI) * gammas[:, np.newaxis] # shape: (nPoints, 3)
            VTOT_m += np.sum(V_induced, axis=0)

            self.panels[m].setTotalInducedVelocity(VTOT_m)

    def getAerodynamicCoefficients(self) -> np.ndarray:
        """Computes and returns aerodynamic force coefficients (CDi, CS, CL)."""
        self._computeVelocityField()
        force = np.zeros(3)
        for panel in self.panels:
            force += panel.getForce()
        return (2.0 / self.surfaceArea) * self.rotationMatrix.dot(force)
