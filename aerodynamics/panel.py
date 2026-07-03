from math import atan2
import numpy as np
from typing import Tuple, List

class VortexPanel:
    """Represents a 3D panel on a lifting surface, including application and control points for VLM."""

    def __init__(self, p1: Tuple[float, float, float], p2: Tuple[float, float, float], 
                 p3: Tuple[float, float, float], p4: Tuple[float, float, float]):
        # Basic Unit Vectors
        self.iVector = np.array([1.0, 0.0, 0.0])
        self.jVector = np.array([0.0, 1.0, 0.0])
        self.kVector = np.array([0.0, 0.0, 1.0])

        # Vertices
        self.X1, self.X2, self.X3, self.X4 = p1[0], p2[0], p3[0], p4[0]
        self.Y1, self.Y2, self.Y3, self.Y4 = p1[1], p2[1], p3[1], p4[1]
        self.Z1, self.Z2, self.Z3, self.Z4 = p1[2], p2[2], p3[2], p4[2]

        # Geometry dimensions
        self.DY = abs(self.Y2 - self.Y1)
        self.DX1 = abs(self.X3 - self.X1)
        self.DX2 = abs(self.X4 - self.X2)
        self.DX = 0.5 * (self.DX1 + self.DX2)

        # Control and application point ratios
        self.tqChord = 0.77
        self.oqChord = 0.25

        # Control Point
        self.Xcontrol = 0.5 * (self.X1 + self.tqChord * (self.X3 - self.X1) + 
                               self.X2 + self.tqChord * (self.X4 - self.X2))
        self.Ycontrol = 0.5 * (self.Y1 + self.Y2)
        self.Zcontrol = 0.5 * (self.Z1 + self.tqChord * (self.Z3 - self.Z1) + 
                               self.Z2 + self.tqChord * (self.Z4 - self.Z2))

        # Application Points
        self.Yapply1 = self.Y1
        self.Yapply2 = self.Y2
        
        self.Xapply1 = self.X1 + self.oqChord * (self.X3 - self.X1)
        self.Xapply2 = self.X2 + self.oqChord * (self.X4 - self.X2)

        self.Zapply1 = self.Z1 + self.oqChord * (self.Z3 - self.Z1)
        self.Zapply2 = self.Z2 + self.oqChord * (self.Z4 - self.Z2)

        # Geometric Angles
        self.phi = 0.5 * atan2(self.Z2 - self.Z1, self.Y2 - self.Y1) + 0.5 * atan2(self.Z4 - self.Z3, self.Y4 - self.Y3)

        tol = 1e-6
        if abs(self.X1 - self.X3) < tol:
            self.delta = 0.5 * atan2(self.Z4 - self.Z2, self.X4 - self.X2)
        elif abs(self.X4 - self.X2) < tol:
            self.delta = 0.5 * atan2(self.Z3 - self.Z1, self.X3 - self.X1)
        else:
            self.delta = 0.5 * atan2(self.Z3 - self.Z1, self.X3 - self.X1) + 0.5 * atan2(self.Z4 - self.Z2, self.X4 - self.X2)

        self.surfaceArea = self.DY * self.DX
        self.gamma = 0.0
        self.inducedVelocity = np.zeros(3)

    def drawPanel3D(self, ax) -> List:
        """Draws the panel in a 3D matplotlib plot."""
        ax.plot([self.X1, self.X2], [self.Y1, self.Y2], [self.Z1, self.Z2], 'k-', lw=1.5)
        ax.plot([self.X2, self.X4], [self.Y2, self.Y4], [self.Z2, self.Z4], 'k-', lw=1.5)
        ax.plot([self.X4, self.X3], [self.Y4, self.Y3], [self.Z4, self.Z3], 'k-', lw=1.5)
        ax.plot([self.X3, self.X1], [self.Y3, self.Y1], [self.Z3, self.Z1], 'k-', lw=1.5)

        ap1, = ax.plot([self.Xapply1], [self.Yapply1], [self.Zapply1], 'b*', lw=1, label='application point 1')
        ap2, = ax.plot([self.Xapply2], [self.Yapply2], [self.Zapply2], 'b+', lw=1, label='application point 2')
        cp, = ax.plot([self.Xcontrol], [self.Ycontrol], [self.Zcontrol], 'ro', lw=1, label='control point')
        return [ap1, ap2, cp]

    def drawPanel2D(self, ax) -> bool:
        """Draws the panel in a 2D matplotlib plot."""
        ax.plot([self.X1, self.X2], [self.Y1, self.Y2], 'k-', lw=1.5)
        ax.plot([self.X2, self.X4], [self.Y2, self.Y4], 'k-', lw=1.5)
        ax.plot([self.X4, self.X3], [self.Y4, self.Y3], 'k-', lw=1.5)
        ax.plot([self.X3, self.X1], [self.Y3, self.Y1], 'k-', lw=1.5)

        ax.plot([self.Xapply1], [self.Yapply1], 'b*', lw=1)
        ax.plot([self.Xapply2], [self.Yapply2], 'b+', lw=1)
        ax.plot([self.Xcontrol], [self.Ycontrol], 'ro', lw=1)
        return True

    def setGamma(self, gamma: float) -> None:
        """Sets the vortex strength of the panel."""
        self.gamma = gamma

    def setTotalInducedVelocity(self, vel: np.ndarray) -> None:
        """Sets the local induced velocity on the panel."""
        self.inducedVelocity = np.array(vel)

    def getForce(self) -> np.ndarray:
        """Computes the aerodynamic force on the panel using Kutta-Joukowski theorem."""
        self.R0 = ((self.Xapply2 - self.Xapply1) * self.iVector + 
                   (self.Yapply2 - self.Yapply1) * self.jVector + 
                   (self.Zapply2 - self.Zapply1) * self.kVector)
        return np.cross(self.inducedVelocity, self.R0) * self.gamma
