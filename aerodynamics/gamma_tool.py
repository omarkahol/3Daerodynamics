import numpy as np
from math import pi
try:
    from scipy.integrate import simpson
except ImportError:
    from scipy.integrate import simps as simpson

class VortexLineTool:
    """Operations on vortex strength distribution along a line (e.g., spanwise gamma)."""

    def __init__(self, gamma: np.ndarray, controlPoints: np.ndarray, Vinf: float):
        self.gamma = np.array(gamma)
        self.Vinf = Vinf
        self.controlPoints = np.array(controlPoints)
        self.DY = abs(self.controlPoints[1] - self.controlPoints[0]) if len(self.controlPoints) > 1 else 0.0
        self.nPoints = len(self.gamma)
        self.inducedAngle = np.zeros(self.nPoints)
        self._compute_dgamma_dy()

    def _compute_dgamma_dy(self) -> None:
        self.DgammaDy = np.gradient(self.gamma, self.controlPoints)

    def getInducedAngle(self) -> np.ndarray:
        """Computes the induced angle of attack at each control point.
        
        Handles the integrand singularity by taking the average of adjacent values.
        """
        self.inducedAngle = []
        multiplier = 1 / (4 * pi * self.Vinf)
        
        for singularity_idx, Y0 in enumerate(self.controlPoints):
            # Compute denominators and replace the singularity point (zero) with 1 to avoid division by zero
            denominators = Y0 - self.controlPoints
            denominators[singularity_idx] = 1.0

            # Integrand function: dGamma/dy / (Y0 - y)
            function = self.DgammaDy / denominators

            # Replace the singular value with average of adjacent values
            if singularity_idx == 0:
                function[singularity_idx] = 0.0
            elif singularity_idx == self.nPoints - 1:
                function[singularity_idx] = 0.0
            else:
                function[singularity_idx] = 0.5 * (function[singularity_idx - 1] + function[singularity_idx + 1])

            # Perform integration using Simpson's rule
            integral = simpson(function, x=self.controlPoints)
            self.inducedAngle.append(multiplier * integral)

        self.inducedAngle = np.array(self.inducedAngle)
        return np.copy(self.inducedAngle)

    def getCLintegral(self) -> float:
        """Integrates gamma across the span."""
        return float(simpson(self.gamma, x=self.controlPoints))

    def getCDintegral(self) -> float:
        """Integrates (gamma * induced_angle) across the span."""
        return float(simpson(self.gamma * self.inducedAngle, x=self.controlPoints))
