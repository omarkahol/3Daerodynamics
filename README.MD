# 3D AERODYNAMICS
3DAerodynamics is a collection of classes and scripts that can calculate the performances of 3D lifting surfaces. Three main approaches are succesfully implemented.

1) Prandtl's lifting line theory
2) Non linear lifting line theory
3) Vortex lattice method

Inside you can find:

## AeroFoil in openAirFoilCoordinates.py
a class that downloads and reads the coordinates and the performances of an airfoil. Every file is downloaded from the [University of Illinois Airfoil Database](http://airfoiltools.com/airfoil/). The parsed name should be the one indicated in their website.
```python
import AeroFoil

name = 'n0012' #OPEN NACA0012 AIRFOIL
folder = 'DATABASE' #LOCATION OF THE FILES
aer = AeroFoil(name,folder)
camberLambda = aer.fitCamberLine(4) #FIT THE MEAN CAMBER LINE, DEGREE = 4
clLambda = aer.fitCL(3) #FIT CL CURVE, DEGREE = 3
```
## vortexPanel
a class that stores the main vertices of a panel. It also computes relevat points such as the control point and the application points of the horseshoe vortex, stores the vortex strength and computes the dihedral and dz/dx angles.
The order is 

P1----P2         |------Y direction
|      |         |
|      |         |
P3----P4         X direction

The control point is located at about 75% of the mean Y distance and the points that generate the vortex at about 25% of the distance. This choice can cause numerical errors since 25 is a divisor of 75 so one panel can easily have the control coordinate equal to the generation coordinate of the other. This is why I have chosen two slightly different numbers. 

## vortexLatticeSolver
a class that implements a non viscous solver for the vortex lattice method. The input should be a panel array. The vortex generates an induced velocity that can be calculated using Biot-Savart's law. The unknown is the vortex strength and the boundary condition is that the normal component velocity at the control point of each panel is zero. This condition generates a linear system of equations that is solved with numpy's linear algebra tools.
```python
import vortexLatticeSolver

solver = vortexLatticeSolver(panels) #SET THE PANEL DISTRIBUTION
solver.setDynamicData(alfa, Vinf) #SET THE ANGLE OF ATTACK AND THE FREESTREAM VELOCITY 
#THE CLASS WILL AUTOMATICALLY BUILD THE VORTICITY MATRIX AND THE RHS

solver.solveSystem() #SOLVE THE SYSTEM OF EQUATIONS AND COMPUTE RELEVANT TERMS
```
The code solves the system of equations and sets each vortex strength to its panel. It also performs a chordwise sum to obatin the leading edge vortex strength used to calculate CL and CDi

## vortexLineTool in gammaTool.py
a class that stores a vortex strength distribution and calculates the induced angle of attack and the integrals necessary to calculate CL and CDi.
```python
import vortexLineTool
vortexTool = vortexLineTool(gamma,controlPoints) #THE INPUTS ARE THE VORTEX STRENGTH AND THE CONTROL POINTS
vortexTool.getInducedAngle()  #RETURNS A VECTOR THAT CONTAINS THE INDUCED ANGLE OF ATTACK AT EACH CONTROL POINT
CL = (2/surfaceArea) * vortexTool.getCLintegral() #COMPUTES quad(gamma)
CDi = (2/surfaceArea) * vortexTool.getCDintegral() #COMPUTES quad(gamma*alfai)
```

## other files
the other files are mainly examples. They include delta wings, non linear methods, flapping wings, ...
