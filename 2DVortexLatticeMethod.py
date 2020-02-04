__author__ = "Omar Kahol"
__copyright__ = "Copyright (C) 2020 Omar Kahol"
__license__ = "Public Domain"
__version__ = "1.0"

from math import *
import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import simps
from vortexPanel import VortexPanel
from vortexLatticeSolver import vortexLatticeSolver
from openAirfoilCoordinates import AeroFoil
from gammaTool import vortexLineTool

#GEOMETRY --> DELTA BORDER
wingSpan = 10
maxChord = 10
wingBorder = lambda y: 2*(maxChord/wingSpan)*abs(y)

#GEOMETRY --> DELTA NOTCH
maxNotch = 0.1*maxChord
notch = lambda y: maxNotch - 2*(maxNotch/wingSpan)*abs(y)

#DISCRETIZE
nPointsY = 21
nPointsX = 11
wingSpanLine = np.linspace(-0.5*wingSpan, 0.5*wingSpan,nPointsY)
panels=[]
for i in range(1,nPointsY):
  Y1 = wingSpanLine[i-1]
  Y2 = wingSpanLine[i]

  Y3, Y4 = Y1, Y2
  lim1 = wingBorder(Y1)
  lim2 = wingBorder(Y2)

  xArray1 = np.linspace(lim1, maxChord-notch(Y1), nPointsX)
  xArray2 = np.linspace(lim2, maxChord-notch(Y2), nPointsX)

  for j in range(1, nPointsX):
    X1 = xArray1[j-1]
    X2 = xArray2[j-1]
    X3 = xArray1[j]
    X4 = xArray2[j]
    panels.append(VortexPanel((X1,Y1,0),(X2,Y2,0),(X3,Y3,0),(X4,Y4,0)))

#DRAW GEOMETRY
fig = plt.figure()
ax = fig.add_subplot(111,ylim=(-0.5*wingSpan,0.5*wingSpan))
for panel in panels:
  panel.drawPanel2D(ax)
ax.set_ylabel('y [m]')
ax.set_xlabel('x [m]')
ax.set_title('WING GEOMETRY AND DISCRETIZATION')
plt.show()

#WING DATA
alfa = 10 #degrees
Vinf = 1 #m/s

#INITIALIZE SOLVER AND SOLVE SYSTEM
solver = vortexLatticeSolver(panels)
solver.setDynamicData(alfa,Vinf)
solver.solveSystem()

#PLOT RESULTS
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter([panel.Xcontrol for panel in panels],[panel.Ycontrol for panel in panels],solver.solution)
ax.plot([maxChord,maxChord-notch(0)],[-0.5*wingSpan, 0],[0,0], 'k-')
ax.plot([maxChord-notch(0), maxChord],[0, 0.5*maxChord],[0,0], 'k-')
ax.plot([0,maxChord],[0, 0.5*wingSpan],[0,0], 'k-')
ax.plot([0,maxChord],[0,-0.5*wingSpan],[0,0], 'k-')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('gamma [m^2/s]')
ax.set_title('VORTEX STRENGHT DISTRIBUTION')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(solver.controlPoints,solver.leadingEdgeGamma,'k-',lw=2)
ax.set_xlabel('wingSpan [m]')
ax.set_ylabel('gamma [m^2/s]')
ax.set_title('LEADING EDGE TOTAL CIRCULATION')
plt.show()

#COMPUTE AERODINAMIC DATA
CDi, CS, CL = solver.getAerodynamicCoefficients()
AR = (wingSpan**2)/solver.surfaceArea
efficiency=(CL**2)/(pi*AR*CDi)
print('''
------------------------
AERODYNAMIC DATA
  surface --> {0} m^2
  AR --> {1}
  CL --> {2}
  CDi --> {3}
  efficiency --> {4}
------------------------'''.format(round(solver.surfaceArea,3),round(AR,4),round(CL,4),round(CDi,4), round(efficiency,3)))