__author__ = "Omar Kahol"
__copyright__ = "Copyright (C) 2020 Omar Kahol"
__license__ = "Public Domain"
__version__ = "1.0"

from math import *
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import simps
from vortexPanel import VortexPanel
from vortexLatticeSolver import vortexLatticeSolver
from openAirfoilCoordinates import AeroFoil
from gammaTool import vortexLineTool


#OPEN THE AIRFOIL PROFILE AND GET THE CAMBER LINE
airFoilName = 'ag18'
airFoil = AeroFoil(airFoilName,'airfoilDATA')
getCamber = airFoil.fitCamberLine(4)

#PLOT THE MEAN CAMBER LINE
fig=plt.figure()
ylim = 5*max(airFoil.yBase)
ax=fig.add_subplot(111,ylim=(-ylim,ylim))
ax.plot(airFoil.xBase,airFoil.yBase,'k-',label='profile')
xAxis = np.linspace(0,1,100)
yCamberData= [getCamber(x) for x in xAxis]
ax.plot(xAxis,yCamberData,'r--', label='fitted camber line')
ax.legend()
ax.set_title('AIRFOIL PROFILE')
plt.show()


#GEOMETRY --> DELTA WING
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
    if abs(maxChord-lim1)<1e-4:
      Z1, Z3 = 0, 0
    else:
      Z1 = getCamber((maxChord-X1)/(maxChord-lim1))
      Z3 = getCamber((maxChord-X3)/(maxChord-lim1))
    if abs(maxChord-lim2)<1e-4:
      Z2, Z4 = 0, 0
    else:
      Z2 = getCamber((maxChord-X2)/(maxChord-lim2))
      Z4 = getCamber((maxChord-X4)/(maxChord-lim2))
    panels.append(VortexPanel((X1,Y1,Z1),(X2,Y2,Z2),(X3,Y3,Z3),(X4,Y4,Z4)))

#DRAW WING
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d',aspect='equal',zlim=(0,0.2))
ax.set_title('WING GEOMETRY')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')
for panel in panels:
  ap1, ap2, cp= panel.drawPanel3D(ax)
ax.legend([ap1,ap2,cp], ['application point 1', 'application point 2', 'control point'])
plt.show()

#WING DATA
Vinf = 1
alfa = 10

#BUILD THE SOLVER
solver = vortexLatticeSolver(panels)
solver.setDynamicData(alfa, Vinf)

#SOLVE THE SYSTEM
solver.solveSystem()

#DISPLAY VORTEX STRENGTH DISTRIBUITION
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('gamma [m^2/s]')
ax.set_title('VORTEX STRENGHT DISTRIBUTION')
ax.scatter([panel.Xcontrol for panel in panels],[panel.Ycontrol for panel in panels],solver.solution)
ax.plot([maxChord,maxChord-notch(0)],[-0.5*wingSpan, 0],[0,0], 'k-')
ax.plot([maxChord-notch(0), maxChord],[0, 0.5*maxChord],[0,0], 'k-')
ax.plot([0,maxChord],[0, 0.5*wingSpan],[0,0], 'k-')
ax.plot([0,maxChord],[0,-0.5*wingSpan],[0,0], 'k-')
plt.show()

#DISPLAY LEADING EDGE VORTEX STRENGTH DISTRIBUTION
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(solver.controlPoints,solver.leadingEdgeGamma,'k-',lw=2)
ax.set_xlabel('wingSpan [m]')
ax.set_ylabel('gamma [m^2/s]')
ax.set_title('LEADING EDGE TOTAL CIRCULATION')
plt.show()

#COMPUTE AERODINAMIC DATA
AR = (wingSpan**2)/solver.surfaceArea
CDi, CS, CL = solver.getAerodynamicCoefficients()
efficiency = (CL**2)/(pi*AR*CDi)
print('''
------------------------
AERODYNAMIC DATA
  surface --> {0} m^2
  AR --> {1}
  CL --> {2}
  CDi --> {3}
  efficiency --> {4}
------------------------'''.format(round(solver.surfaceArea,3),round(AR,4),round(CL,4),round(CDi,4),round(efficiency,3)))