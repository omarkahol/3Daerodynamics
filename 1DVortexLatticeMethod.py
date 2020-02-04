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

#WINGA DATA
wingSpan = 20
chord = 2
baseChord = 1
nPoints = 51
wingSpanLine = np.linspace(-0.5*wingSpan, 0.5*wingSpan,nPoints)

#DEFINE WING GEOMETRY
trailingEdgeSlope= 0.5
leadingEdgeSlope = 3
trailingEdgeBound = lambda y: trailingEdgeSlope*(baseChord/wingSpan)*abs(y) + chord
leadingEdgeBound = lambda y: leadingEdgeSlope*(baseChord/wingSpan)*abs(y)

panels=[]
for i in range(1,nPoints):
  y1 = wingSpanLine[i-1]
  y2 = wingSpanLine[i]
  y3=y1
  y4=y2

  x1 = leadingEdgeBound(y1)
  x3 = trailingEdgeBound(y3)

  x2 = leadingEdgeBound(y2)
  x4 = trailingEdgeBound(y4)
 
  panel = VortexPanel((x1,y1,0),(x2,y2,0),(x3,y3,0),(x4,y4,0))
  panels.append(panel)


#DRAW GEOMETRY
fig = plt.figure()
ax = fig.add_subplot(111,ylim=(-0.5*wingSpan,0.5*wingSpan), xlim=(0,2*chord))
for panel in panels:
  panel.drawPanel2D(ax)
ax.set_ylabel('y [m]')
ax.set_xlabel('x [m]')
ax.set_title('WING GEOMETRY AND DISCRETIZATION')
plt.show()

#WING DATA
alfa = 5
Vinf = 1

#INITIALIZE SOLVER AND SOLVE SYSTEM
solver = vortexLatticeSolver(panels)
solver.setDynamicData(alfa, Vinf)
solver.solveSystem()

#PLOT RESULTS
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('wingSpan [m]')
ax.set_ylabel('gamma [m^2/s]')
ax.set_title('VORTEX STRENGTH')
ax.plot(solver.controlPoints, solver.leadingEdgeGamma, 'k-',lw=2)
plt.show()

#COMPUTE THE TAPER RATION
ctChord = trailingEdgeBound(0)-leadingEdgeBound(0)
crChord = trailingEdgeBound(0.5*wingSpan)-leadingEdgeBound(0.5*wingSpan)
taperRatio = ctChord/crChord

#COMPUTE THE SWEEP ANGLE AT THE HALF CHORD LINE
xEnd = leadingEdgeBound(0.5*wingSpan) + 0.5*(trailingEdgeBound(0.5*wingSpan)-leadingEdgeBound(0.5*wingSpan))
xStart = leadingEdgeBound(0) + 0.5*(trailingEdgeBound(0)-leadingEdgeBound(0))
lambdaAngle = atan((xEnd-xStart)/(0.5*wingSpan))

#COMPUTE AERODINAMIC DATA
CDi, CS, CL = solver.getAerodynamicCoefficients()
AR = (wingSpan**2)/solver.surfaceArea
print('''
------------------------
AERODYNAMIC DATA
  surface --> {0} m^2
  AR --> {1}
  CL --> {2}
  CDi --> {3}
------------------------'''.format(round(solver.surfaceArea,3),round(AR,4),round(CL,4),round(CDi,4)))

#VALIDATE AERODYNAMIC DATA WITH KUCHEMANN'S FORMULA
a0 = 2*pi
dCLda = a0*cos(lambdaAngle)/(sqrt(1+(a0*cos(lambdaAngle)/(pi*AR))) + (a0*cos(lambdaAngle)/(pi*AR)))
print('''
------------------------
DATA VALIDATION
  taper ratio --> {0}
  sweep angle --> {1} [°deg]
  dCLda --> {2}
  CL --> {3}
------------------------'''.format(round(taperRatio,3),round(degrees(lambdaAngle),3), round(dCLda,3), round(dCLda*radians(alfa),5)))