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

#GEOMETRY AND WING DISCRETIZATION
wingSpan = 10
chord = 2
baseChord = 1
nPoints = 21
wingSpanLine = np.linspace(-0.5*wingSpan, 0.5*wingSpan,nPoints)

trailingEdgeBound = lambda y: leadingEdgeBound(y)+chord
leadingEdgeBound = lambda y: 2*(baseChord/wingSpan)*abs(y)

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
ax = fig.add_subplot(111,ylim=(-0.5*wingSpan,0.5*wingSpan), xlim=(0,trailingEdgeBound(0.5*wingSpan)))
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


#COMPUTE AERODINAMIC DATA
vortexTool = vortexLineTool(solver.leadingEdgeGamma,solver.controlPoints)
vortexTool.getInducedAngle()
AR = (wingSpan**2)/solver.surfaceArea
CL = (2/solver.surfaceArea) * vortexTool.getCLintegral()
CDi = (2/solver.surfaceArea) * vortexTool.getCDintegral()

print('''
------------------------
AERODYNAMIC DATA
  surface --> {0} m^2
  AR --> {1}
  CL --> {2}
  CDi --> {3}
------------------------'''.format(round(solver.surfaceArea,3),round(AR,4),round(CL,4),round(CDi,4)))