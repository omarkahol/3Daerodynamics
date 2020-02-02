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
import matplotlib.animation as animation

#CHOOSE PROFILE
airFoilName = 'n0012'
airFoil = AeroFoil(airFoilName,'airfoilDATA')
getCamber = airFoil.fitCamberLine(4)

#GEOMETRY DELTA WING AND DISCRETIZATION
wingSpan = 10
maxChord = 10
wingBorder = lambda y: 2*(maxChord/wingSpan)*abs(y)
nPointsY = 5
nPointsX = 6
wingSpanLine = np.linspace(-0.5*wingSpan, 0.5*wingSpan,nPointsY)
def getPanels(theta):
  trailingEdgeBound = lambda y: leadingEdgeBound(y)+chord
  leadingEdgeBound = lambda y: 2*(baseChord/wingSpan)*abs(y)
  panels = []
  for i in range(2,nPointsY-1):
    Y10 = wingSpanLine[i-1]
    Y20 = wingSpanLine[i]

    Y30, Y40 = Y10, Y20

    Y1 = Y10*cos(theta)
    Y2 = Y20*cos(theta)
    Y3 = Y30*cos(theta)
    Y4 = Y40*cos(theta)

    lim1 = wingBorder(Y10)
    lim2 = wingBorder(Y20)

    xArray1 = np.linspace(lim1, maxChord, nPointsX)
    xArray2 = np.linspace(lim2, maxChord, nPointsX)

    for j in range(1, nPointsX):
      X1 = xArray1[j-1]
      X2 = xArray2[j-1]
      X3 = xArray1[j]
      X4 = xArray2[j]
      Z1 = getCamber((maxChord-X1)/(maxChord-lim1)) + abs(Y10)*sin(theta)
      Z2 = getCamber((maxChord-X2)/(maxChord-lim2)) + abs(Y20)*sin(theta)
      Z3 = getCamber((maxChord-X3)/(maxChord-lim1)) + abs(Y30)*sin(theta)
      Z4 = getCamber((maxChord-X4)/(maxChord-lim2)) + abs(Y40)*sin(theta)
      panels.append(VortexPanel((X1,Y1,Z1),(X2,Y2,Z2),(X3,Y3,Z3),(X4,Y4,Z4)))
  return panels

panels = getPanels(radians(20))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d',zlim=(-5,5), ylim=(-5,5), xlim=(0,10))

lineCollection = []
for panel in panels:
  line, = ax.plot([],[],[],'k-',lw=2)
  line1, = ax.plot([],[],[],'k-',lw=2)
  line2, = ax.plot([],[],[],'k-',lw=2)
  line3, = ax.plot([],[],[],'k-',lw=2)
  lineCollection.append([line,line1,line2,line3])

currentTime = 0
dt = 0.1
timeArray = []
CLarray = []
thetaArray=[]
Vinf = 1
alfa = 17
def animate(i):
  global dt, currentTime
  theta = (pi/4)*sin(currentTime)
  panels = getPanels(theta)
  timeArray.append(currentTime)
  thetaArray.append(theta)

  #INITIALIZE SOLVER
  solver = vortexLatticeSolver(panels)
  solver.setDynamicData(alfa,Vinf)
  solver.solveSystem()
  CL = (2/(Vinf*solver.surfaceArea))*simps(solver.leadingEdgeGamma, solver.controlPoints)

  CLarray.append(CL)
  for panel, lines in zip(panels, lineCollection):
    lines[0].set_data([panel.X1,panel.X2],[panel.Y1,panel.Y2])
    lines[1].set_data([panel.X2,panel.X4],[panel.Y2,panel.Y4])
    lines[2].set_data([panel.X4,panel.X3],[panel.Y4,panel.Y3])
    lines[3].set_data([panel.X3,panel.X1],[panel.Y3,panel.Y1])

    lines[0].set_3d_properties([panel.Z1,panel.Z2])
    lines[1].set_3d_properties([panel.Z2,panel.Z4])
    lines[2].set_3d_properties([panel.Z4,panel.Z3])
    lines[3].set_3d_properties([panel.Z3,panel.Z1])
  currentTime += dt
  return None

ani = animation.FuncAnimation(fig, animate, interval=1, blit=False)
plt.show(ani)
plt.plot(timeArray,CLarray)
plt.show()

plt.plot(thetaArray,CLarray)
plt.show()
