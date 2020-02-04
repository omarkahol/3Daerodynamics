__author__ = "Omar Kahol"
__copyright__ = "Copyright (C) 2020 Omar Kahol"
__license__ = "Public Domain"
__version__ = "1.0"

from math import *
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import quad, simps
from vortexPanel import VortexPanel
from vortexLatticeSolver import vortexLatticeSolver
from openAirfoilCoordinates import AeroFoil
from gammaTool import vortexLineTool

#READ PROFILE FILE
airFoilName = 'n2h15'
airFoil = AeroFoil(airFoilName,'airfoilDATA')
airFoil.downloadPerformance()
getCL = airFoil.fitCL(5)

#PLOT CL ALFA CURVE
fig = plt.figure()
ax = fig.add_subplot(111)
maxAngle = 15
xaxis = np.linspace(-maxAngle,maxAngle,100)
ax.plot(airFoil.alfa,airFoil.CL,'k-',lw=2, label='CL-alfa curve')
ax.plot(xaxis,[getCL(alfa) for alfa in np.radians(xaxis)],'r--',lw=1, label='fitted data')
ax.set_title('CL-ALPHA CURVE')
ax.set_xlabel('alpha (deg)')
ax.set_ylabel('CL')
plt.show()

#WING DATA 
wingSpan = 20
maxChord = 2
wingSlope=1.6
chord = lambda y: maxChord - wingSlope*(maxChord/wingSpan)*abs(y)
alfaWing = radians(5)
Vinf = 1

#PANELIZE THE WING COMPUTE THE SURFACE AREA AND DRAW RESULTS
nPanels = 50
panels = []
wingSpanLine = np.linspace(-0.5*wingSpan,0.5*wingSpan, nPanels)
surfaceArea = 0
for i in range(1,nPanels):
  Y1 = wingSpanLine[i-1]
  Y2 = wingSpanLine[i]
  Y3, Y4 = Y1, Y2
  X1 = 0.5*chord(Y1)
  X2 = 0.5*chord(Y2)
  X3, X4 = -X1, -X2
  panel = VortexPanel([X1,Y1,0],[X2,Y2,0],[X3,Y3,0],[X4,Y4,0])
  surfaceArea += panel.surfaceArea
  panels.append(panel)

fig = plt.figure()
ax = fig.add_subplot(111,ylim=(-0.5*wingSpan,0.5*wingSpan), xlim=(-0.5*maxChord,0.5*maxChord))
for panel in panels:
  panel.drawPanel2D(ax)
ax.set_ylabel('y [m]')
ax.set_xlabel('x [m]')
ax.set_title('WING DESIGN')
plt.show()


#SET ITERATION DATA
dampingFactor = 0.02
nPoints = 100
wingSpanLine = np.linspace(-0.5*wingSpan,0.5*wingSpan,nPoints)
chords = np.array([chord(Y0) for Y0 in wingSpanLine])
gammaOld = [chord(Y0) for Y0 in wingSpanLine]
tol = 1e-2
error = tol + 1
itmax = 1000
iteration = 0

#START ITERATIONS
while (error > tol) and (iteration < itmax):
  iteration = iteration + 1

  #START THE VORTEX TOOL CLASS
  vortexTool = vortexLineTool(gammaOld,wingSpanLine,Vinf)
  alfaI=vortexTool.getInducedAngle()

  #COMPUTE NEW CIRCULATION
  alfaEff = alfaWing - alfaI
  clWing = np.array([getCL(alfa) for alfa in alfaEff])
  gammaNew = 0.5*chords*clWing
  gammaNew[0]=0
  gammaNew[-1] = 0
  error = np.linalg.norm(gammaNew-gammaOld)
  gammaOld = gammaOld + dampingFactor*(gammaNew-gammaOld)

print('''----------------------------
ITERATION RESULTS
  iterations --> {0}
  final error --> {1}
----------------------------'''.format(iteration,round(error,5)))

#COMPUTE AERODINAMIC DATA
vortexTool = vortexLineTool(gammaOld,wingSpanLine,Vinf)
vortexTool.getInducedAngle()

AR = (wingSpan**2)/surfaceArea
CL = (2/surfaceArea) * vortexTool.getCLintegral()
CDi = (2/surfaceArea) * vortexTool.getCDintegral()
efficiency = (CL**2)/(pi*AR*CDi)

print('''
------------------------
AERODYNAMIC DATA
  surface --> {0} [m^2]
  alfa = {1} [°deg]
  Vinf = {2} [m/s]
  AR --> {3}
  CL --> {4}
  CDi --> {5}
  efficiency = {6}
------------------------'''.format(round(surfaceArea,3),round(degrees(alfaWing),3),round(Vinf),round(AR,4),round(CL,4),round(CDi,4), round(efficiency,4)))

#PLOT CIRCULATION
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wingSpanLine,gammaNew,'b-',lw=2)
ax.set_title('CIRCULATION')
ax.set_xlabel('wingSpan')
ax.set_ylabel('gamma')
plt.show()