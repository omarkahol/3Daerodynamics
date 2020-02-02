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
airFoilName = 'ag18'
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
wingSpan = 16
maxChord = 2
chord = lambda y: maxChord*sqrt(1-(2*y/wingSpan)**2)
alphaWing = radians(5)

#PLOT WING
fig = plt.figure()
lim = (-0.5*wingSpan,0.5*wingSpan)
ax = fig.add_subplot(111,aspect='equal',xlim=lim,ylim=lim)
nPoints = 100
t = np.linspace(0,2*pi,nPoints)
ax.plot(0.5*wingSpan*np.cos(t),0.5*maxChord*np.sin(t))
ax.set_title('WING GEOMETRY')
plt.show()

#SET ITERATION DATA
dampingFactor = 0.05
nPoints = 150
wingSpanLine = np.linspace(-0.5*wingSpan,0.5*wingSpan,nPoints)
chords = np.array([chord(Y0) for Y0 in wingSpanLine])
gammaOld = 0.5*np.ones(nPoints)
tol = 1e-4
error = tol + 1
itmax = 1000
iteration = 0

#START ITERATIONS
while (error > tol) and (iteration < itmax):
  iteration = iteration + 1

  #START THE VORTEX TOOL CLASS
  vortexTool = vortexLineTool(gammaOld,wingSpanLine)
  alfaI=vortexTool.getInducedAngle()

  #COMPUTE NEW CIRCULATION
  alfaEff = alphaWing - alfaI
  clWing = np.array([getCL(alfa) for alfa in alfaEff])
  gammaNew = 0.5*chords*clWing
  error = np.linalg.norm(gammaNew-gammaOld)
  gammaOld = gammaOld + dampingFactor*(gammaNew-gammaOld)

print('''----------------------------
ITERATION RESULTS
  iterations --> {0}
  final error --> {1}
----------------------------'''.format(iteration,round(error,5)))

#COMPUTE AERODINAMIC DATA
vortexTool = vortexLineTool(gammaOld,wingSpanLine)
vortexTool.getInducedAngle()

surface = pi*0.25*wingSpan*maxChord
AR = (wingSpan**2)/surface
CL = (2/surface) * vortexTool.getCLintegral()
CDi = (2/surface) * vortexTool.getCDintegral()

print('''
------------------------
AERODYNAMIC DATA
  surface --> {0} m^2
  AR --> {1}
  CL --> {2}
  CDi --> {3}
------------------------'''.format(round(surface,3),round(AR,4),round(CL,4),round(CDi,4)))

#PLOT CIRCULATION
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wingSpanLine,gammaNew,'b-',lw=2)
ax.set_title('CIRCULATION')
ax.set_xlabel('wingSpan')
ax.set_ylabel('gamma')
plt.show()