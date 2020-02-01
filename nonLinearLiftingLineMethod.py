from math import *
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import quad, simps

#READ PROFILE FILE
fileName = 'NACA0015.dat'
airFoilData = []
with open(fileName,'r') as file:
  for line in file:
    rawLine = line.split(' ')
    line = [float(el) for el in rawLine if not el=='']
    airFoilData.append(line)

airFoilData = np.array(airFoilData)
alfas = np.radians(airFoilData[:,0])
cl = airFoilData[:,1]

#FIT DATA
coeff = np.polyfit(alfas,cl,5)
getCL = lambda alfa: np.polyval(coeff, alfa)

#PLOT CL ALFA CURVE
fig = plt.figure()
ax = fig.add_subplot(111)
xaxis = np.linspace(-pi/6,pi/6,100)
ax.plot(airFoilData[:,0],cl,'k-',lw=2)
ax.plot(np.degrees(xaxis),[getCL(alfa) for alfa in xaxis],'r--',lw=1)
ax.set_title('CL-ALPHA CURVE FOR NACA0015 AT RE=1E+6')
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
ellipticalCirculation = lambda y: sqrt(1-(2*y/wingSpan)**2)
gammaOld = [ellipticalCirculation(Y0) for Y0 in wingSpanLine]
tol = 1e-4
error = tol + 1
itmax = 1000
iteration = 0

#GAMMA FUNCTION
def getAlphaIvector(dGamma_dy):
  dy = wingSpanLine[1]-wingSpanLine[0]
  singularityIndex=0
  alfaIvector = []
  #CYCLE THROUGH EVERY SECTION OF THE WING
  for Y0 in wingSpanLine:
    integral = 0
    #COMPUTE THE INTEGRAL FOR EVERY SECTION OF THE WING
    for i in range(1,nPoints-1,1):
      multiplier = (dy/(12*pi))
      if i-1==singularityIndex:
        secondTerm = 4*dGamma_dy[i]/(Y0-wingSpanLine[i])
        thirdTerm = dGamma_dy[i+1]/(Y0-wingSpanLine[i+1])
        firstTerm = 0.5*(secondTerm+thirdTerm)
      elif i==singularityIndex:
        firstTerm = dGamma_dy[i-1]/(Y0-wingSpanLine[i-1])
        thirdTerm = dGamma_dy[i+1]/(Y0-wingSpanLine[i+1])
        secondTerm = 0.5*(firstTerm+thirdTerm)
      elif i+1==singularityIndex:
        firstTerm = dGamma_dy[i-1]/(Y0-wingSpanLine[i-1])
        secondTerm = 4*dGamma_dy[i]/(Y0-wingSpanLine[i])
        thirdTerm = 0.5*(firstTerm+secondTerm)
      else:
        firstTerm = dGamma_dy[i-1]/(Y0-wingSpanLine[i-1])
        secondTerm = 4*dGamma_dy[i]/(Y0-wingSpanLine[i])
        thirdTerm = dGamma_dy[i+1]/(Y0-wingSpanLine[i+1])
      integral += multiplier*(firstTerm+secondTerm+thirdTerm)
    alfaIvector.append(integral)
    singularityIndex += 1
  return alfaIvector


#START ITERATIONS
while (error > tol) and (iteration < itmax):
  iteration = iteration + 1
  dGamma_dy = np.gradient(gammaOld,wingSpanLine[1]-wingSpanLine[0])
  alfaI = getAlphaIvector(dGamma_dy)
  alfaEff = alphaWing - np.array(alfaI)
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
surface = pi*0.25*wingSpan*maxChord
AR = (wingSpan**2)/surface
CL = (2/surface) * simps(gammaNew,wingSpanLine)
CDi = (2/surface) * simps(gammaNew*alfaI,wingSpanLine)

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