from math import *
import numpy as np 
import matplotlib.pyplot as plt 

#WING DATA 
wingSpan = 16
maxChord = 2
chord = lambda y: maxChord*sqrt(1-(2*y/wingSpan)**2)

#PLOT WING
fig = plt.figure()
lim = (-0.5*wingSpan,0.5*wingSpan)
ax = fig.add_subplot(111,aspect='equal',xlim=lim,ylim=lim)
nPoints = 1000
t = np.linspace(0,2*pi,nPoints)
ax.plot(0.5*wingSpan*np.cos(t),0.5*maxChord*np.sin(t))
ax.set_title('WING GEOMETRY')
plt.show()

#SET ANGLE OF ATTACK AND FREESTREAM VELOCITY
alfa = radians(5)
alfaL0 = radians(0)
Vinf = 1

#BUILD SYSTEM OF EQUATIONS
nCoeff = 50
points = np.linspace(-0.5*wingSpan + 1/wingSpan, 0.5*wingSpan-1/wingSpan,nCoeff)
A = np.ones((nCoeff,nCoeff))
b = []
tol = 1e-4
for i in range(nCoeff):
  y = points[i]
  c = chord(y)
  theta = acos(-2*y/wingSpan)
  b.append(alfa-alfaL0)
  for j in range(nCoeff):
    A[i,j] = (2*wingSpan/(pi*c))*sin((j+1)*theta) + (j+1)*sin((j+1)*theta)/sin(theta)

#SOLVE THE SYSTEM
coeff = np.linalg.solve(A,b)

#OBTAIN GAMMA
getGamma = lambda theta: 2*wingSpan*Vinf*sum([a*sin((n+1)*theta) for a,n in zip(coeff, range(nCoeff))])
gamma = []
xaxis = np.linspace(-0.5*wingSpan,0.5*wingSpan,nPoints)
for j in range(nPoints):
  y = xaxis[j]
  theta = acos(-2*y/wingSpan)
  gamma.append(getGamma(theta))

#COMPUTE LIFT AND DRAG
surface = pi*0.5*wingSpan*0.5*maxChord
AR = (wingSpan**2)/surface
CL = pi*AR*coeff[0]

delta = sum([i*(coeff[i]/coeff[0])**2 for i in range(2,nCoeff-2,1)])
efficiency = 1/(1+delta)
CDi = (CL**2)/(pi*AR*efficiency)

#PRINT RESULTS
print('''
----------------------------
WING DATA
  wingspan --> {0} m
  maximum chord --> {1} m
  design --> ellipse
  surface --> {2} m^2
  aspect ratio --> {3}

AERODYNAMIC DATA
  efficiency --> {4}
  calculated CL --> {5}
  calculated CDi --> {6}
-----------------------------
 '''.format(wingSpan,maxChord,round(surface,3),round(AR,4),round(efficiency,4),round(CL,4),round(CDi,4)))

#PLOT SOLUTION
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(xaxis,gamma,'b-',lw=2)
ax.set_title('CIRCULATION')
ax.set_xlabel('wingSpan')
ax.set_ylabel('gamma')
plt.show()