from math import *
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import quad, simps
import urllib.request
import os

#OPEN THE AIRFOIL PROFILE
airFoilURL = 'https://m-selig.ae.illinois.edu/ads/coord/ag18.dat'
airFoilName = 'ag18.dat'
if not os.path.isfile(airFoilName):
  print('downloading file...')
  urllib.request.urlretrieve(airFoilURL, airFoilName)

x,y = [],[]
import pandas as pd 
df = pd.read_table(airFoilName)
for line in df.values:
  values = line[0].split(' ')
  xDone = False
  for value in values:
    if not len(value.strip())>4:
      continue
    elif not xDone:
      x.append(float(value)) 
      xDone = True
    else:
      y.append(float(value))
      break
xBase = np.array(x)
yBase = np.array(y)

#CALCULATE AND FIT THE CAMBER LINE AND SHOW RESULTS
yCamber = 0.5*(yBase + yBase[::-1])
coeff = np.polyfit(xBase, yCamber,4)
getCamber = lambda x: np.polyval(coeff,x)
fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(xBase,yBase,'k-')
xAxis = np.linspace(0,1,100)
yCamberData= [getCamber(x) for x in xAxis]
ax.plot(xAxis,yCamberData,'r--')
plt.show()

#VORTEX PANEL CLASS
class VortexPanel:
  def __init__(self,p1,p2,p3,p4):
    self.X1, self.X2, self.X3, self.X4 = p1[0],p2[0],p3[0],p4[0]
    self.Y1, self.Y2, self.Y3, self.Y4 = p1[1],p2[1],p3[1],p4[1]
    self.Z1, self.Z2, self.Z3, self.Z4 = p1[2],p2[2],p3[2],p4[2]

    self.DY = self.Y2-self.Y1
    self.DX =self.X3-self.X1

    self.tqChord = 0.75
    self.oqChord = 0.23

    self.Xcontrol = 0.5*(self.X1 + self.tqChord*(self.X3-self.X1) + self.X2+self.tqChord*(self.X4-self.X2))
    self.Ycontrol = 0.5*(self.Y1+self.Y2)
    self.Zcontrol = 0.5*(self.Z1 + self.tqChord*(self.Z3-self.Z1) + self.Z2+self.tqChord*(self.Z4-self.Z2))

    self.Yapply1 = self.Y1
    self.Yapply2 = self.Y2

    self.Xapply1 = self.X1 +self.oqChord*(self.X3-self.X1)
    self.Xapply2 = self.X2 +self.oqChord*(self.X4-self.X2)

    self.Zapply1 = self.Z1 +self.oqChord*(self.Z3-self.Z1)
    self.Zapply2 = self.Z2 +self.oqChord*(self.Z4-self.Z2)

    self.phi = 0.5*atan2(self.Z2-self.Z1,self.Y2-self.Y1) + 0.5*atan2(self.Z4-self.Z3, self.Y4-self.Y3)
    self.delta = 0.5*atan2(self.Z3-self.Z1,self.X3-self.X1) + 0.5*atan2(self.Z4-self.Z2, self.X4-self.X2)

    self.surfaceArea = self.DY*self.DX
    return None

  def drawPanel(self,ax):
    ax.plot([self.X1,self.X2],[self.Y1,self.Y2], [self.Z1,self.Z2], 'k-', lw =2)
    ax.plot([self.X2,self.X4],[self.Y2,self.Y4], [self.Z2,self.Z4], 'k-', lw =2)
    ax.plot([self.X4,self.X3],[self.Y4,self.Y3], [self.Z4,self.Z3], 'k-', lw =2)
    ax.plot([self.X3,self.X1],[self.Y3,self.Y1], [self.Z3,self.Z1], 'k-', lw =2)

    ax.plot([self.Xapply1],[self.Yapply1],[self.Zapply1],'b*',lw=1)
    ax.plot([self.Xapply2],[self.Yapply2],[self.Zapply2],'b+',lw=1)

    ax.plot([self.Xcontrol],[self.Ycontrol],[self.Zcontrol],'ro',lw=1)
    return True
  def setGamma(self, gamma):
    self.gamma = gamma

#GEOMETRY DELTA WING AND DISCRETIZATION
wingSpan = 10
maxChord = 10
wingBorder = lambda y: 2*(maxChord/wingSpan)*abs(y)
nPointsY = 11
nPointsX = 5
wingSpanLine = np.linspace(-0.5*wingSpan, 0.5*wingSpan,nPointsY)
panels=[]
for i in range(2,nPointsY-1):
  Y1 = wingSpanLine[i-1]
  Y2 = wingSpanLine[i]

  Y3, Y4 = Y1, Y2
  lim1 = wingBorder(Y1)
  lim2 = wingBorder(Y2)

  xArray1 = np.linspace(lim1, maxChord, nPointsX)
  xArray2 = np.linspace(lim2, maxChord, nPointsX)

  for j in range(1, nPointsX):
    X1 = xArray1[j-1]
    X2 = xArray2[j-1]
    X3 = xArray1[j]
    X4 = xArray2[j]
    Z1 = getCamber((maxChord-X1)/(maxChord-lim1))
    Z2 = getCamber((maxChord-X2)/(maxChord-lim2))
    Z3 = getCamber((maxChord-X3)/(maxChord-lim1))
    Z4 = getCamber((maxChord-X4)/(maxChord-lim2))
    panels.append(VortexPanel((X1,Y1,Z1),(X2,Y2,Z2),(X3,Y3,Z3),(X4,Y4,Z4)))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d',aspect='equal',zlim=(0,0.1))
for panel in panels:
  panel.drawPanel(ax)
plt.show()

#WING DATA
Vinf = 1
alfa = radians(2)

nPoints = len(panels)
velocityMatrix = np.ones((nPoints,nPoints))
rhs = []
#BUILD THE SYSTEM
for m in range(nPoints):
  rhs.append(-Vinf*sin(alfa-panels[m].delta)*cos(panels[m].phi))
  for n in range(nPoints):
    xm = panels[m].Xcontrol
    ym = panels[m].Ycontrol
    zm = panels[m].Zcontrol

    y2n = panels[n].Yapply2
    y1n = panels[n].Yapply1
    x1n = panels[n].Xapply1
    x2n = panels[n].Xapply2
    z1n = panels[n].Zapply1
    z2n = panels[n].Zapply2

    iVector = np.array([1,0,0])
    jVector = np.array([0,1,0])
    kVector = np.array([0,0,1])

    r0 = (x2n-x1n)*iVector + (y2n-y1n)*jVector + (z2n-z1n)*kVector
    r1 = (xm-x1n)*iVector + (ym-y1n)*jVector + (zm-z1n)*kVector
    r2 = (xm-x2n)*iVector + (ym-y2n)*jVector + (zm-z2n)*kVector

    CP1 = (ym-y1n)*(zm-z2n) - (ym -y2n)*(zm-z1n)
    CP2 = (xm-x1n)*(zm-z2n) - (xm-x2n)*(zm-z1n)
    CP3 = (xm-x1n)*(ym-y2n) - (xm-x2n)*(ym-y1n)

    F10 = (0.25/pi)*(CP1*iVector - CP2*jVector + CP3*kVector)/(CP1**2 + CP2**2 + CP3**2)

    T1 = (x2n-x1n)*(xm-x1n) + (y2n-y1n)*(ym-y1n) + (z2n-z1n)*(zm-z1n)
    T2 = (x2n-x1n)*(xm-x2n) + (y2n-y1n)*(ym-y2n) + (z2n-z1n)*(zm-z2n)
    F20 = (T1/sqrt((xm-x1n)**2 + (ym-y1n)**2 + (zm-z1n)**2 )) - (T2/sqrt((xm-x2n)**2 + (ym-y2n)**2 + (zm-z2n)**2))
    VAB = F10*F20

    F11 = (0.25/pi)*((zm-z1n)*jVector +(y1n-ym)*kVector)/((zm-z1n)**2 + (y1n-ym)**2)
    F21 = 1 + (xm-x1n)/sqrt((xm-x1n)**2 + (ym-y1n)**2 + (zm-z1n)**2 )
    VAI = F11*F21

    F12 = (0.25/pi)*((zm-z2n)*jVector +(y2n-ym)*kVector)/((zm-z2n)**2 + (y2n-ym)**2)
    F22 = 1 + (xm-x2n)/sqrt((xm-x2n)**2 + (ym-y2n)**2 + (zm-z2n)**2 )
    VBI = -F12*F22

    VTOT = VAB + VAI + VBI

    phi = panels[n].phi
    delta = panels[n].delta

    velocityMatrix[m,n] = -VTOT[0]*sin(delta)*cos(phi) - VTOT[1]*cos(delta)*sin(phi) + VTOT[2]*cos(phi)*cos(delta)

gammas = np.linalg.solve(velocityMatrix,rhs)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter([panel.Xcontrol for panel in panels],[panel.Ycontrol for panel in panels],gammas)
ax.plot([maxChord,maxChord],[-0.5*wingSpan, 0.5*wingSpan],[0,0], 'k-')
ax.plot([0,maxChord],[0, 0.5*wingSpan],[0,0], 'k-')
ax.plot([0,maxChord],[0,-0.5*wingSpan],[0,0], 'k-')
plt.show()

for gamma, panel in zip(gammas, panels):
  panel.setGamma(gamma)

gammaDict = {}
surfaceArea = 0
for panel in panels:
  key = str(round(panel.Ycontrol,5))
  surfaceArea += panel.surfaceArea
  if key in gammaDict:
    gammaDict[key].append(panel.gamma)
  else:
    gammaDict[key] = [panel.gamma]

sumGamma = []
controlPoints = []
for key in gammaDict:
  sumGamma.append(sum(gammaDict[key]))
  controlPoints.append(float(key))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(controlPoints,sumGamma)
plt.show()

CL = (2/(Vinf*surfaceArea))*simps(sumGamma, controlPoints)
print(CL) 