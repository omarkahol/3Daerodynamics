from math import *
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import quad, simps
import urllib.request
import os
import matplotlib.patches as patches

class VortexPanel:
  def __init__(self,p1,p2,p3,p4):
    self.X1, self.X2, self.X3, self.X4 = p1[0],p2[0],p3[0],p4[0]
    self.Y1, self.Y2, self.Y3, self.Y4 = p1[1],p2[1],p3[1],p4[1]
    

    self.DY = self.Y2-self.Y1
    self.DX =self.X3-self.X1

    self.Xcontrol = 0.5*(self.X1 + 0.75*(self.X3-self.X1) + self.X2+0.75*(self.X4-self.X2))
    self.Ycontrol = 0.5*(self.Y1+self.Y2)

    self.Yapply1 = self.Y1
    self.Yapply2 = self.Y2

    self.Xapply1 = self.X1 +0.23*(self.X3-self.X1)
    self.Xapply2 = self.X2 +0.23*(self.X4-self.X2)

    self.surfaceArea = self.DY*self.DX
    return None

  def drawPanel(self,ax):
    ax.plot([self.X1,self.X2],[self.Y1,self.Y2], 'k-', lw =2)
    ax.plot([self.X2,self.X4],[self.Y2,self.Y4], 'k-', lw =2)
    ax.plot([self.X4,self.X3],[self.Y4,self.Y3], 'k-', lw =2)
    ax.plot([self.X3,self.X1],[self.Y3,self.Y1], 'k-', lw =2)

    ax.plot([self.Xapply1],[self.Yapply1],'b*',lw=1)
    ax.plot([self.Xapply2],[self.Yapply2],'b+',lw=1)

    ax.plot([self.Xcontrol],[self.Ycontrol],'ro',lw=1)
    return True
  def setGamma(self, gamma):
    self.gamma = gamma

#GEOMETRY DELTA WING AND DISCRETIZATION
wingSpan = 10
maxChord = 10
wingBorder = lambda y: 2*(maxChord/wingSpan)*abs(y)
nPointsY = 101
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
    panels.append(VortexPanel((X1,Y1),(X2,Y2),(X3,Y3),(X4,Y4)))

#DRAW GEOMETRY
fig = plt.figure()
ax = fig.add_subplot(111,ylim=(-0.5*wingSpan,0.5*wingSpan), xlim=(0,maxChord))
for panel in panels:
  panel.drawPanel(ax)
plt.show()

#WING DATA
alfa = radians(2)
Vinf = 1

#COMPUTE SYSTEM
nPoints = len(panels)
velocityMatrix = np.ones((nPoints,nPoints))
rhs = []
for m in range(nPoints):
  rhs.append(-Vinf*sin(alfa))
  for n in range(nPoints):
    xm = panels[m].Xcontrol
    ym = panels[m].Ycontrol

    y2n = panels[n].Yapply2
    y1n = panels[n].Yapply1
    x1n = panels[n].Xapply1
    x2n = panels[n].Xapply2

    multiplier1 = 1/((xm-x1n)*(ym-y2n)-(xm-x2n)*(ym-y1n))
    firstTerm1 = ((x2n-x1n)*(xm-x1n) + (y2n - y1n)*(ym-y1n))/sqrt((xm-x1n)**2 + (ym-y1n)**2)
    firstTerm2 = ((x2n-x1n)*(xm-x2n) + (y2n - y1n)*(ym-y2n))/sqrt((xm-x2n)**2 + (ym-y2n)**2)

    multiplier2 = 1/(y1n-ym)
    secondTerm = 1 + (xm-x1n)/sqrt((xm-x1n)**2 + (ym-y1n)**2)

    multiplier3 = 1/(y2n-ym)
    thirdTerm = 1 + (xm-x2n)/sqrt((xm-x2n)**2 + (ym-y2n)**2)

    wMN = (0.25/pi)*(multiplier1*(firstTerm1-firstTerm2) + multiplier2*secondTerm - multiplier3*thirdTerm)
    velocityMatrix[m,n] = wMN

gammas = np.linalg.solve(velocityMatrix, rhs)

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