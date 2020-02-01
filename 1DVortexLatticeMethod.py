from math import *
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import quad, simps
import urllib.request
import os
import matplotlib.patches as patches

#CLASS PANEL
#|   p2---p4 
#y   |     |
#|   p1---p3
#
# ------------- x 
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

    self.Xapply1 = self.X1 +0.25*(self.X3-self.X1)
    self.Xapply2 = self.X2 +0.25*(self.X4-self.X2)

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

#GEOMETRY AND WING DISCRETIZATION
wingSpan = 10
chord = 2
baseChord = 1
nPoints = 13
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
 
  panel = VortexPanel((x1,y1),(x2,y2),(x3,y3),(x4,y4))
  panels.append(panel)


#DRAW GEOMETRY
fig = plt.figure()
ax = fig.add_subplot(111,ylim=(-0.5*wingSpan,0.5*wingSpan), xlim=(0,trailingEdgeBound(0.5*wingSpan)))
for panel in panels:
  panel.drawPanel(ax)
plt.show()

#WING DATA
alfa = radians(5)
Vinf = 1

#COMPUTE SYSTEM
velocityMatrix = np.ones((nPoints-1,nPoints-1))
rhs = []
for m in range(nPoints-1):
  rhs.append(-Vinf*sin(alfa))
  for n in range(nPoints-1):
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
ySpanLine = [panel.Ycontrol for panel in panels]
plt.plot(ySpanLine, gammas)
plt.show()

surfaceArea = sum([panel.surfaceArea for panel in panels])
CL = (2/(surfaceArea*Vinf))*simps(gammas,ySpanLine)
print(CL)
