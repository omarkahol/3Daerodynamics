__author__ = "Omar Kahol"
__copyright__ = "Copyright (C) 2020 Omar Kahol"
__license__ = "Public Domain"
__version__ = "1.0"

from math import *
import numpy as np 
from gammaTool import vortexLineTool

class vortexLatticeSolver:
  def __init__(self,panels):
    self.panels = panels
    self.surfaceArea = 0
    self.alfa = 0
    self.Vinf = 0
    self.nPoints = len(self.panels)
    self.velocityMatrix = np.ones((self.nPoints,self.nPoints))
    self.RHS = np.ones(self.nPoints)
    self.iVector = np.array([1,0,0])
    self.jVector = np.array([0,1,0])
    self.kVector = np.array([0,0,1])
    self.controlPoints=np.ones(self.nPoints)
    self.leadingEdgeGamma = np.ones(self.nPoints)

  def setDynamicData(self,alfa, Vinf, isDegrees=True):
    self.alfa = radians(alfa) if isDegrees else alfa
    self.Vinf = abs(Vinf)
    self.__buildSystem__()
    return True

  def resetPanels(self, panels):
    self.panels = panels
    self.nPoints = len(self.panels)
    self.controlPoints=np.ones(self.nPoints)
    self.leadingEdgeGamma = np.ones(self.nPoints)
    self.__buildSystem__()
    return True

  def __buildSystem__(self):
    self.controlPoints=[]
    self.leadingEdgeGamma = []
    for m in range(self.nPoints):
      self.RHS[m]=-self.Vinf*sin(self.alfa-self.panels[m].delta)*cos(self.panels[m].phi)
      for n in range(self.nPoints):
        xm = self.panels[m].Xcontrol
        ym = self.panels[m].Ycontrol
        zm = self.panels[m].Zcontrol

        y2n = self.panels[n].Yapply2
        y1n = self.panels[n].Yapply1
        x1n = self.panels[n].Xapply1
        x2n = self.panels[n].Xapply2
        z1n = self.panels[n].Zapply1
        z2n = self.panels[n].Zapply2

        r0 = (x2n-x1n)*self.iVector + (y2n-y1n)*self.jVector + (z2n-z1n)*self.kVector
        r1 = (xm-x1n)*self.iVector + (ym-y1n)*self.jVector + (zm-z1n)*self.kVector
        r2 = (xm-x2n)*self.iVector + (ym-y2n)*self.jVector + (zm-z2n)*self.kVector

        # VORTEX STRUCTURE
        #
        # A-------B
        # |       |
        # |       |
        # |       |
        #INF     INF

        #CALCULATE THE INDUCED VELOCITY BY AB SEGMENT 
        F10 = np.cross(r1,r2)/(np.linalg.norm(np.cross(r1,r2))**2)
        F20 = np.dot(r0, (r1/np.linalg.norm(r1))-(r2/np.linalg.norm(r2)))
        VAB = (0.25/pi)*F10*F20

        #CALCULATE INDUCED VELOCITY BY A-INF SEGMENT
        F11 = (0.25/pi)*((zm-z1n)*self.jVector +(y1n-ym)*self.kVector)/((zm-z1n)**2 + (y1n-ym)**2)
        F21 = 1 + (xm-x1n)/np.linalg.norm(r1)
        VAI = F11*F21

        #CALCULATE INDUCED VELOCITY BY B-INF SEGMENT
        F12 = (0.25/pi)*((zm-z2n)*self.jVector +(y2n-ym)*self.kVector)/((zm-z2n)**2 + (y2n-ym)**2)
        F22 = 1 + (xm-x2n)/np.linalg.norm(r2)
        VBI = -F12*F22

        #TOTAL INDUCED VELOCITY
        VTOT = VAB + VAI + VBI

        #GEOMETRIC ANGLES
        phi = self.panels[n].phi
        delta = self.panels[n].delta
        nVector = np.array([-sin(delta)*cos(phi),-cos(delta)*sin(phi),cos(phi)*cos(delta)])

        self.velocityMatrix[m,n] = np.dot(VTOT, nVector)
    return True

  def solveSystem(self):
    self.solution = np.linalg.solve(self.velocityMatrix, self.RHS)
    gammaDict = {}
    for i, gamma in enumerate(self.solution):
      self.panels[i].setGamma(gamma)
      self.surfaceArea += self.panels[i].surfaceArea
      key = str(round(self.panels[i].Ycontrol,5))
      if key in gammaDict:
        gammaDict[key].append(self.panels[i].gamma)
      else:
        gammaDict[key] = [self.panels[i].gamma]
    for key in gammaDict:
      self.leadingEdgeGamma.append(sum(gammaDict[key]))
      self.controlPoints.append(float(key))
    return True