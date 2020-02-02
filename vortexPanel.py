__author__ = "Omar Kahol"
__copyright__ = "Copyright (C) 2020 Omar Kahol"
__license__ = "Public Domain"
__version__ = "1.0"

from math import *
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#VORTEX PANEL CLASS
#CLASS PANEL
#|   p2---p4 
#y   |     |
#|   p1---p3
#
# ------------- x 
#INPUTS --> 4 points, i.e. the vertices of a panel in 3D
class VortexPanel:
  def __init__(self,p1,p2,p3,p4):
    #SET UP THE COORDINATES OF THE POINTS
    self.X1, self.X2, self.X3, self.X4 = p1[0],p2[0],p3[0],p4[0]
    self.Y1, self.Y2, self.Y3, self.Y4 = p1[1],p2[1],p3[1],p4[1]
    self.Z1, self.Z2, self.Z3, self.Z4 = p1[2],p2[2],p3[2],p4[2]

    #APPROXIMATION OF WIDTH AND HEIGHT OF THE PANEL
    self.DY = self.Y2-self.Y1
    self.DX =self.X3-self.X1

    #THE VORTEX APPLICATION POINT SHOULD BE AT ABOUT 25% OF THE CHORD AND THE CONTROL POINT AT 75%
    #IF THE DISCRETIZATION IS VERY DENSE, NUMERICAL ERRORS CAN ARISE.
    #CHOOSING TWO NUMBERS THAT DO NOT HAVE COMMON DIVISORS SEEMS TO WORK VERY WELL
    self.tqChord = 0.77
    self.oqChord = 0.25

    #CALCULATION OF THE CONTROL POINT AS AN AVERAGE
    self.Xcontrol = 0.5*(self.X1 + self.tqChord*(self.X3-self.X1) + self.X2+self.tqChord*(self.X4-self.X2))
    self.Ycontrol = 0.5*(self.Y1+self.Y2)
    self.Zcontrol = 0.5*(self.Z1 + self.tqChord*(self.Z3-self.Z1) + self.Z2+self.tqChord*(self.Z4-self.Z2))

    #CALCULATION OF THE APPLICATION POINTS
    self.Yapply1 = self.Y1
    self.Yapply2 = self.Y2

    self.Xapply1 = self.X1 +self.oqChord*(self.X3-self.X1)
    self.Xapply2 = self.X2 +self.oqChord*(self.X4-self.X2)

    self.Zapply1 = self.Z1 +self.oqChord*(self.Z3-self.Z1)
    self.Zapply2 = self.Z2 +self.oqChord*(self.Z4-self.Z2)

    #CALCULATE THE ANGE RELATIVE TO THE Y-AXIS AS AN AVERAGE
    self.phi = 0.5*atan2(self.Z2-self.Z1,self.Y2-self.Y1) + 0.5*atan2(self.Z4-self.Z3, self.Y4-self.Y3)

    #CALCULATION OF THE DIHEDRAL ANGLE
    tol = 1e-6
    if abs(self.X1 -self.X3) < tol:
        self.delta = 0.5*atan2(self.Z4-self.Z2, self.X4-self.X2)
    elif abs(self.X4 -self.X2) < tol:
        self.delta = 0.5*atan2(self.Z3-self.Z1,self.X3-self.X1)
    else:
        self.delta = 0.5*atan2(self.Z3-self.Z1,self.X3-self.X1) + 0.5*atan2(self.Z4-self.Z2, self.X4-self.X2)

    self.surfaceArea = self.DY*self.DX
    return None

  def drawPanel3D(self,ax):
    #DRAW THE PANEL
    ax.plot([self.X1,self.X2],[self.Y1,self.Y2], [self.Z1,self.Z2], 'k-', lw =2)
    ax.plot([self.X2,self.X4],[self.Y2,self.Y4], [self.Z2,self.Z4], 'k-', lw =2)
    ax.plot([self.X4,self.X3],[self.Y4,self.Y3], [self.Z4,self.Z3], 'k-', lw =2)
    ax.plot([self.X3,self.X1],[self.Y3,self.Y1], [self.Z3,self.Z1], 'k-', lw =2)

    ax.plot([self.Xapply1],[self.Yapply1],[self.Zapply1],'b*',lw=1)
    ax.plot([self.Xapply2],[self.Yapply2],[self.Zapply2],'b+',lw=1)

    ax.plot([self.Xcontrol],[self.Ycontrol],[self.Zcontrol],'ro',lw=1)
    return True
  
  def drawPanel2D(self,ax):
    #DRAW THE PANEL
    ax.plot([self.X1,self.X2],[self.Y1,self.Y2], 'k-', lw =2)
    ax.plot([self.X2,self.X4],[self.Y2,self.Y4], 'k-', lw =2)
    ax.plot([self.X4,self.X3],[self.Y4,self.Y3], 'k-', lw =2)
    ax.plot([self.X3,self.X1],[self.Y3,self.Y1],'k-', lw =2)

    ax.plot([self.Xapply1],[self.Yapply1],'b*',lw=1)
    ax.plot([self.Xapply2],[self.Yapply2],'b+',lw=1)

    ax.plot([self.Xcontrol],[self.Ycontrol],'ro',lw=1)
    return True

  def setGamma(self, gamma):
    #SET THE VORTEX STRENGTH
    self.gamma = gamma