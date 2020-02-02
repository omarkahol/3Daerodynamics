__author__ = "Omar Kahol"
__copyright__ = "Copyright (C) 2020 Omar Kahol"
__license__ = "Public Domain"
__version__ = "1.0"

from math import *
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import simps

#VORTEX LINE TOOL
#A small class used to store methods to perform operations on gamma(Y)
class vortexLineTool:
  def __init__(self, gamma, controlPoints):
    self.gamma = gamma
    self.controlPoints = controlPoints
    self.DY = abs(self.controlPoints[1]-self.controlPoints[0])
    self.nPoints = len(self.gamma)
    self.inducedAngle = []
    self.__DgammaDy__()
    return None
  
  #GET THE DERIVATIVE OF GAMMA
  def __DgammaDy__(self):
    self.DgammaDy=np.gradient(self.gamma,self.DY)
    return False

  #AN IMPLEMENTATION OF SIMPSON'S QUADRATURE FORMULA THAT INTEGRATES THE VORTICITY DISTRIBUTION
  #integral((DgammaDy/(y0-y)))
  #THE INTEGRAL IS UNDEFINED WHEN Y=Y0 
  #FOR THIS REASON IF ONE OF THE THREE TERMS OF SIMPSON'S FORMULA CONTAINS A SINGULARITY
  #THE CODE REPLACES IT WITH THE AVERAGE OF THE OTHER 2
  def getInducedAngle(self):
    alfaIvector=[]
    singularityIndex = 0
    for Y0 in self.controlPoints:
      integral = 0
      #COMPUTE THE INTEGRAL FOR EVERY SECTION OF THE WING
      for i in range(1,self.nPoints-1,1):
        multiplier = (self.DY/(12*pi))
        if i-1==singularityIndex:
          secondTerm = 4*self.DgammaDy[i]/(Y0-self.controlPoints[i])
          thirdTerm = self.DgammaDy[i+1]/(Y0-self.controlPoints[i+1])
          firstTerm = 0.5*(secondTerm+thirdTerm)
        elif i==singularityIndex:
          firstTerm = self.DgammaDy[i-1]/(Y0-self.controlPoints[i-1])
          thirdTerm = self.DgammaDy[i+1]/(Y0-self.controlPoints[i+1])
          secondTerm = 0.5*(firstTerm+thirdTerm)
        elif i+1==singularityIndex:
          firstTerm = self.DgammaDy[i-1]/(Y0-self.controlPoints[i-1])
          secondTerm = 4*self.DgammaDy[i]/(Y0-self.controlPoints[i])
          thirdTerm = 0.5*(firstTerm+secondTerm)
        else:
          firstTerm = self.DgammaDy[i-1]/(Y0-self.controlPoints[i-1])
          secondTerm = 4*self.DgammaDy[i]/(Y0-self.controlPoints[i])
          thirdTerm = self.DgammaDy[i+1]/(Y0-self.controlPoints[i+1])
        integral += multiplier*(firstTerm+secondTerm+thirdTerm)
      alfaIvector.append(integral)
      singularityIndex += 1
    self.inducedAngle = np.array(alfaIvector)
    return np.array(alfaIvector)

  def getCLintegral(self):
    return simps(self.gamma,self.controlPoints)
  
  def getCDintegral(self):
    return simps(self.gamma*self.inducedAngle,self.controlPoints)