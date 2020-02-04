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
  def __init__(self, gamma, controlPoints, Vinf):
    self.gamma = gamma
    self.Vinf = Vinf
    self.controlPoints = controlPoints
    self.DY = abs(self.controlPoints[1]-self.controlPoints[0])
    self.nPoints = len(self.gamma)
    self.inducedAngle = []
    self.__DgammaDy__()
    return None
  
  #GET THE DERIVATIVE OF GAMMA
  def __DgammaDy__(self):
    self.DgammaDy=np.gradient(self.gamma,self.controlPoints)
    return False

  #AN IMPLEMENTATION OF SIMPSON'S QUADRATURE FORMULA THAT INTEGRATES THE VORTICITY DISTRIBUTION
  #integral((DgammaDy/(y0-y)))
  #THE INTEGRAL IS UNDEFINED WHEN Y=Y0 
  #THAT TERM WILL BE REPLACE WITH THE SUM OF THE ADIACENT TWO AS SUGGESTED BY JOHN ANDERSONS'S FOUNDAMENTALS OF AERODYNAMICS BOOK
  def getInducedAngle(self):
    self.inducedAngle = []
    multiplier = 1/(4*pi*self.Vinf)
    singularityIndex = 0
    for Y0 in self.controlPoints:

      #COMPUTE DENOMINATORS AND REPLACE THE NULL WITH 1
      denominators = Y0 - np.array(self.controlPoints)
      denominators[singularityIndex]=1

      #NOW COMPUTE THE INTEGRAND FUNCTION
      function = self.DgammaDy/denominators

      #REPLACE WRONG TERM
      if singularityIndex == 0:
        function[singularityIndex]=0
      elif singularityIndex == self.nPoints-1:
        function[singularityIndex]=0
      else:
        function[singularityIndex] = 0.5*(function[singularityIndex-1] + function[singularityIndex+1])

      #COMPUTE THE INTEGRAL USING SISMPSON'S RULE
      integral = simps(function,self.controlPoints)
      
      #INCREASE THE INDEX AT WHICH THE SINGULARITY OCCURS
      singularityIndex += 1

      #APPEND RESULT
      self.inducedAngle.append(multiplier*integral)

    #COVERT TO NUMPY ARRAY
    self.inducedAngle = np.array(self.inducedAngle)
    return np.copy(self.inducedAngle)

  def getCLintegral(self):
    return simps(self.gamma,self.controlPoints)
  
  def getCDintegral(self):
    return simps(self.gamma*self.inducedAngle,self.controlPoints)