__author__ = "Omar Kahol"
__copyright__ = "Copyright (C) 2020 Omar Kahol"
__license__ = "Public Domain"
__version__ = "1.0"

import urllib.request
import os
import pandas as pd 
from math import *
import numpy as np
import matplotlib.pyplot as plt

class AeroFoil:
  def __init__(self, airFoilName, dataFolder):
    self.dataFolder = dataFolder
    self.airFoilName = airFoilName
    self.airFoilFile = self.dataFolder+'/'+self.airFoilName+'_coordinates.dat'
    self.airFoilURL = 'https://m-selig.ae.illinois.edu/ads/coord_updates/'+airFoilName.upper()+'.DAT'

    if not os.path.isfile(self.airFoilFile):
      print('downloading airfoil file...')
      try:
        urllib.request.urlretrieve(self.airFoilURL, self.airFoilFile)
      except:
        self.airFoilURL = 'https://m-selig.ae.illinois.edu/ads/coord/'+airFoilName+'.dat'
        urllib.request.urlretrieve(self.airFoilURL, self.airFoilFile)
    self.__readCoords__()
    return None

  def __readCoords__(self):
    self.xBase = []
    self.yBase = []
    with open(self.airFoilFile,'r') as f:
      for line in f:
        values=line.split(' ')
        numbers = [el for el in values if len(el)>3]
        if len(numbers) < 2:
          continue
        try:
          for n in numbers:
            float(n.strip())
        except:
          continue
        if not float(numbers[0].strip()) > 1:
          self.xBase.append(float(numbers[0].strip()))
          self.yBase.append(float(numbers[1].strip()))
    self.xBase=np.array(self.xBase)
    self.yBase=np.array(self.yBase)
    return True

  def fitCamberLine(self, degree):
    yCamber = 0.5*(self.yBase + self.yBase[::-1])
    coeff = np.polyfit(self.xBase, yCamber,degree)
    return (lambda x: np.polyval(coeff,x))

  def downloadPerformance(self):
    self.urlStringPerf = 'http://airfoiltools.com/polar/text?polar=xf-'+self.airFoilName+'-il-1000000'
    self.perfFileName = self.dataFolder+'/'+self.airFoilName+'_performance.txt'

    if not os.path.isfile(self.perfFileName):
      print('downloading airfoil file...')
      urllib.request.urlretrieve(self.urlStringPerf, self.perfFileName)
    self.__readPerf__()
    return True

  def __readPerf__(self):
    self.alfa, self.CL, self.CD = [],[], []
    with open(self.perfFileName,'r') as f:
      for line in f:
        values=line.split(' ')
        numbers = [el for el in values if len(el)>3]
        if len(numbers) < 3:
          continue
        try:
          for n in numbers:
            float(n)
        except:
          continue
        self.alfa.append(float(numbers[0].strip()))
        self.CL.append(float(numbers[1].strip()))
        self.CD.append(float(numbers[2].strip()))
    return True

  def fitCL(self,degree, convertRadians=True):
    alfaFit = np.radians(self.alfa) if convertRadians else self.alfa 
    coeff = np.polyfit(alfaFit, self.CL,degree)
    return (lambda x: np.polyval(coeff,x))

  def fitCD(self,degree, convertRadians=True):
    alfaFit = np.radians(self.alfa) if convertRadians else self.alfa 
    coeff = np.polyfit(alfaFit, self.CD,degree)
    return (lambda x: np.polyval(coeff,x))





