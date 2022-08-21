#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 10:36:23 2018

@author: sergey
"""

import sys;
from PyQt5 import QtGui, QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import pylab as PL
import numpy as NP
import srhArray

class SrhPairCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.subplot = fig.add_subplot(111)
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def imshow(self, array):
        self.subplot.imshow(array)
        self.draw()
    
class SrhAntennaShadowing(QtWidgets.QMainWindow):
    
    def onCalc(self):
        hourAngle = NP.deg2rad(float(self.hourAngle.toPlainText()))
        declination = NP.deg2rad(float(self.declination.toPlainText()))
        frequency = float(self.frequency.toPlainText()) * 1e6
        diameter = float(self.diameter.toPlainText())
        antennaA = int(self.antennaA.toPlainText())
        antennaB = int(self.antennaB.toPlainText())
        self.fillDish(diameter, frequency)
        
        uvw = self.SRH.baseline2uvw(hourAngle, declination, antennaA, antennaB) / (3e8 / frequency)
        u0 = self.M//2
        v0 = self.M//2
        u1 = self.M//2 +int(uvw[1]/self.uvPerPix + .5)
        v1 = self.M//2 -int(uvw[0]/self.uvPerPix + .5)
        
        self.uvPlainSP[:,:] = complex(0,0)
        self.uvPlainCP[:,:] = complex(0,0)
        self.uvPlainSP[u0 - self.N//2:u0 + self.N//2,v0 - self.N//2:v0 + self.N//2] = self.dish*complex(0,-1)
        self.uvPlainSP[u1 - self.N//2:u1 + self.N//2,v1 - self.N//2:v1 + self.N//2] += self.dish
        self.uvPlainCP[u0 - self.N//2:u0 + self.N//2,v0 - self.N//2:v0 + self.N//2] = self.dish
        self.uvPlainCP[u1 - self.N//2:u1 + self.N//2,v1 - self.N//2:v1 + self.N//2] += self.dish*complex(0,-1)
        beamPattern = self.fftBeam(self.fftConvolution(self.uvPlainSP, self.uvPlainCP).real).real
        self.apertureCanvas.imshow(NP.abs(self.uvPlainSP[self.M//2 - self.N:self.M//2 + self.N, self.M//2 - self.N:self.M//2 + self.N]))
        self.patternCanvas.imshow(beamPattern[self.M//2 - self.N:self.M//2 + self.N, self.M//2 - self.N:self.M//2 + self.N] + 5e6*self.qSun[self.M//2 - self.N:self.M//2 + self.N, self.M//2 - self.N:self.M//2 + self.N])
        
    def fillDish(self, diameter, frequency):
        radius = diameter / 2 / (3e8 / frequency) / self.uvPerPix
        for i in range(self.N):
            for j in range(self.N):
                x=i - self.N/2
                y=j - self.N/2
                self.dish[i, j] = 0.
                if (NP.sqrt(x**2 + y**2) < radius):
                    self.dish[i, j] = 1.
                    if (x > 0):
                        self.shadowedDish[i, j] = 1.

    def fillQSun(self, diameter):
        for i in range(self.N):
            for j in range(self.N):
                x=i - self.N/2
                y=j - self.N/2
                self.qSun[i, j] = 0.
                if (NP.sqrt(x**2 + y**2) < diameter/2):
                    self.qSunSmall[i, j] = 1.
        self.qSun[self.M//2 - self.N//2:self.M//2 + self.N//2,self.M//2 - self.N//2:self.M//2 + self.N//2] = self.qSunSmall

    def shift2D(self, arr):
        return NP.roll(arr, (arr.shape[0]//2, arr.shape[0]//2), axis=(0,1))
                       
    def fftConvolution(self, arr1, arr2):
        size = arr1.shape[0]
        return NP.roll((NP.fft.ifft2(NP.fft.fft2(arr1) * NP.conjugate(NP.fft.fft2(arr2)))),(size//2,size//2),axis=(0,1))
        
    def fftBeam(self, uvArr):
        return self.shift2D(NP.fft.fft2(self.shift2D(uvArr)))
        
    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self,parent)
        
        self.N = 128
        self.M = 26 * self.N
        self.meterPerPix = 0.05
        self.pixPerMeter = int(1 / self.meterPerPix + .5)
        self.uvPerPix = 1.0
        self.radPerPix = 1./(self.uvPerPix * self.M)
        self.dish = NP.zeros((self.N, self.N))
        self.shadowedDish = NP.zeros((self.N, self.N))
        self.SRH = srhArray.SrhArray()
        self.uvPlainSP = NP.zeros((self.M, self.M),dtype='complex')
        self.uvPlainCP = NP.zeros((self.M, self.M),dtype='complex')
        self.qSun = NP.zeros((self.M, self.M))
        self.qSunSmall = NP.zeros((self.N, self.N))
        
        self.fillQSun(NP.deg2rad(2160/3600)/self.radPerPix)
        
        self.setGeometry(100,100,800,600)
        self.hourAngleLabel = QtWidgets.QLabel('Hour Angle', self)
        self.declinationLabel = QtWidgets.QLabel('Declination', self)
        self.antennaALabel = QtWidgets.QLabel('AntennaA',self)
        self.antennaBLabel = QtWidgets.QLabel('AntennaB',self)
        self.frequencyLabel = QtWidgets.QLabel('Frequency [MHz]', self)
        self.diameterLabel = QtWidgets.QLabel('Dish diameter', self)
        self.calcButton = QtWidgets.QPushButton('Calc', self)
        
        self.hourAngle = QtWidgets.QTextEdit(self)
        self.declination = QtWidgets.QTextEdit(self)
        self.frequency = QtWidgets.QTextEdit(self)
        self.diameter = QtWidgets.QTextEdit(self)
        self.antennaA = QtWidgets.QTextEdit(self)
        self.antennaB = QtWidgets.QTextEdit(self)

        x0 = 10
        y0 = 10
        dX = 90
        dW = 85
        dY = 15
        dH = 25

        self.hourAngleLabel.setGeometry(x0,y0,dW,20)        
        self.declinationLabel.setGeometry(x0 + dX,y0,dW,20)
        self.frequencyLabel.setGeometry(x0 + 2*dX,y0,dW,20)
        self.diameterLabel.setGeometry(x0 + 3*dX,y0,dW,20)
        self.antennaALabel.setGeometry(x0 + 4*dX,y0,dW,20)
        self.antennaBLabel.setGeometry(x0 + 5*dX,y0,dW,20)
        self.calcButton.setGeometry(x0 + 6*dX,y0,60,40)
        self.calcButton.clicked.connect(self.onCalc)
        
        self.hourAngle.setGeometry(x0,y0 + dY,dW,dH)        
        self.declination.setGeometry(x0 + dX,y0 + dY,dW,dH)
        self.frequency.setGeometry(x0 + 2*dX,y0 + dY,dW,dH)
        self.diameter.setGeometry(x0 + 3*dX,y0 + dY,dW,dH)
        self.antennaA.setGeometry(x0 + 4*dX,y0 + dY,dW,dH)
        self.antennaB.setGeometry(x0 + 5*dX,y0 + dY,dW,dH)
        
        self.hourAngle.setText('0.0')
        self.declination.setText('0.0')
        self.frequency.setText('5900.0')
        self.diameter.setText('1.8')
        self.antennaA.setText('64')
        self.antennaB.setText('192')

        self.apertureCanvas = SrhPairCanvas(self)
        self.apertureCanvas.setGeometry(x0,100,400,400)
        self.patternCanvas = SrhPairCanvas(self)
        self.patternCanvas.setGeometry(x0 + 400,100,400,400)
        
        self.patternCanvas.imshow(self.qSun[self.M//2 - self.N:self.M//2 + self.N, self.M//2 - self.N:self.M//2 + self.N])
        
application = QtWidgets.QApplication.instance();
if not application:
    application = QtWidgets.QApplication(sys.argv);
    
if sys.platform == 'linux':
    font = QtGui.QFont()
    application.setFont(QtGui.QFont(font.defaultFamily(),8));

SRHantShad=SrhAntennaShadowing();
SRHantShad.setWindowTitle('SRH antenna shadowing')
SRHantShad.show();
sys.exit(application.exec_());



