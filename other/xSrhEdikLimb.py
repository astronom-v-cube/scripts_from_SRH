#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 10:36:23 2018

@author: sergey
"""

import sys;
from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtWidgets import QFileDialog
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as NP
from srhFitsFile import SrhFitsFile
from skimage.transform import warp, AffineTransform
from astropy.io import fits
from astropy.time import Time, TimeDelta
import casacore.tables as T
import base2uvw as bl2uvw
import srhMS2
from mpl_toolkits.mplot3d import Axes3D
import pylab as PL
from scipy.signal import argrelextrema
from sunpy import coordinates
import matplotlib.animation as animation

class ResponseCanvas(FigureCanvas):
    mouseSignal = QtCore.pyqtSignal(float, float, name = 'xyChanged')
    
    def mouse_moved(self, event):
        1
#        self.mouseSignal.emit(event.xdata, event.ydata)
        
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.subplot = self.fig.add_subplot(111)
        self.subplot.xaxis.set_visible(False)
        self.subplot.yaxis.set_visible(False)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.mpl_connect('motion_notify_event', self.mouse_moved)
#        self.mpl_connect('button_release_event', self.mouse_moved)
        self.cmap = 'ocean'
        
    def setData(self, array):
        self.imageObject.set_data(array)
        self.draw()

    def setCurveXY(self, curveX, curveY):
        self.plotObject[0].setx_data(curveX)
        self.plotObject[0].sety_data(curveY)
        self.draw()

    def imshow(self, array, arrayMin, arrayMax):
        self.imageObject = self.subplot.imshow(array, vmin = arrayMin, vmax = arrayMax, cmap=self.cmap, origin='lower')
        self.draw()

    def contour(self, array, levels):
        self.contourObject = self.subplot.contour(array, levels)
        self.draw()

    def plot(self, data):
        self.plotObject = self.subplot.plot(data)
        self.draw()
    
    def scatter(self, data):
        self.plotObject = self.subplot.plot(data, '.')
        self.draw()
    
    def clear(self):
        self.subplot.cla()
        self.draw()
    
    def redraw(self):
        self.draw()
    
    def setColormap(self, cmap):
        self.cmap = cmap

class SrhEdik(QtWidgets.QMainWindow):
#    def buildEwPhase(self):
#        self.ewLcpPhaseCorrection[:] = 0.
#        self.ewRcpPhaseCorrection[:] = 0.
#        for j in range(16):
#            for i in range(16):
#                self.ewLcpPhaseCorrection[16 + i] +=  NP.deg2rad(self.ewPhaseCoefsLcp[self.currentFrequencyChannel, j] * (-1)**(i // (j + 1))) 
#                self.ewLcpPhaseCorrection[15 - i]  += -NP.deg2rad(self.ewPhaseCoefsLcp[self.currentFrequencyChannel, j] * (-1)**(i // (j + 1))) 
#                self.ewRcpPhaseCorrection[16 + i] +=  NP.deg2rad(self.ewPhaseCoefsRcp[self.currentFrequencyChannel, j] * (-1)**(i // (j + 1))) 
#                self.ewRcpPhaseCorrection[15 - i]  += -NP.deg2rad(self.ewPhaseCoefsRcp[self.currentFrequencyChannel, j] * (-1)**(i // (j + 1))) 
#        for j in range(32):
#                self.ewLcpPhaseCorrection[j] +=  NP.deg2rad(self.ewLcpPhaseSlope[self.currentFrequencyChannel] * (j - 16)) 
#                self.ewRcpPhaseCorrection[j] +=  NP.deg2rad(self.ewRcpPhaseSlope[self.currentFrequencyChannel] * (j - 16)) 
#        self.srhFits.changeEastWestPhase(self.ewLcpPhaseCorrection, self.ewRcpPhaseCorrection)
#        
#    def buildSPhase(self):
#        self.sLcpPhaseCorrection[:] = 0.
#        self.sRcpPhaseCorrection[:] = 0.
#        for j in range(16):
#            for i in range(16):
#                self.sLcpPhaseCorrection[i] +=  NP.deg2rad(self.sPhaseCoefsLcp[self.currentFrequencyChannel, j] * (-1)**(i // (j + 1))) 
#                self.sRcpPhaseCorrection[i] +=  NP.deg2rad(self.sPhaseCoefsRcp[self.currentFrequencyChannel, j] * (-1)**(i // (j + 1))) 
#        for j in range(16):
#                self.sLcpPhaseCorrection[j] +=  NP.deg2rad(self.sLcpPhaseSlope[self.currentFrequencyChannel] * j) 
#                self.sRcpPhaseCorrection[j] +=  NP.deg2rad(self.sRcpPhaseSlope[self.currentFrequencyChannel] * j) 
#        self.srhFits.changeSouthPhase(self.sLcpPhaseCorrection, self.sRcpPhaseCorrection)
#
#
    def buildEwPhase(self):
        self.ewLcpPhaseCorrection[:] = 0.
        self.ewRcpPhaseCorrection[:] = 0.
        for j in range(16):
            for i in range(16):
                self.ewLcpPhaseCorrection[16 + i] +=  NP.deg2rad(self.ewPhaseCoefsLcp[self.currentFrequencyChannel, j] * (-1)**(i // (j + 1))) 
                self.ewLcpPhaseCorrection[15 - i]  += -NP.deg2rad(self.ewPhaseCoefsLcp[self.currentFrequencyChannel, j] * (-1)**(i // (j + 1))) 
                self.ewRcpPhaseCorrection[16 + i] +=  NP.deg2rad(self.ewPhaseCoefsRcp[self.currentFrequencyChannel, j] * (-1)**(i // (j + 1))) 
                self.ewRcpPhaseCorrection[15 - i]  += -NP.deg2rad(self.ewPhaseCoefsRcp[self.currentFrequencyChannel, j] * (-1)**(i // (j + 1))) 
        for j in range(32):
                self.ewLcpPhaseCorrection[j] += NP.deg2rad(self.ewLcpPhaseSlope[self.currentFrequencyChannel] * (j - 15.5)) 
                self.ewRcpPhaseCorrection[j] += NP.deg2rad(self.ewRcpPhaseSlope[self.currentFrequencyChannel] * (j - 15.5))
                self.ewLcpPhaseCorrection[j] += NP.deg2rad(self.ewLcpPhaseAnt[j])
                self.ewRcpPhaseCorrection[j] += NP.deg2rad(self.ewRcpPhaseAnt[j])
        self.srhFits.changeEastWestPhase(self.ewLcpPhaseCorrection, self.ewRcpPhaseCorrection)
        
    def buildSPhase(self):
        self.sLcpPhaseCorrection[:] = 0.
        self.sRcpPhaseCorrection[:] = 0.
        for j in range(16):
            for i in range(16):
                self.sLcpPhaseCorrection[i] +=  NP.deg2rad(self.sPhaseCoefsLcp[self.currentFrequencyChannel, j] * (-1)**(i // (j + 1))) 
                self.sRcpPhaseCorrection[i] +=  NP.deg2rad(self.sPhaseCoefsRcp[self.currentFrequencyChannel, j] * (-1)**(i // (j + 1))) 
        for j in range(16):
                self.sLcpPhaseCorrection[j] += NP.deg2rad(self.sLcpPhaseSlope[self.currentFrequencyChannel] * (j + .5)) 
                self.sRcpPhaseCorrection[j] += NP.deg2rad(self.sRcpPhaseSlope[self.currentFrequencyChannel] * (j + .5)) 
                self.sLcpPhaseCorrection[j] += NP.deg2rad(self.sLcpPhaseAnt[j])
                self.sRcpPhaseCorrection[j] += NP.deg2rad(self.sRcpPhaseAnt[j])
        self.srhFits.changeSouthPhase(self.sLcpPhaseCorrection, self.sRcpPhaseCorrection)
 
    def onFindPhase(self):
        self.lcpMaxTrace = []
        self.rcpMaxTrace = []
        self.srhFits.setSizeOfUv(256)
        for sPhaseInd in range(-18,18):
            self.sLcpPhaseCorrection[:] = NP.deg2rad(sPhaseInd*10)
            self.sRcpPhaseCorrection[:] = NP.deg2rad(sPhaseInd*10)
            self.srhFits.changeSouthPhase(self.sLcpPhaseCorrection, self.sRcpPhaseCorrection)
            self.srhFits.vis2uv(self.currentScan, phaseCorrect=self.phaseCorrect, amplitudeCorrect=self.amplitudeCorrect);
            self.srhFits.uv2lmImage()
            self.lcpMaxTrace.append(self.srhFits.lcp.real[128-32:128+32,128-32:128+32].mean())
            self.rcpMaxTrace.append(self.srhFits.rcp.real[128-32:128+32,128-32:128+32].mean())
        self.srhFits.setSizeOfUv(self.uvSize)
        
        phaseIndLcp = int(10*(NP.argmax(self.lcpMaxTrace) - 18) + .5)
        phaseIndRcp = int(10*(NP.argmax(self.rcpMaxTrace) - 18) + .5)
  
        self.sPhaseCoefsLcp[self.currentFrequencyChannel, 15] = phaseIndLcp
        self.sPhaseCoefsRcp[self.currentFrequencyChannel, 15] = phaseIndRcp
        self.sPhaseStairLcp.setValue(phaseIndLcp)
        self.sPhaseStairRcp.setValue(phaseIndRcp)
        self.lcpMaxCanvas.plot(self.lcpMaxTrace)
        self.rcpMaxCanvas.plot(self.rcpMaxTrace)
        
    def onFindMinLcpRcp(self):
        self.lcpMaxTrace = []
        self.srhFits.setSizeOfUv(256)
        for sPhaseSlopeInd in range(-60,60):
            self.sRcpPhaseSlope[self.currentFrequencyChannel] = sPhaseSlopeInd
            self.buildSPhase()
#            self.srhFits.vis2uv(self.currentScan, phaseCorrect=self.phaseCorrect, amplitudeCorrect=self.amplitudeCorrect, PSF=False);
            self.srhFits.vis2uv(self.currentScan, phaseCorrect=self.phaseCorrect, amplitudeCorrect=self.amplitudeCorrect);
            self.srhFits.uv2lmImage()
            self.lcpMaxTrace.append(((self.srhFits.lcp.real - self.srhFits.rcp.real)**2).mean())
        self.sRcpPhaseSlope[self.currentFrequencyChannel] = NP.array(self.lcpMaxTrace).argmin() - 60
        self.sRcpPhaseSlopeSpin.setValue(self.sRcpPhaseSlope[self.currentFrequencyChannel])

#        self.lcpMaxTrace = []
#        for ewPhaseSlopeInd in range(-60,60):
#            self.ewRcpPhaseSlope[self.currentFrequencyChannel] = ewPhaseSlopeInd
#            self.buildEwPhase()
#            self.srhFits.vis2uv(self.currentScan, phaseCorrect=self.phaseCorrect, amplitudeCorrect=self.amplitudeCorrect, PSF=False);
#            self.srhFits.uv2lmImage()
#            self.lcpMaxTrace.append(((self.srhFits.lcp.real - self.srhFits.rcp.real)**2).mean())
#        self.ewRcpPhaseSlope[self.currentFrequencyChannel] = NP.array(self.lcpMaxTrace).argmin() - 60
#        self.ewRcpPhaseSlopeSpin.setValue(self.ewRcpPhaseSlope[self.currentFrequencyChannel])

        self.srhFits.setSizeOfUv(self.uvSize)       

    def buildImage(self):
        if self.indexOfImageType == 3:
            self.srhFits.vis2uv(self.currentScan, phaseCorrect=self.phaseCorrect, amplitudeCorrect=self.amplitudeCorrect, PSF=True);
        elif self.indexOfImageType == 4:
            self.srhFits.vis2uv(self.currentScan, phaseCorrect=self.phaseCorrect, amplitudeCorrect=self.amplitudeCorrect, PSF=True);
            self.srhFits.uv2lmImage()
            self.psfData = self.srhFits.lcp.real
            self.srhFits.vis2uv(self.currentScan, phaseCorrect=self.phaseCorrect, amplitudeCorrect=self.amplitudeCorrect, PSF=False);
        else:
            self.srhFits.vis2uv(self.currentScan, phaseCorrect=self.phaseCorrect, amplitudeCorrect=self.amplitudeCorrect, PSF=False);
        self.srhFits.uv2lmImage()

        scaling = self.srhFits.getPQScale(self.uvSize, NP.deg2rad(self.arcsecPerPixel * (self.uvSize - 1)/3600.))
        scale = AffineTransform(scale=(self.uvSize/scaling[0], self.uvSize/scaling[1]))
        shift = AffineTransform(translation=(-self.uvSize/2,-self.uvSize/2))
        rotate = AffineTransform(rotation = self.pAngle)
        matrix = AffineTransform(matrix = self.srhFits.getPQ2HDMatrix())
        back_shift = AffineTransform(translation=(self.uvSize/2,self.uvSize/2))

        if self.indexOfFrameType == 0:
            self.lcpData = self.srhFits.lcp.real
            self.rcpData = self.srhFits.rcp.real
        else:
            dataResult0 = warp(self.srhFits.lcp.real,(shift + (scale + back_shift)).inverse)
            self.lcpData = warp(dataResult0,(shift + (matrix + back_shift)).inverse)
            dataResult0 = warp(self.srhFits.rcp.real,(shift + (scale + back_shift)).inverse)
            self.rcpData = warp(dataResult0,(shift + (matrix + back_shift)).inverse)
            if self.indexOfImageType == 4:
                dataResult0 = warp(self.psfData,(shift + (scale + back_shift)).inverse)
                self.psfData = warp(dataResult0,(shift + (matrix + back_shift)).inverse)
            dataResult0 = 0
            if self.indexOfFrameType == 2:
                self.lcpData = warp(self.lcpData,(shift + (rotate + back_shift)).inverse)
                self.rcpData = warp(self.rcpData,(shift + (rotate + back_shift)).inverse)

        self.lcpData = NP.flip(self.lcpData,0)
        self.rcpData = NP.flip(self.rcpData,0)
        
        if self.indexOfImageType == 4:
            self.lcpData += self.psfData  / NP.max(self.psfData) * NP.max(self.lcpData) * .5
            self.rcpData += self.psfData  / NP.max(self.psfData) * NP.max(self.rcpData) * .5

        if self.imageCalibrateButton.isChecked():
            lcpHist = NP.histogram(self.lcpData, bins=100)
            sky_ind = NP.argmax(lcpHist[0])
            sun_ind = (argrelextrema(lcpHist[0][sky_ind:], NP.greater, order=5)[0] + sky_ind)[0]
            sky_level = lcpHist[1][sky_ind + 1]
            sun_level = lcpHist[1][sun_ind + 1]
            print('Levels')
            print(sky_level, sun_level)
            self.lcpData = (self.lcpData - sky_level)/(sun_level - sky_level) * 15
            rcpHist = NP.histogram(self.rcpData, bins=100)
            sky_ind = NP.argmax(rcpHist[0])
            sun_ind = (argrelextrema(rcpHist[0][sky_ind:], NP.greater, order=5)[0] + sky_ind)[0]
            sky_level = rcpHist[1][sky_ind + 1]
            sun_level = rcpHist[1][sun_ind + 1]
            self.rcpData = (self.rcpData - sky_level)/(sun_level - sky_level) * 15
            
        if self.indexOfImageType == 0 or self.indexOfImageType == 3 or self.indexOfImageType == 4:
            self.lcpCanvas.setData(self.lcpData*self.imageScale + self.imageOffset)
            self.rcpCanvas.setData(self.rcpData*self.imageScale + self.imageOffset)
            if self.indexOfImageType == 4:
                self.lcpCanvas.contour(self.psfData*self.imageScale + self.imageOffset, NP.max(self.psfData)*.5*self.imageScale + self.imageOffset)
        elif self.indexOfImageType == 1:
            self.lcpCanvas.setData((self.lcpData + self.lcpRcpRel * self.rcpData)*.5*self.imageScale + self.imageOffset)
            self.rcpCanvas.setData((self.lcpData -  self.lcpRcpRel * self.rcpData)*.5*self.imageScale + self.imageOffset)
        else:
            self.lcpCanvas.setData(NP.abs(self.srhFits.uvLcp)**.5 *self.imageScale + self.imageOffset)
            self.rcpCanvas.setData(NP.abs(self.srhFits.uvRcp)**.5 *self.imageScale + self.imageOffset)
#            self.lcpCanvas.setData(curI*.5*self.imageScale + self.imageOffset)
#            self.rcpCanvas.setData((NP.abs( curI))*.5*self.imageScale + self.imageOffset)
        if (self.imageLimbButton.isChecked()):
            self.lcpCanvas.contour(self.qSun, levels=[0.5])

        self.lcpMaxTrace.append(NP.max(self.srhFits.lcp.real))
        self.rcpMaxTrace.append(NP.max(self.srhFits.rcp.real))
#        self.lcpMaxTrace.append(NP.max(self.srhFits.lcp.real) - NP.min(self.srhFits.lcp.real))
#        self.rcpMaxTrace.append(NP.max(self.srhFits.rcp.real) - NP.min(self.srhFits.rcp.real))
#        self.lcpMaxTrace.append(self.srhFits.lcp.real[128-32:128+32,128-32:128+32].mean())
#        self.rcpMaxTrace.append(self.srhFits.rcp.real[128-32:128+32,128-32:128+32].mean())

        self.lcpMaxCanvas.clear()
        self.rcpMaxCanvas.clear()
        if self.imageCenteringButton.isChecked():
            self.lcpMaxCanvas.scatter(self.lcpMaxTrace)
            self.rcpMaxCanvas.scatter(self.rcpMaxTrace)
        else:
            self.lcpMaxCanvas.plot(self.lcpMaxTrace)
            self.rcpMaxCanvas.plot(self.rcpMaxTrace)
        
    def onEastWestLcpPhaseSlopeChanged(self, value):
        self.ewLcpPhaseSlope[self.currentFrequencyChannel] = value
        self.buildEwPhase()
        if (self.imageUpdate):
            self.buildImage()

    def onEastWestRcpPhaseSlopeChanged(self, value):
        self.ewRcpPhaseSlope[self.currentFrequencyChannel] = value
        self.buildEwPhase()
        if (self.imageUpdate):
            self.buildImage()

    def onSouthLcpPhaseSlopeChanged(self, value):
        self.sLcpPhaseSlope[self.currentFrequencyChannel] = value
        self.buildSPhase()
        if (self.imageUpdate):
            self.buildImage()

    def onSouthRcpPhaseSlopeChanged(self, value):
        self.sRcpPhaseSlope[self.currentFrequencyChannel] = value
        self.buildSPhase()
        if (self.imageUpdate):
            self.buildImage()

    def onEastWestPhaseStairLcpChanged(self, value):
        self.ewPhaseCoefsLcp[self.currentFrequencyChannel, self.ewStairLength - 1] = value
        self.buildEwPhase()
        if (self.imageUpdate):
            self.buildImage()

    def onEastWestPhaseStairRcpChanged(self, value):
        self.ewPhaseCoefsRcp[self.currentFrequencyChannel, self.ewStairLength - 1] = value
        self.buildEwPhase()
        if (self.imageUpdate):
            self.buildImage()

    def onSouthPhaseStairLcpChanged(self, value):
        self.sPhaseCoefsLcp[self.currentFrequencyChannel, self.sStairLength - 1] = value
        self.buildSPhase()
        if (self.imageUpdate):
            self.buildImage()
        
    def onSouthPhaseStairRcpChanged(self, value):
        self.sPhaseCoefsRcp[self.currentFrequencyChannel, self.sStairLength - 1] = value
        self.buildSPhase()
        if (self.imageUpdate):
            self.buildImage()
        
    def onEwPhaseStairLengthChanged(self, value):
        self.ewStairLength = value
        self.ewPhaseStairLcp.setValue(self.ewPhaseCoefsLcp[self.currentFrequencyChannel, self.ewStairLength - 1])
        self.ewPhaseStairRcp.setValue(self.ewPhaseCoefsRcp[self.currentFrequencyChannel, self.ewStairLength - 1])

    def onSPhaseStairLengthChanged(self, value):
        self.sStairLength = value
        self.sPhaseStairLcp.setValue(self.sPhaseCoefsLcp[self.currentFrequencyChannel, self.sStairLength - 1])
        self.sPhaseStairRcp.setValue(self.sPhaseCoefsRcp[self.currentFrequencyChannel, self.sStairLength - 1])

    def onClear(self):
        self.lcpMaxTrace = []
        self.rcpMaxTrace = []
        self.lcpMaxCanvas.clear()
        self.rcpMaxCanvas.clear()
        
    def onImageUpdate(self, value):
        self.imageUpdate = value
        if (self.imageUpdate):
            self.buildImage()
        
    def onImageLimb(self, value):
        self.buildImage()
        
    def onImageCentering(self, centering):
        if (centering):
            self.figure3d = PL.figure()
            self.plot3d = self.figure3d.add_subplot(111, projection='3d')
            self.plot3d.set_zlim(0,2)
        else:
            PL.close()

    def onImageCalibrate(self, calibrate):
        self.buildImage()

    def onImageAnimate(self, value):
        self.imageAnimate = value
        if (self.imageAnimate):
            self.animTimer.start(100)
        else:
            self.animTimer.stop()

    def onPhaseCorrect(self, value):
        self.phaseCorrect = value
        if (self.imageUpdate):
            self.buildImage()

    def onAmplitudeCorrect(self, value):
        self.amplitudeCorrect = value
        if (self.imageUpdate):
            self.buildImage()
 
    def onFrequencyChannelChanged(self, value):
        self.currentFrequencyChannel = value
        self.srhFits.setFrequencyChannel(value)
        self.frequencyList.setCurrentIndex(self.currentFrequencyChannel)
        self.ewLcpPhaseSlopeSpin.setValue(self.ewLcpPhaseSlope[self.currentFrequencyChannel])
        self.ewRcpPhaseSlopeSpin.setValue(self.ewRcpPhaseSlope[self.currentFrequencyChannel])
        self.sLcpPhaseSlopeSpin.setValue(self.sLcpPhaseSlope[self.currentFrequencyChannel])
        self.sRcpPhaseSlopeSpin.setValue(self.sRcpPhaseSlope[self.currentFrequencyChannel])
        if (self.imageUpdate):
            self.buildSPhase()
            self.buildEwPhase()
            self.buildImage()

    def onScanChanged(self, value):
        self.currentScan = value
        self.srhFits.getHourAngle(self.currentScan)
        self.timeList.setCurrentIndex(self.currentScan)
        if (self.imageUpdate):
            self.buildImage()
            if self.imageCenteringButton.isChecked():
                X0 = 128
                Y0 = 128
                dN = 32
                XX = NP.linspace(X0-dN,X0+dN,2*dN)
                YY = NP.linspace(Y0-dN,Y0+dN,2*dN)
                X, Y = NP.meshgrid(XX,YY)
                F2O = {}
                F2O[2] = 2.0
                F2O[5] = 1.5
                F2O[8] = 1.0
                F2O[11] = 0.5
                F2O[14] = 0.0

                F2L = {}
#                F2L[2] =   [3.0,4.0,5.0,6.0]
#                F2L[5] =   [2.0,3.0,4.0,5.0]
#                F2L[8] =   [1.5,2.0,2.5,3.0]
#                F2L[11] = [1.0,1.5,2.0,2.5]
#                F2L[14] = [1.0,1.5,2.0,2.5]
                F2L[2] =  [20.0,30.0,80.0,100.0]
                F2L[5] =  [20.0,30.0,80.0,100.0]
                F2L[8] =  [20.0,30.0,80.0,100.0]
                F2L[11] = [20.0,30.0,80.0,100.0]
                F2L[14] = [20.0,30.0,80.0,100.0]

                self.plot3d.clear()
                self.plot3d.set_zlim(0,2)
                self.imageUpdate = False
                for freq in [2, 5, 8, 11, 14]:
                    self.currentFrequencyChannel = freq
                    self.srhFits.setFrequencyChannel(freq)
                    self.buildSPhase()
                    self.buildEwPhase()
                    self.buildImage()
                    lcpData02 = self.lcpData.copy()
                    max_ind = NP.unravel_index( NP.argmax(lcpData02), lcpData02.shape)
                    lcpData02 = NP.roll(lcpData02, lcpData02.shape[0] // 2 - max_ind[0], axis=0)
                    lcpData02 = NP.roll(lcpData02, lcpData02.shape[1] // 2 - max_ind[1], axis=1)
                    self.plot3d.contourf(X,Y,NP.transpose(lcpData02[Y0-dN:Y0+dN,X0-dN:X0+dN]), F2L[freq], offset=F2O[freq])
                self.figure3d.canvas.draw()
                PL.show()
                self.imageUpdate = True
                XX, X, YY, Y, lcpData02 = 0

    def onEwAntennaChanged(self, value):
        self.ewAntennaNumber = value
        self.ewLcpAntennaPhase.setValue(self.ewLcpPhaseAnt[self.ewAntennaNumber - 49])
        self.ewRcpAntennaPhase.setValue(self.ewRcpPhaseAnt[self.ewAntennaNumber - 49])
    
    def onSAntennaChanged(self, value):
        self.sAntennaNumber = value
        self.sLcpAntennaPhase.setValue(self.sLcpPhaseAnt[15 - (self.sAntennaNumber - 177)])
        self.sRcpAntennaPhase.setValue(self.sRcpPhaseAnt[15 - (self.sAntennaNumber - 177)])
    
    def onEwLcpAntennaPhaseChanged(self, value):
        self.ewLcpPhaseAnt[self.ewAntennaNumber - 49] = value
        self.buildEwPhase()
        if (self.imageUpdate):
            self.buildImage()
    
    def onSLcpAntennaPhaseChanged(self, value):
        self.sLcpPhaseAnt[15 - (self.sAntennaNumber - 177)] = value
        self.buildSPhase()
        if (self.imageUpdate):
            self.buildImage()
    
    def onEwRcpAntennaPhaseChanged(self, value):
        self.ewRcpPhaseAnt[self.ewAntennaNumber - 49] = value
        self.buildEwPhase()
        if (self.imageUpdate):
            self.buildImage()
    
    def onSRcpAntennaPhaseChanged(self, value):
        self.sRcpPhaseAnt[15 - (self.sAntennaNumber - 177)] = value
        self.buildSPhase()
        if (self.imageUpdate):
            self.buildImage()
    
    def onCalibScanChanged(self, value):
        self.srhFits.setCalibIndex(value)
        if (self.imageUpdate):
            self.buildImage()

    def onImageOffsetSlider(self, value):
        self.imageOffset = value*0.1
        if (self.imageUpdate):
            self.buildImage()

    def onImageScaleSlider(self, value):
        self.imageScale = value*0.1
        if (self.imageUpdate):
            self.buildImage()
        
    def onAnimTimer(self):
        if (self.currentScan < self.srhFits.dataLength):
            self.currentScan += 1
            self.scan.setValue(self.currentScan)
        else:
            self.animTimer.stop()
            self.imageAnimateButton.setChecked(False)
            self.scan.setValue(0)
        
    def onFrequencyListSelected(self):
        self.frequencyChannel.setValue(self.frequencyList.currentIndex())
        
    def onTimeListSelected(self):
        self.scan.setValue(self.timeList.currentIndex())
        
    def onCanvasXyChanged(self, x, y):
        self.xInd = int(x)
        self.yInd = int(y)
        if self.srhFits.isOpen:
            self.lcpMaxCanvas.clear()
            if self.indexTypeOfImage == 0:
                self.lcpMaxCanvas.plot(self.lcpData[self.yInd,:])
                self.lcpMaxCanvas.plot(self.lcpData[:,self.xInd])
            else:
                self.lcpMaxCanvas.plot((self.lcpData[self.yInd,:] + self.rcpData[self.yInd,:]*self.lcpRcpRel)*.5)
                self.lcpMaxCanvas.plot((self.lcpData[:,self.xInd] + self.rcpData[:,self.xInd]*self.lcpRcpRel)*.5)
        
            self.rcpMaxCanvas.clear()
            if self.indexTypeOfImage == 0:
                self.rcpMaxCanvas.plot(self.rcpData[self.yInd,:])
                self.rcpMaxCanvas.plot(self.rcpData[:,self.xInd])
            else:
                self.rcpMaxCanvas.plot((self.lcpData[self.yInd,:] - self.rcpData[self.yInd,:]*self.lcpRcpRel)*.5)
                self.rcpMaxCanvas.plot((self.lcpData[:,self.xInd] - self.rcpData[:,self.xInd]*self.lcpRcpRel)*.5)

    def onTypeOfImage(self, index):
        self.indexOfImageType = index
        if (self.imageUpdate):
            self.buildImage()
        
    def onTypeOfFrame(self, index):
        self.indexOfFrameType = index
        if (self.imageUpdate):
            self.buildImage()

    def onLcpRcpRelationChanged(self, value):
        self.lcpRcpRel = value
        if (self.imageUpdate):
            self.buildImage()
        
    def closeEvent(self, event):
        close = QtWidgets.QMessageBox()
        close.setText("Are you sure to exit?")
        close.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.Cancel)
        close = close.exec()

        if close == QtWidgets.QMessageBox.Yes:
            settings = QtCore.QSettings('xSrhEdik.conf', QtCore.QSettings.IniFormat)
            settings.setValue('currentFrequencyChannel', self.currentFrequencyChannel)
            settings.setValue('currentScan', self.currentScan)
            settings.setValue('imageOffset', self.imageOffset)
            settings.setValue('imageScale', self.imageScale)
            settings.setValue('ewLcpPhaseSlope', self.ewLcpPhaseSlope[self.currentFrequencyChannel])
            settings.setValue('sLcpPhaseSlope', self.sLcpPhaseSlope[self.currentFrequencyChannel])
            settings.setValue('phaseCorrect', self.phaseCorrect)
            settings.setValue('amplitudeCorrect', self.amplitudeCorrect)
            settings.setValue('ewStairLength', self.ewStairLength)
            settings.setValue('sStairLength', self.sStairLength)
            settings.setValue('imageUpdate', self.imageUpdate)
            settings.setValue('ewLcpPhaseSlope', self.ewLcpPhaseSlope[self.currentFrequencyChannel])
            settings.setValue('sLcpPhaseSlope', self.sLcpPhaseSlope[self.currentFrequencyChannel])

            event.accept()
        else:
            event.ignore()        
        
    def resizeEvent(self, event):
        width = event.size().width() // 2
        height = event.size().height()
        self.lcpCanvas.setGeometry(0, 50, width, width)
        self.rcpCanvas.setGeometry(width, 50, width, width)
        self.lcpMaxCanvas.setGeometry(0,width + 50,width, height - width - 50)
        self.rcpMaxCanvas.setGeometry(width,width + 50,width, height - width - 50)
        
    def antennasPhases(self):
        ewPhases = NP.zeros((self.srhFits.dataLength,2, 32))
        sPhases = NP.zeros((self.srhFits.dataLength,2, 16))
        for scan in range(self.srhFits.dataLength):
            self.srhFits.calibIndex = scan
            self.srhFits.updateAntennaPhase()
            ewPhases[scan, 0, :] = self.srhFits.ewAntPhaLcp[1:]
            ewPhases[scan, 1, :] = self.srhFits.ewAntPhaRcp[1:]
            sPhases[scan, 0] = self.srhFits.sAntPhaLcp[1:]
            sPhases[scan, 1] = self.srhFits.sAntPhaRcp[1:]
#            sPhases[scan, 0, 0:6] = self.srhFits.sAntPhaLcp6[1:]
#            sPhases[scan, 1, 0:6] = self.srhFits.sAntPhaRcp6[1:]
#            sPhases[scan, 0, 11]  = self.srhFits.sAntPhaLcp5[1]
#            sPhases[scan, 1, 11]  = self.srhFits.sAntPhaRcp5[1]
#            sPhases[scan, 0, 12:] = self.srhFits.sAntPhaLcp5[2:]
#            sPhases[scan, 1, 12:] = self.srhFits.sAntPhaRcp5[2:]
        
        return ewPhases, sPhases

    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self,parent)
        
#        settings = QtCore.QSettings('xSrhEdik.conf', QtCore.QSettings.IniFormat)
#        
#        self.currentFrequencyChannel = settings.value('currentFrequencyChannel')
#        self.currentScan = settings.value('currentScan')
#        self.imageOffset = settings.value('imageOffset')
#        self.imageScale = settings.value('imageScale')
#        self.ewLcpPhaseSlope = settings.value('ewLcpPhaseSlope')
#        self.sLcpPhaseSlope = settings.value('sLcpPhaseSlope')
#        self.phaseCorrect = settings.value('phaseCorrect')
#        self.amplitudeCorrect = settings.value('amplitudeCorrect')
#        self.ewStairLength = settings.value('ewStairLength')
#        self.sStairLength = settings.value('sStairLength')
#        self.imageUpdate = settings.value('imageUpdate')

        self.currentFrequencyChannel = 5
        self.currentScan = 0
        self.imageOffset = 0
        self.imageScale = 1
        self.phaseCorrect = True
        self.amplitudeCorrect = False
        self.indexOfFrameType = 0
        self.pAngle = 0.
        self.ewStairLength = 16
        self.sStairLength = 16
        self.imageUpdate = True
        self.lcpRcpRel = 1.
        self.arcsecPerPixel = 4.91104

        self.lcpMaxTrace = []
        self.rcpMaxTrace = []
        self.ewLcpPhaseCorrection = NP.zeros(32)
        self.ewRcpPhaseCorrection = NP.zeros(32)
        self.sLcpPhaseCorrection = NP.zeros(16)
        self.sRcpPhaseCorrection = NP.zeros(16)
        self.ewLcpPhaseAnt = NP.zeros(32, dtype=int)
        self.ewRcpPhaseAnt = NP.zeros(32, dtype=int)
        self.sLcpPhaseAnt = NP.zeros(16, dtype=int)
        self.sRcpPhaseAnt = NP.zeros(16, dtype=int)
        self.ewAntennaNumber = 49
        self.sAntennaNumber = 177
        self.imageUpdate = True
        self.imageAnimate = False
        self.animTimer = QtCore.QTimer(self)
        self.animTimer.timeout.connect(self.onAnimTimer)
        
        self.setGeometry(100,100,1024,700)
        self.openButton = QtWidgets.QPushButton('Open...', self)
        self.openButton.clicked.connect(self.onOpen)

        self.saveButton = QtWidgets.QPushButton('Save as...', self)
        self.saveButton.clicked.connect(self.onSaveAs)

        self.findPhaseButton = QtWidgets.QPushButton('Find phase', self)
        self.findPhaseButton.clicked.connect(self.onFindPhase)

        self.ewPhaseStairLcp = QtWidgets.QSpinBox(self, prefix='EWL ')
        self.ewPhaseStairLcp.setRange(-180,180)
        self.ewPhaseStairLcp.valueChanged.connect(self.onEastWestPhaseStairLcpChanged)
        self.ewPhaseStairRcp = QtWidgets.QSpinBox(self, prefix='EWR ')
        self.ewPhaseStairRcp.setRange(-180,180)
        self.ewPhaseStairRcp.valueChanged.connect(self.onEastWestPhaseStairRcpChanged)

        self.ewPhaseStairLength = QtWidgets.QSpinBox(self, prefix='period ')
        self.ewPhaseStairLength.setRange(1,16)
        self.ewPhaseStairLength.setValue(16)
        self.ewPhaseStairLength.valueChanged.connect(self.onEwPhaseStairLengthChanged)

        self.sPhaseStairLcp = QtWidgets.QSpinBox(self, prefix='SL ')
        self.sPhaseStairLcp.setRange(-180,180)
        self.sPhaseStairLcp.valueChanged.connect(self.onSouthPhaseStairLcpChanged)
        self.sPhaseStairRcp = QtWidgets.QSpinBox(self, prefix='SR ')
        self.sPhaseStairRcp.setRange(-180,180)
        self.sPhaseStairRcp.valueChanged.connect(self.onSouthPhaseStairRcpChanged)

        self.sPhaseStairLength = QtWidgets.QSpinBox(self, prefix = 'period ')
        self.sPhaseStairLength.setRange(1,16)
        self.sPhaseStairLength.setValue(16)
        self.sPhaseStairLength.valueChanged.connect(self.onSPhaseStairLengthChanged)

        self.lcpCanvas = ResponseCanvas(self)
        self.lcpCanvas.mouseSignal.connect(self.onCanvasXyChanged)
        self.lcpMaxCanvas = ResponseCanvas(self)
        self.rcpCanvas = ResponseCanvas(self)
        self.rcpCanvas.mouseSignal.connect(self.onCanvasXyChanged)
        self.rcpMaxCanvas = ResponseCanvas(self)

        self.clearButton = QtWidgets.QPushButton('Clear trace', self)
        self.clearButton.clicked.connect(self.onClear)
        
        self.imageUpdateButton = QtWidgets.QPushButton('Update', self)
        self.imageUpdateButton.setCheckable(True)
        self.imageUpdateButton.setChecked(True)
        self.imageUpdateButton.clicked.connect(self.onImageUpdate)
        
        self.imageLimbButton = QtWidgets.QPushButton('Limb', self)
        self.imageLimbButton.setCheckable(True)
        self.imageLimbButton.clicked.connect(self.onImageLimb)
        
        self.findMinLcpRcpButton = QtWidgets.QPushButton('Find min', self)
        self.findMinLcpRcpButton.clicked.connect(self.onFindMinLcpRcp)

        self.imageCenteringButton = QtWidgets.QPushButton('Centering', self)
        self.imageCenteringButton.setCheckable(True)
        self.imageCenteringButton.setChecked(False)
        self.imageCenteringButton.clicked.connect(self.onImageCentering)
        
        self.imageCalibrateButton = QtWidgets.QPushButton('Calibrate', self)
        self.imageCalibrateButton.setCheckable(True)
        self.imageCalibrateButton.setChecked(False)
        self.imageCalibrateButton.clicked.connect(self.onImageCalibrate)

        self.frequencyChannel = QtWidgets.QSpinBox(self, prefix='channel ')
        self.frequencyChannel.setRange(0,0)
        self.frequencyChannel.valueChanged.connect(self.onFrequencyChannelChanged)

        self.frequencyList = QtWidgets.QComboBox(self)
        self.frequencyList.currentIndexChanged.connect(self.onFrequencyListSelected)
        self.timeList = QtWidgets.QComboBox(self)
        self.timeList.currentIndexChanged.connect(self.onTimeListSelected)
        
        self.scan = QtWidgets.QSpinBox(self, prefix='scan ')
        self.scan.setRange(0,0)
        self.scan.valueChanged.connect(self.onScanChanged)

        self.calibScan = QtWidgets.QSpinBox(self, prefix='c_scan ')
        self.calibScan.setRange(0,0)
        self.calibScan.valueChanged.connect(self.onCalibScanChanged)

        self.ewLcpPhaseSlopeSpin = QtWidgets.QSpinBox(self, prefix = 'EWLS ')
        self.ewLcpPhaseSlopeSpin.setRange(-180.,180.)
        self.ewLcpPhaseSlopeSpin.valueChanged.connect(self.onEastWestLcpPhaseSlopeChanged)

        self.ewRcpPhaseSlopeSpin = QtWidgets.QSpinBox(self, prefix='EWRS ')
        self.ewRcpPhaseSlopeSpin.setRange(-180.,180.)
        self.ewRcpPhaseSlopeSpin.valueChanged.connect(self.onEastWestRcpPhaseSlopeChanged)

        self.sLcpPhaseSlopeSpin = QtWidgets.QSpinBox(self, prefix='SLS ')
        self.sLcpPhaseSlopeSpin.setRange(-180.,180.)
        self.sLcpPhaseSlopeSpin.valueChanged.connect(self.onSouthLcpPhaseSlopeChanged)

        self.sRcpPhaseSlopeSpin = QtWidgets.QSpinBox(self, prefix='SRS ')
        self.sRcpPhaseSlopeSpin.setRange(-180.,180.)
        self.sRcpPhaseSlopeSpin.valueChanged.connect(self.onSouthRcpPhaseSlopeChanged)

        self.phaseCorrectButton = QtWidgets.QPushButton('Phase', self)
        self.phaseCorrectButton.setCheckable(True)
        self.phaseCorrectButton.setChecked(True)
        self.phaseCorrectButton.clicked.connect(self.onPhaseCorrect)
        
        self.amplitudeCorrectButton = QtWidgets.QPushButton('Amplitude', self)
        self.amplitudeCorrectButton.setCheckable(True)
        self.amplitudeCorrectButton.clicked.connect(self.onAmplitudeCorrect)
        
        self.typeOfFrame = QtWidgets.QComboBox(self)
        self.typeOfFrame.currentIndexChanged.connect(self.onTypeOfFrame)
        self.typeOfFrame.addItem('l, m')
        self.typeOfFrame.addItem('h,d')
        self.typeOfFrame.addItem('h-P,d-P')

        self.imageAnimateButton = QtWidgets.QPushButton('Animate', self)
        self.imageAnimateButton.setCheckable(True)
        self.imageAnimateButton.setChecked(self.imageAnimate)
        self.imageAnimateButton.clicked.connect(self.onImageAnimate)
        
        self.imageOffsetSlider = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.imageOffsetSlider.setRange(-100,100)
        self.imageOffsetSlider.setValue(0)
        self.imageOffsetSlider.valueChanged.connect(self.onImageOffsetSlider)
        
        self.imageScaleSlider = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.imageScaleSlider.setRange(1,100)
        self.imageScaleSlider.setValue(10)
        self.imageScaleSlider.valueChanged.connect(self.onImageScaleSlider)
        
        self.typeOfImage = QtWidgets.QComboBox(self)
        self.typeOfImage.currentIndexChanged.connect(self.onTypeOfImage)
        self.typeOfImage.addItem('LCP, RCP')
        self.typeOfImage.addItem('I,V')
        self.typeOfImage.addItem('uvLCP,uvRCP')
        self.typeOfImage.addItem('PSF')
        self.typeOfImage.addItem('IV + PSF')
        
        self.lcpRcpRelation = QtWidgets.QDoubleSpinBox (self)
        self.lcpRcpRelation.setSingleStep(0.01)
        self.lcpRcpRelation.setRange(0.,2.)
        self.lcpRcpRelation.valueChanged.connect(self.onLcpRcpRelationChanged)
        self.lcpRcpRelation.setValue(1.)
        
        self.ewAntenna = QtWidgets.QSpinBox(self, prefix='EW_ant ')
        self.ewAntenna.setRange(49,80)
        self.ewAntenna.valueChanged.connect(self.onEwAntennaChanged)

        self.ewLcpAntennaPhase = QtWidgets.QSpinBox(self)
        self.ewLcpAntennaPhase.setRange(-180,180)
        self.ewLcpAntennaPhase.valueChanged.connect(self.onEwLcpAntennaPhaseChanged)

        self.sLcpAntennaPhase = QtWidgets.QSpinBox(self)
        self.sLcpAntennaPhase.setRange(-180,180)
        self.sLcpAntennaPhase.valueChanged.connect(self.onSLcpAntennaPhaseChanged)

        self.ewRcpAntennaPhase = QtWidgets.QSpinBox(self)
        self.ewRcpAntennaPhase.setRange(-180,180)
        self.ewRcpAntennaPhase.valueChanged.connect(self.onEwRcpAntennaPhaseChanged)

        self.sRcpAntennaPhase = QtWidgets.QSpinBox(self)
        self.sRcpAntennaPhase.setRange(-180,180)
        self.sRcpAntennaPhase.valueChanged.connect(self.onSRcpAntennaPhaseChanged)

        self.sAntenna = QtWidgets.QSpinBox(self, prefix='S_ant ')
        self.sAntenna.setRange(177,192)
        self.sAntenna.valueChanged.connect(self.onSAntennaChanged)

        self.openButton.setGeometry(0,0,60,25)
        self.frequencyChannel.setGeometry(70,0,80,25)
        self.scan.setGeometry(70,25,80,25)
        self.calibScan.setGeometry(150,0,80,25)
        self.clearButton.setGeometry(150,25,80,25)
        self.ewPhaseStairLcp.setGeometry(235,0,70,25)
        self.ewPhaseStairRcp.setGeometry(305,0,70,25)
        self.ewPhaseStairLength.setGeometry(235,25,140,25)

        self.sPhaseStairLcp.setGeometry(375,0,70,25)
        self.sPhaseStairRcp.setGeometry(445,0,70,25)
        self.sPhaseStairLength.setGeometry(375,25,140,25)

        self.ewLcpPhaseSlopeSpin.setGeometry(500,0,80,25)
        self.sLcpPhaseSlopeSpin.setGeometry(500,25,80,25)
        self.ewRcpPhaseSlopeSpin.setGeometry(580,0,80,25)
        self.sRcpPhaseSlopeSpin.setGeometry(580,25,80,25)

        self.phaseCorrectButton.setGeometry(660,0,60,25)
        self.amplitudeCorrectButton.setGeometry(660,25,60,25)
        
        self.typeOfFrame.setGeometry(720,0,60,25)
        
        self.imageAnimateButton.setGeometry(720,25,60,25)

        self.imageOffsetSlider.setGeometry(780,0,60,25)
        self.imageScaleSlider.setGeometry(780,25,60,25)
        self.findPhaseButton.setGeometry(840,0,60,25)
        self.imageUpdateButton.setGeometry(840,25,60,25)
        self.frequencyList.setGeometry(900, 0, 100, 25)
        self.timeList.setGeometry(900, 25, 100, 25)
        
        self.typeOfImage.setGeometry(1000, 0, 100, 25)
        self.lcpRcpRelation.setGeometry(1000, 25, 100, 25)

        self.imageCenteringButton.setGeometry(1100,0,60,25)
        self.imageCalibrateButton.setGeometry(1100,25,60,25)

        self.ewAntenna.setGeometry(1160, 0, 80, 25)
        self.sAntenna.setGeometry(1160, 25, 80, 25)
        self.ewLcpAntennaPhase.setGeometry(1240, 0, 50, 25)
        self.sLcpAntennaPhase.setGeometry(1240, 25, 50, 25)
        self.ewRcpAntennaPhase.setGeometry(1290, 0, 50, 25)
        self.sRcpAntennaPhase.setGeometry(1290, 25, 50, 25)
        
        self.findMinLcpRcpButton.setGeometry(1340, 0, 50, 25)
        self.imageLimbButton.setGeometry(1340, 25, 50, 25)

        self.saveButton.setGeometry(0,25,60,25)
        
        self.lcpCanvas.setGeometry(0,50,512,512)
        self.lcpMaxCanvas.setGeometry(0,560,512,100)
        self.rcpCanvas.setGeometry(512,50,512,512)
        self.rcpMaxCanvas.setGeometry(512,560,512,100)
        
        
#        self.lcpCanvas.setGeometry(0,50,1000,1000)
#        self.lcpMaxCanvas.setGeometry(0,1060,512,100)
#        self.rcpCanvas.setGeometry(512,1050,512,512)
#        self.rcpMaxCanvas.setGeometry(512,1560,512,100)

    def onOpen(self):
        fitsNames, _ = QtWidgets.QFileDialog.getOpenFileNames(self)        
        if fitsNames[0]:
            self.uvSize = 512
            self.currentScan = 0
            self.scan.setValue(0)
            self.srhFits = SrhFitsFile(fitsNames[0], self.uvSize)
            self.srhFits.getHourAngle(self.currentScan)
            self.pAngle = NP.deg2rad(coordinates.get_sun_P(self.srhFits.dateObs).to_value())
            self.srhFits.setCalibIndex(0)
            self.srhFits.vis2uv(self.currentScan, phaseCorrect=True, PSF=False);
            self.srhFits.uv2lmImage()

            self.frequencyList.clear()
            for freq in self.srhFits.freqList:
                self.frequencyList.addItem(str(freq))

            self.frequencyChannel.setRange(0, self.srhFits.freqListLength)
            self.frequencyChannel.setValue(self.currentFrequencyChannel)
            self.scan.setRange(0, self.srhFits.dataLength)
            self.scan.setValue(self.currentScan)
            self.calibScan.setRange(0, self.srhFits.dataLength)

            self.ewPhaseCoefsLcp = NP.zeros((self.srhFits.freqListLength, 16))
            self.ewPhaseCoefsRcp = NP.zeros((self.srhFits.freqListLength, 16))
            self.sPhaseCoefsLcp = NP.zeros((self.srhFits.freqListLength, 16))
            self.sPhaseCoefsRcp = NP.zeros((self.srhFits.freqListLength, 16))
            
            self.ewLcpPhaseSlope = NP.zeros((self.srhFits.freqListLength))
            self.ewRcpPhaseSlope = NP.zeros((self.srhFits.freqListLength))
            self.sLcpPhaseSlope = NP.zeros((self.srhFits.freqListLength))
            self.sRcpPhaseSlope = NP.zeros((self.srhFits.freqListLength))
            
            self.qSun = NP.zeros((self.uvSize, self.uvSize))
            sunRadius = 16 * 60 / self.arcsecPerPixel
            for i in range(self.uvSize):
                x = i - self.uvSize/2
                for j in range(self.uvSize):
                    y = j - self.uvSize/2
                    if (NP.sqrt(x*x + y*y) < sunRadius):
                        self.qSun[i , j] = 1.
                        
            data = NP.flip(self.srhFits.lcp.real,0)
            self.lcpCanvas.imshow(data, NP.min(data), NP.max(data))
            self.lcpCanvas.plot([0,1])
            self.lcpMaxTrace.append(NP.max(data))
#            self.lcpMaxCanvas.plot(self.lcpMaxTrace)

            data = NP.flip(self.srhFits.rcp.real,0)
            self.rcpCanvas.imshow(data, 1.*NP.min(data), 1.*NP.max(data))
            self.rcpMaxTrace.append(NP.max(data))
#            self.rcpMaxCanvas.plot(self.rcpMaxTrace)

            self.setWindowTitle('SRH phase edit:' + fitsNames[0])

            for fitsName in fitsNames[1:]:
                self.srhFits.append(fitsName)
                self.scan.setRange(0, self.srhFits.dataLength)
                self.calibScan.setRange(0, self.srhFits.dataLength)

        for tim in self.srhFits.freqTime[0,:]:
            fTime = QtCore.QTime(0,0)
            fTime = fTime.addMSecs(tim * 1000)
            self.timeList.addItem(fTime.toString('hh:mm:ss'))
                
    def saveAsFits(self, saveName):
        for scan in range(self.srhFits.dataLength):
            self.scan.setValue(scan)
            nameParts = saveName.split('.') 
            fitsName = ('%s_%03d.%s' % (nameParts[0], scan, nameParts[1]))
            print(fitsName)
            if self.indexOfImageType == 3:
                self.lcpData /= NP.max(self.lcpData)
                self.rcpData /= NP.max(self.rcpData)
    
            pHeader = fits.Header();
            
            pHeader['DATE-OBS']     = self.srhFits.hduList[0].header['DATE-OBS']
            pHeader['T-OBS']     = self.srhFits.hduList[0].header['DATE-OBS'] + 'T' + phaseEdit.timeList.itemText(scan)#self.srhFits.hduList[0].header['TIME-OBS']
            pHeader['INSTRUME']     = self.srhFits.hduList[0].header['INSTRUME']
            pHeader['ORIGIN']       = self.srhFits.hduList[0].header['ORIGIN']
            pHeader['OBS-LAT']      = self.srhFits.hduList[0].header['OBS-LAT']
            pHeader['OBS-LONG']     = self.srhFits.hduList[0].header['OBS-LONG']
            pHeader['OBS-ALT']      = self.srhFits.hduList[0].header['OBS-ALT']
            pHeader['FR_CHAN']      = self.srhFits.hduList[0].header['FR_CHAN']
            pHeader['FREQUENC']     = ('%3.3f') % (self.srhFits.freqList[self.currentFrequencyChannel]*1e6)
            pHeader['CDELT1']       = self.arcsecPerPixel
            pHeader['CDELT2']       = self.arcsecPerPixel
            pHeader['CRPIX1']       = self.uvSize // 2
            pHeader['CRPIX2']       = self.uvSize // 2
            pHeader['CTYPE1']       = 'HPLN-TAN'
            pHeader['CTYPE2']       = 'HPLT-TAN'
            pHeader['CUNIT1']       = 'arcsec'
            pHeader['CUNIT2']       = 'arcsec'
            

            savePath, saveExt = fitsName.split('.')
            
            pHdu = fits.PrimaryHDU(self.lcpData + self.rcpData, header=pHeader);
            hduList = fits.HDUList([pHdu]);
            hduList.writeto(savePath + '_I.' + saveExt, clobber=True);
            hduList.close();
            
            pHdu = fits.PrimaryHDU(self.lcpData - self.rcpData, header=pHeader);
            hduList = fits.HDUList([pHdu]);
            hduList.writeto(savePath + '_V.' + saveExt, clobber=True);
            
            hduList.close();
        
    def updateCanvas(self, scan):
        self.scan.setValue(scan)
        self.buildImage()
        return self.rcpCanvas.imageObject,
        
    def saveAsMp4(self, saveName):
        ani = animation.FuncAnimation(self.rcpCanvas.fig, self.updateCanvas, frames=self.srhFits.dataLength, blit=True, repeat=False)
        ani.save(saveName)
        

    def saveAsMs2(self, saveName):
        ms2Table = srhMS2.SrhMs2(saveName)
        ms2Table.initDataTable(self.srhFits, self.currentFrequencyChannel, self.ewLcpPhaseCorrection, self.ewRcpPhaseCorrection, self.sLcpPhaseCorrection, self.sRcpPhaseCorrection)
        ms2Table.initAntennaTable(self.srhFits)
        ms2Table.initSpectralWindowTable(self.srhFits, self.currentFrequencyChannel)
        ms2Table.initDataDescriptionTable()
        ms2Table.initPolarizationTable()
        ms2Table.initSourceTable()
        ms2Table.initFieldTable()
        ms2Table.initFeedTable(self.srhFits, self.currentFrequencyChannel)
        ms2Table.initObservationTable(self.srhFits)
        
    def onSaveAs(self):
        saveName, _ = QtWidgets.QFileDialog.getSaveFileName(self)        
        if saveName:
            saveExt = saveName.split('.')[1]
            if saveExt == 'fit' or saveExt == 'fits':
                self.saveAsFits(saveName)
            elif saveExt == 'ms':
                self.saveAsMs2(saveName)
            elif saveExt == 'mp4':
                self.saveAsMp4(saveName)
                               
application = QtWidgets.QApplication.instance();
if not application:
    application = QtWidgets.QApplication(sys.argv);
    
if sys.platform == 'linux':
    font = QtGui.QFont()
    application.setFont(QtGui.QFont(font.defaultFamily(),8));

phaseEdit = SrhEdik();
phaseEdit.setWindowTitle('SRH editor')
phaseEdit.show();
sys.exit(application.exec_());
