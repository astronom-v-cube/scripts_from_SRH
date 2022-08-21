#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:21:36 2020

@author: maria
"""

import sys;
from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtWidgets import QFileDialog, QTabWidget, QVBoxLayout, QGridLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as NP
from astropy.io import fits
import struct

class ResponseCanvas(FigureCanvas):
    mouseSignal = QtCore.pyqtSignal(float, float, name = 'xyChanged')
    
    def mouse_moved(self, event):
        1
#        self.mouseSignal.emit(event.xdata, event.ydata)
        
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.subplot = self.fig.add_subplot(111)
        self.subplot.xaxis.set_visible(True)
        self.subplot.yaxis.set_visible(True)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.mpl_connect('motion_notify_event', self.mouse_moved)
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
    
    def plotXY(self, x, y):
        self.plotObject = self.subplot.plot(x, y)
        self.draw()
        
    def setXLabel(self, label):
        self.subplot.set_xlabel(label)
        self.draw()
        
    def setYLabel(self, label):
        self.subplot.set_ylabel(label)
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

class FitsView(QtWidgets.QWidget):#MainWindow):

    def visIndex(self, hh, vv):
        h = hh if hh>vv else vv
        v = vv if hh>vv else hh
        ind = (h//4 - v//4) * 16 + (2 * 52 - v//4 + 1) * (v//4) * 16 // 2 + v%4 * 4 + h%4
        return ind

    def antennaIndex(self, ant):
        ind = []
        if ant >= 81 and ant <=128:
            ind = NP.where(self.receiver1 == ant)[0][0]
        elif ant >= 1 and ant <=48:
            ind = NP.where(self.receiver3 == ant)[0][0] + 32
        elif ant >= 129 and ant <=176:
            ind = NP.where(self.receiver5 == ant)[0][0] + 64
        elif ant >= 197 and ant <=227:
            ind = NP.where(self.receiver7 == ant)[0][0] + 96
        return ind
    
    def buildImageLeft(self):
        if not self.leftOverplotButton.isChecked():
            self.leftCanvas.clear()
        if self.leftPlotX.currentIndex() == 0: # TIME
            xLabel = 'Time'
#            x = self.time[0]
            x = NP.arange(0, len(self.time[0]))
            if self.leftPolarization.currentText == 'LCP':
                y = self.visibilitiesLcp[:, self.visIndex(self.antennaIndex(self.antPairs[self.l_pair,0]), self.antennaIndex(self.antPairs[self.l_pair,1]))]
                amp = self.ampLcp
            else:
                y = self.visibilitiesRcp[:, self.visIndex(self.antennaIndex(self.antPairs[self.l_pair,0]), self.antennaIndex(self.antPairs[self.l_pair,1]))]
                amp = self.ampRcp
        elif self.leftPlotX.currentIndex() == 1: # ANTENNA
            xLabel = 'Antenna'
            x = NP.arange(0,64)
            ind = NP.array([NP.arange(0,16), NP.arange(32,48), NP.arange(64,80), NP.arange(96,112)]).reshape(64)
            if self.leftPolarization.currentText == 'LCP':
                y = self.ampLcp[self.l_scan, ind]
            else:
                y = self.ampRcp[self.l_scan, ind]
        elif self.leftPlotX.currentIndex() == 2: # BASELINE
            xLabel = 'Baseline'
            x = NP.arange(0,60)
            if self.leftPolarization.currentText == 'LCP':
                y = self.visibilitiesLcp[self.l_scan, self.redundantVisIndexes.astype(int)]
            else:
                y = self.visibilitiesRcp[self.l_scan, self.redundantVisIndexes.astype(int)]

            
        if self.leftPlotY.currentIndex() == 0: # ABS
            yLabel = 'Abs'
            y = NP.abs(y)
        elif self.leftPlotY.currentIndex() == 1: # PHASE
            yLabel = 'Phase'
            y = NP.angle(y)
        elif self.leftPlotY.currentIndex() == 2: # REAL
            yLabel = 'Real'
            y = y.real
        elif self.leftPlotY.currentIndex() == 3: # IMAG
            yLabel = 'Imag'
            y = y.imag
        elif self.leftPlotY.currentIndex() == 4: # ANT AMP
            yLabel = 'Amplitude'
            y = amp[:, self.l_antenna]
        elif self.leftPlotY.currentIndex() == 5: # MEAN ANT AMP
            yLabel = 'Amplitude'
            y = NP.sum(amp, axis = 1)/50.
            
        self.leftCanvas.plotXY(x, y)
        self.leftCanvas.setXLabel(xLabel)
        self.leftCanvas.setYLabel(yLabel)
        
    def buildImageRight(self):
        if not self.rightOverplotButton.isChecked():
            self.rightCanvas.clear()
        if self.rightPlotX.currentIndex() == 0: # TIME
            xLabel = 'Time'
#            x = self.time[0]
            x = NP.arange(0, len(self.time[0]))
            if self.rightPolarization.currentText == 'LCP':
                y = self.visibilitiesLcp[:, self.visIndex(self.antennaIndex(self.antPairs[self.r_pair,0]), self.antennaIndex(self.antPairs[self.r_pair,1]))]
                amp = self.ampLcp
            else:
                y = self.visibilitiesRcp[:, self.visIndex(self.antennaIndex(self.antPairs[self.r_pair,0]), self.antennaIndex(self.antPairs[self.r_pair,1]))]
                amp = self.ampRcp
        elif self.rightPlotX.currentIndex() == 1: # ANTENNA
            xLabel = 'Antenna'
            x = NP.arange(0,64)
            ind = NP.array([NP.arange(0,16), NP.arange(32,48), NP.arange(64,80), NP.arange(96,112)]).reshape(64)
            if self.rightPolarization.currentText == 'LCP':
                y = self.ampLcp[self.r_scan, ind]
            else:
                y = self.ampRcp[self.r_scan, ind]
        elif self.rightPlotX.currentIndex() == 2: # BASELINE
            xLabel = 'Baseline'
            x = NP.arange(0,60)
            if self.rightPolarization.currentText == 'LCP':
                y = self.visibilitiesLcp[self.r_scan, self.redundantVisIndexes.astype(int)]
            else:
                y = self.visibilitiesRcp[self.r_scan, self.redundantVisIndexes.astype(int)]

            
        if self.rightPlotY.currentIndex() == 0: # ABS
            yLabel = 'Abs'
            y = NP.abs(y)
        elif self.rightPlotY.currentIndex() == 1: # PHASE
            yLabel = 'Phase'
            y = NP.angle(y)
        elif self.rightPlotY.currentIndex() == 2: # REAL
            yLabel = 'Real'
            y = y.real
        elif self.rightPlotY.currentIndex() == 3: # IMAG
            yLabel = 'Imag'
            y = y.imag
        elif self.rightPlotY.currentIndex() == 4: # ANT AMP
            yLabel = 'Amplitude'
            y = amp[:, self.r_antenna]
        elif self.rightPlotY.currentIndex() == 5: # MEAN ANT AMP
            yLabel = 'Amplitude'
            y = NP.sum(amp, axis = 1)/50.
            
        self.rightCanvas.plotXY(x, y)
        self.rightCanvas.setXLabel(xLabel)
        self.rightCanvas.setYLabel(yLabel)
        
    def onImageUpdate(self, value):
        self.imageUpdate = value
        if (self.imageUpdate):
            self.buildImageLeft()
            self.buildImageRight()

#    def onPhaseCorrect(self, value):
#        self.phaseCorrect = value
#        if (self.imageUpdate):
#            self.buildImage()
#
#    def onAmplitudeCorrect(self, value):
#        self.amplitudeCorrect = value
#        if (self.imageUpdate):
#            self.buildImage()
#            
 
#    def onFrequencyChannelChanged(self, value):
#        self.currentFrequencyChannel = value
#        self.srhFits.setFrequencyChannel(value)
#        self.frequencyList.setCurrentIndex(self.currentFrequencyChannel)
#        self.ewLcpPhaseSlopeSpin.setValue(self.ewLcpPhaseSlope[self.currentFrequencyChannel])
#        self.ewRcpPhaseSlopeSpin.setValue(self.ewRcpPhaseSlope[self.currentFrequencyChannel])
#        self.sLcpPhaseSlopeSpin.setValue(self.sLcpPhaseSlope[self.currentFrequencyChannel])
#        self.sRcpPhaseSlopeSpin.setValue(self.sRcpPhaseSlope[self.currentFrequencyChannel])
#        if (self.imageUpdate):
#            self.buildSPhase()
#            self.buildEwPhase()
#            self.buildImage()

#    
#    def onCalibScanChanged(self, value):
#        self.srhFits.setCalibIndex(value)
#        if (self.imageUpdate):
#            self.buildImage()
#
#    def onImageOffsetSlider(self, value):
#        self.imageOffset = value*0.1
#        if (self.imageUpdate):
#            self.buildImage()
#
#    def onImageScaleSlider(self, value):
#        self.imageScale = value*0.1
#        if (self.imageUpdate):
#            self.buildImage()

#    def onFrequencyListSelected(self):
#        self.frequencyChannel.setValue(self.frequencyList.currentIndex())
        
#    def onTimeListSelected(self):
#        self.scan.setValue(self.timeList.currentIndex())
        
    def onLeftPlotXSelected(self):
        if (self.imageUpdate):
            self.buildImageLeft()
        
    def onLeftPlotYSelected(self):
        if (self.imageUpdate):
            self.buildImageLeft()
        
    def onRightPlotXSelected(self):
        if (self.imageUpdate):
            self.buildImageRight()
        
    def onRightPlotYSelected(self):
        if (self.imageUpdate):
            self.buildImageRight()
            
    def onLeftOverplot(self):
        if (self.imageUpdate):
            self.buildImageLeft()
    
    def onRightOverplot(self):
        if (self.imageUpdate):
            self.buildImageRight()
    
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

    def onLeftPolarizationSelected(self):
        if (self.imageUpdate):
            self.buildImageLeft()
    
    def onRightPolarizationSelected(self):
        if (self.imageUpdate):
            self.buildImageRight()
            
    def onLeftAntChanged(self, value):
        self.leftAnt.setPrefix('Antenna \''+str(self.antNumbers[value])+'\' ')
        if value >= 16 and value < 32: 
            self.l_antenna = value + 16
        elif value >= 32 and value < 48: 
            self.l_antenna = value + 32
        elif value >= 48 and value < 64: 
            self.l_antenna = value + 48
        else:
            self.l_antenna = value
        if (self.imageUpdate):
            self.buildImageLeft()
            
    def onLeftBaselineChanged(self, value):
        self.leftBaseline.setPrefix('Baseline \''+str(self.antPairs[value,0])+'-'+str(self.antPairs[value,1])+'\' ')
        self.l_pair = value
        if (self.imageUpdate):
            self.buildImageLeft()
        
    def onLeftScanChanged(self, value):
        self.l_scan = value
        if (self.imageUpdate):
            self.buildImageLeft()
            
    def onRightAntChanged(self, value):
        self.rightAnt.setPrefix('Antenna \''+str(self.antNumbers[value])+'\' ')
        if value >= 16 and value < 32: 
            self.r_antenna = value + 16
        elif value >= 32 and value < 48: 
            self.r_antenna = value + 32
        elif value >= 48 and value < 64: 
            self.r_antenna = value + 48
        else:
            self.r_antenna = value
        if (self.imageUpdate):
            self.buildImageRight()
            
    def onRightBaselineChanged(self, value):
        self.rightBaseline.setPrefix('Baseline \''+str(self.antPairs[value,0])+'-'+str(self.antPairs[value,1])+'\' ')
        self.r_pair = value
        if (self.imageUpdate):
            self.buildImageRight()
            
    def onRightScanChanged(self, value):
        self.r_scan = value
        if (self.imageUpdate):
            self.buildImageRight()
        
    def closeEvent(self, event):
        close = QtWidgets.QMessageBox()
        close.setText("Are you sure to exit?")
        close.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.Cancel)
        close = close.exec()

        if close == QtWidgets.QMessageBox.Yes:
#            settings = QtCore.QSettings('xSrhEdik.conf', QtCore.QSettings.IniFormat)
#            settings.setValue('currentFrequencyChannel', self.currentFrequencyChannel)
#            settings.setValue('currentScan', self.currentScan)
#            settings.setValue('imageOffset', self.imageOffset)
#            settings.setValue('imageScale', self.imageScale)
#            settings.setValue('ewLcpPhaseSlope', self.ewLcpPhaseSlope[self.currentFrequencyChannel])
#            settings.setValue('sLcpPhaseSlope', self.sLcpPhaseSlope[self.currentFrequencyChannel])
#            settings.setValue('phaseCorrect', self.phaseCorrect)
#            settings.setValue('amplitudeCorrect', self.amplitudeCorrect)
#            settings.setValue('ewStairLength', self.ewStairLength)
#            settings.setValue('sStairLength', self.sStairLength)
#            settings.setValue('imageUpdate', self.imageUpdate)
#            settings.setValue('ewLcpPhaseSlope', self.ewLcpPhaseSlope[self.currentFrequencyChannel])
#            settings.setValue('sLcpPhaseSlope', self.sLcpPhaseSlope[self.currentFrequencyChannel])

            event.accept()
        else:
            event.ignore()        
        

    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self,parent)
        self.receiver1 = NP.array([81,83,85,87,89,91,95,98,101,104,107,112,116,120,124,128])
        self.receiver3 = NP.array([48,46,44,42,40,38,34,31,28,25,22,17,13,9,5,1])
        self.receiver5 = NP.array([176,174,172,170,168,166,162,159,156,153,150,145,141,137,133,129])
        self.receiver7 = NP.array([197,199,201,203,205,207,209,211,213,215,217,219,221,223,225,227])
        
        self.imageUpdate = False
        self.l_scan = 0
        self.l_pair = 0
        self.l_antenna = 0
        self.r_scan = 0
        self.r_pair = 0
        self.r_antenna = 0
        self.antNumbers = NP.concatenate((self.receiver1, self.receiver3, self.receiver5, self.receiver7))
        self.antPairs = NP.concatenate((NP.array((self.receiver1[:-1],self.receiver1[1:])).T, NP.array((self.receiver3[:-1],self.receiver3[1:])).T, 
                                       NP.array((self.receiver5[:-1],self.receiver5[1:])).T, NP.array((self.receiver7[:-1],self.receiver7[1:])).T), axis=0)
        self.redundantVisIndexes = NP.zeros(46)
        for i in range(46):
            self.redundantVisIndexes[i] = self.visIndex(self.antennaIndex(self.antPairs[i][0]), self.antennaIndex(self.antPairs[i][1]))
        
#        self.setGeometry(100,100,1024,700)
        self.openButton = QtWidgets.QPushButton('Open...', self)
        self.openButton.clicked.connect(self.onOpen)

        self.leftCanvas = ResponseCanvas(self)
        self.rightCanvas = ResponseCanvas(self)
#
#        self.timeList = QtWidgets.QComboBox(self)
#        self.timeList.currentIndexChanged.connect(self.onTimeListSelected)

        self.leftPlotXLabel = QtWidgets.QLabel(self)
        self.leftPlotXLabel.setText('X: ')
        self.leftPlotX = QtWidgets.QComboBox(self)
        self.leftPlotX.addItem('Time')
        self.leftPlotX.addItem('Antenna')
        self.leftPlotX.addItem('Baseline')
        self.leftPlotX.currentIndexChanged.connect(self.onLeftPlotXSelected)
        self.leftPlotYLabel = QtWidgets.QLabel(self)
        self.leftPlotYLabel.setText('Y: ')
        self.leftPlotY = QtWidgets.QComboBox(self)
        self.leftPlotY.currentIndexChanged.connect(self.onLeftPlotYSelected)
        self.leftPlotY.addItem('Amp (Vis)')
        self.leftPlotY.addItem('Phase')
        self.leftPlotY.addItem('Real')
        self.leftPlotY.addItem('Imag')
        self.leftPlotY.addItem('Amp (Ant)')
        self.leftPlotY.addItem('Mean Amp (Ant)')
        
        self.leftPolarization = QtWidgets.QComboBox(self)
        self.leftPolarization.addItem('LCP')
        self.leftPolarization.addItem('RCP')
        self.leftPolarization.currentIndexChanged.connect(self.onLeftPolarizationSelected)
        
        self.leftOverplotButton = QtWidgets.QPushButton('Overplot', self)
        self.leftOverplotButton.setCheckable(True)
        self.leftOverplotButton.setChecked(False)
        self.leftOverplotButton.clicked.connect(self.onLeftOverplot)
        
        self.rightPolarization = QtWidgets.QComboBox(self)
        self.rightPolarization.addItem('LCP')
        self.rightPolarization.addItem('RCP')
        self.rightPolarization.currentIndexChanged.connect(self.onRightPolarizationSelected)
        
        self.rightOverplotButton = QtWidgets.QPushButton('Overplot', self)
        self.rightOverplotButton.setCheckable(True)
        self.rightOverplotButton.setChecked(False)
        self.rightOverplotButton.clicked.connect(self.onRightOverplot)
        
        self.leftAnt = QtWidgets.QSpinBox(self)
        self.leftAnt.setPrefix('Antenna \'81\' ')
        self.leftAnt.setRange(0,63)
        self.leftAnt.setValue(0)
        self.leftAnt.valueChanged.connect(self.onLeftAntChanged)
        
        self.leftBaseline = QtWidgets.QSpinBox(self)
        self.leftBaseline.setPrefix('Baseline \'81-83\'')
        self.leftBaseline.setRange(0,60)
        self.leftBaseline.setValue(0)
        self.leftBaseline.valueChanged.connect(self.onLeftBaselineChanged)
        
        self.leftScan = QtWidgets.QSpinBox(self)
        self.leftScan.setPrefix('Scan: ')
        self.leftScan.setValue(0)
        self.leftScan.valueChanged.connect(self.onLeftScanChanged)
        
        self.rightAnt = QtWidgets.QSpinBox(self)
        self.rightAnt.setPrefix('Antenna \'81\' ')
        self.rightAnt.setRange(0,63)
        self.rightAnt.setValue(0)
        self.rightAnt.valueChanged.connect(self.onRightAntChanged)
        
        self.rightBaseline = QtWidgets.QSpinBox(self)
        self.rightBaseline.setPrefix('Baseline \'81-83\'')
        self.rightBaseline.setRange(0,60)
        self.rightBaseline.setValue(0)
        self.rightBaseline.valueChanged.connect(self.onRightBaselineChanged)
        
        self.rightScan = QtWidgets.QSpinBox(self)
        self.rightScan.setPrefix('Scan: ')
        self.rightScan.setValue(0)
        self.rightScan.valueChanged.connect(self.onRightScanChanged)
        
        
        self.rightPlotXLabel = QtWidgets.QLabel(self)
        self.rightPlotXLabel.setText('X: ')
        self.rightPlotX = QtWidgets.QComboBox(self)
        self.rightPlotX.addItem('Time')
        self.rightPlotX.addItem('Antenna')
        self.rightPlotX.addItem('Baseline')
        self.rightPlotX.currentIndexChanged.connect(self.onRightPlotXSelected)
        self.rightPlotYLabel = QtWidgets.QLabel(self)
        self.rightPlotYLabel.setText('Y: ')
        self.rightPlotY = QtWidgets.QComboBox(self)
        self.rightPlotY.currentIndexChanged.connect(self.onRightPlotYSelected)
        self.rightPlotY.addItem('Amp (Vis)')
        self.rightPlotY.addItem('Phase')
        self.rightPlotY.addItem('Real')
        self.rightPlotY.addItem('Imag')
        self.rightPlotY.addItem('Amp (Ant)')
        self.rightPlotY.addItem('Mean Amp (Ant)')
        
        self.openButton.setGeometry(0,0,80,25)
        self.leftPlotXLabel.setGeometry(10,25,10,25)
        self.leftPlotX.setGeometry(25,25,100,25)
        self.leftPlotYLabel.setGeometry(10,50,10,25)
        self.leftPlotY.setGeometry(25,50,100,25)
        self.leftPolarization.setGeometry(130,25,80,25)
        self.leftOverplotButton.setGeometry(130,50,80,25)
        self.leftAnt.setGeometry(220,25,150,25)
        self.leftBaseline.setGeometry(220,50,150,25)
        self.leftScan.setGeometry(380,25,150,25)
        self.leftCanvas.setGeometry(0,75,1000,270)
        
        self.rightPlotXLabel.setGeometry(10,350,10,25)
        self.rightPlotX.setGeometry(25,350,100,25)
        self.rightPlotYLabel.setGeometry(10,375,10,25)
        self.rightPlotY.setGeometry(25,375,100,25)
        self.rightPolarization.setGeometry(130,350,80,25)
        self.rightOverplotButton.setGeometry(130,375,80,25)
        self.rightAnt.setGeometry(220,350,150,25)
        self.rightBaseline.setGeometry(220,375,150,25)
        self.rightScan.setGeometry(380,350,150,25)
        self.rightCanvas.setGeometry(0,400,1000,270)
        



    def onOpen(self):
        fitsNames, _ = QtWidgets.QFileDialog.getOpenFileNames(self) 
        print(fitsNames)
        if fitsNames[0]:
            nfits = fits.open(fitsNames[0])
            self.time = nfits[1].data['time']
            samplesNumber = self.time.shape[1]
            amplcp = nfits[1].data['amp_lcp']
            amplitudeNumber = amplcp.shape[1]//samplesNumber
            self.ampLcp = amplcp.reshape(1,samplesNumber,amplitudeNumber)
            rawlcp = nfits[1].data['vis_lcp']
            visibilityNumber = rawlcp.shape[1]//samplesNumber
            unpackFormat = ">%dl" % (rawlcp.shape[1]*2)
            ff = NP.array(struct.unpack(unpackFormat,rawlcp))
            fff = ff.reshape(rawlcp.shape[1],2)
            fff_real_l = fff[:,0].reshape(1,samplesNumber,visibilityNumber)
            fff_imag_l = fff[:,1].reshape(1,samplesNumber,visibilityNumber)
            
            amprcp = nfits[1].data['amp_rcp']
            self.ampRcp = amprcp.reshape(1,samplesNumber,amplitudeNumber)
            rawrcp = nfits[1].data['vis_rcp']
            visibilityNumber = rawrcp.shape[1]//samplesNumber
            unpackFormat = ">%dl" % (rawrcp.shape[1]*2)
            ff = NP.array(struct.unpack(unpackFormat,rawrcp))
            fff = ff.reshape(rawrcp.shape[1],2)
            fff_real_r = fff[:,0].reshape(1,samplesNumber,visibilityNumber)
            fff_imag_r = fff[:,1].reshape(1,samplesNumber,visibilityNumber)
            
            self.currentScan = 0
            self.currentFrequencyChannel = nfits[1].data['frequency'][0]
#            self.frequencyChannel.setValue(self.currentFrequencyChannel)
#            self.scan.setRange(0, samplesNumber-1)
#            self.scan.setValue(self.currentScan)
#            self.calibScan.setRange(0, samplesNumber-1)

            self.setWindowTitle('SRH Fits View 3-6: ' + fitsNames[0])

            for fitsName in fitsNames[1:]:
                nfits = fits.open(fitsName)
                self.time = NP.concatenate((self.time, nfits[1].data['time']), axis = 1)
                samplesNumber = nfits[1].data['time'].shape[1]
                rawlcp = nfits[1].data['vis_lcp']
                visibilityNumber = rawlcp.shape[1]//samplesNumber
                unpackFormat = ">%dl" % (rawlcp.shape[1]*2)
                ff = NP.array(struct.unpack(unpackFormat,rawlcp))
                fff = ff.reshape(rawlcp.shape[1],2)
                fff_real_l = NP.concatenate((fff_real_l, fff[:,0].reshape(1,samplesNumber,visibilityNumber)), axis = 1)
                fff_imag_l = NP.concatenate((fff_imag_l, fff[:,1].reshape(1,samplesNumber,visibilityNumber)), axis = 1)
                
                amplcp = nfits[1].data['amp_lcp']
                amplitudeNumber = amplcp.shape[1]//samplesNumber
                self.ampLcp = NP.concatenate((self.ampLcp, amplcp.reshape(1,samplesNumber,amplitudeNumber)), axis = 1)
                
                rawrcp = nfits[1].data['vis_rcp']
                unpackFormat = ">%dl" % (rawrcp.shape[1]*2)
                ff = NP.array(struct.unpack(unpackFormat,rawrcp))
                fff = ff.reshape(rawrcp.shape[1],2)
                fff_real_r = NP.concatenate((fff_real_r, fff[:,0].reshape(1,samplesNumber,visibilityNumber)), axis = 1)
                fff_imag_r = NP.concatenate((fff_imag_r, fff[:,1].reshape(1,samplesNumber,visibilityNumber)), axis = 1)
                
                amprcp = nfits[1].data['amp_rcp']
                self.ampRcp = NP.concatenate((self.ampRcp, amprcp.reshape(1,samplesNumber,amplitudeNumber)), axis = 1)
                
            self.leftScan.setRange(0, self.time.shape[1])
            self.rightScan.setRange(0, self.time.shape[1])
#            self.calibScan.setRange(0, self.time.shape[1])
        
        self.visibilitiesLcp = fff_real_l[0] + fff_imag_l[0] * 1j
        self.ampLcp = self.ampLcp[0]
        self.visibilitiesRcp = fff_real_r[0] + fff_imag_r[0] * 1j
        self.ampRcp = self.ampRcp[0]
#        for tim in self.time[0]:
#            fTime = QtCore.QTime(0,0)
#            fTime = fTime.addMSecs(tim * 1000)
#            self.timeList.addItem(fTime.toString('hh:mm:ss'))
        self.imageUpdate = True
        self.buildImageLeft()
        self.buildImageRight()

    def updateCanvas(self, scan):
        self.scan.setValue(scan)
        self.buildImage()
        return self.rcpCanvas.imageObject,

                               
#application = QtWidgets.QApplication.instance();
#if not application:
#    application = QtWidgets.QApplication(sys.argv);
#    
#if sys.platform == 'linux':
#    font = QtGui.QFont()
#    application.setFont(QtGui.QFont(font.defaultFamily(),8));

fv = FitsView();
fv.setWindowTitle('SRH Fits View 3-6')
fv.show();
#sys.exit(application.exec_());