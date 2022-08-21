#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 03:49:09 2019

@author: svlesovoi
"""

import sys;
from PyQt5 import QtGui, QtCore, QtWidgets, QtNetwork
from PyQt5.QtWidgets import QFileDialog
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as NP
from srhFitsFile import SrhFitsFile
from skimage.transform import warp, AffineTransform
from astropy.io import fits
from astropy.time import Time, TimeDelta
import base2uvw as bl2uvw
import pylab as PL
from scipy.signal import argrelextrema
import ctypes
import struct

class SrhRTVCanvas(FigureCanvas):
    mouseSignal = QtCore.pyqtSignal(float, float, name = 'xyChanged')
    
    def EW_format(self, a, pos):
        return '%2d' % (49 + a)

    def S_format(self, a, pos):
        return '%03d' % (192 - a*100)
    
    def mouse_moved(self, event):
        1
#        self.mouseSignal.emit(event.xdata, event.ydata)
        
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.subplot = self.fig.add_subplot(111)
        self.subplot.xaxis.set_visible(True)
        self.subplot.yaxis.set_visible(True)
        self.subplot.xaxis.set_major_formatter(PL.FuncFormatter(self.EW_format))
        self.subplot.xaxis.set_major_locator(PL.MultipleLocator(1))
        self.subplot.yaxis.set_major_formatter(PL.FuncFormatter(self.S_format))
        self.subplot.yaxis.set_major_locator(PL.MultipleLocator(0.01))
        self.subplot.set_ylim(0,0.2)
        self.subplot.grid()
        
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

    def setCurve(self, curve):
        self.plotObject[0].sety_data(curve)
        self.draw()

    def imshow(self, array, arrayMin, arrayMax):
        self.imageObject = self.subplot.imshow(array, vmin = arrayMin, vmax = arrayMax, cmap=self.cmap)
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

class SrhRealTimeView(QtWidgets.QMainWindow):
    def onDataRequest(self, request):
        if (request):
            self.dataRequestTimer.start(280)
        else:
            self.dataRequestTimer.stop()
            
    def onDataRequestTimer(self):
        self.dataSocket.writeDatagram(b'SRH_DATA_REQUEST', QtNetwork.QHostAddress('192.168.0.251'), 16565)

    def onDataSocketReadyToRead(self):
        uvBytes = self.dataSocket.readDatagram(512 * ctypes.sizeof(ctypes.c_float))
        uvData = struct.unpack('<512f', uvBytes[0])
        visArr = NP.array(uvData)
        visArr = NP.reshape(visArr, (16,32))
       
        if (self.firstTimePlot == False):
            for i in range(len(self.uvCanvas.subplot.axes.lines)):
                self.uvCanvas.subplot.axes.lines[0].remove()
        for i in range(16):
            self.uvCanvas.subplot.plot(visArr[i] + 0.01*i)
            self.uvCanvas.subplot.plot(visArr[i] + 0.01*i, '.')
        self.firstTimePlot = False
        self.uvCanvas.redraw()

    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self,parent)
        self.firstTimePlot = True
        self.setGeometry(100,100,1024,700)
        self.uvCanvas = SrhRTVCanvas(self)
        self.uvCanvas.setGeometry(10,50,1000,512)
        self.dataRequestButton = QtWidgets.QPushButton('Data request', self)
        self.dataRequestButton.setCheckable(True)
        self.dataRequestButton.clicked.connect(self.onDataRequest)
        self.dataRequestTimer = QtCore.QTimer(self)
        self.dataRequestTimer.timeout.connect(self.onDataRequestTimer)
        self.dataSocket=QtNetwork.QUdpSocket()
        self.dataSocket.bind(QtNetwork.QHostAddress.Any,0)
        self.dataSocket.readyRead.connect(self.onDataSocketReadyToRead)

application = QtWidgets.QApplication.instance();
if not application:
    application = QtWidgets.QApplication(sys.argv);
    
if sys.platform == 'linux':
    font = QtGui.QFont()
    application.setFont(QtGui.QFont(font.defaultFamily(),8));

phaseEdit = SrhRealTimeView();
phaseEdit.setWindowTitle('SRH RTV')
phaseEdit.show();
sys.exit(application.exec_());
