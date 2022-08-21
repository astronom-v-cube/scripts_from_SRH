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
from scipy.signal import detrend
import ctypes
import struct
from scipy import signal
from BadaryRAO import BadaryRAO
import datetime
from astropy import coordinates
from scipy.stats import linregress
from AntennaDelayMatrix import AntennaDelayMatrix

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

    def imshowAutoscale(self, array):
        self.imageObject = self.subplot.imshow(array, cmap=self.cmap)
        self.draw()

    def plot(self, data):
        self.plotObject = self.subplot.plot(data)
        self.draw()
    
    def plotXY(self, curveX, curveY):
        self.subplot.plot(curveX, curveY)
        self.draw()
    
    def scatter(self, data):
        self.plotObject = self.subplot.plot(data, '.', markersize=3)
        self.draw()
    
    def scatterXY(self, X, Y):
        self.plotObject = self.subplot.plot(X, Y, '.', markersize=3)
        self.draw()
    
    def clear(self):
        self.subplot.cla()
        self.draw()
    
    def redraw(self):
        self.draw()
    
    def setColormap(self, cmap):
        self.cmap = cmap

def hhmmssuu_format(t):
  hh = int(t / 3600.)
  t -= hh*3600.
  mm = int(t / 60.)
  t -= mm*60.
  return '%02d:%02d:%02.3f' % (hh,mm,t);

class SrhRealTimeView(QtWidgets.QMainWindow):
    def onDataRequest(self, request):
        if (request):
            self.dataSocket.writeDatagram(b'SRH_DATA_REQUEST', QtNetwork.QHostAddress('192.168.0.251'), 16565)
            self.freqDrSocket.writeDatagram(b'START_FREQ      ', QtNetwork.QHostAddress('192.168.0.191'), 32)
            self.freqCorrSocket.writeDatagram(b'START_FREQ      ', QtNetwork.QHostAddress('192.168.0.190'), 32)
        else:
            self.dataSocket.writeDatagram(b'SRH_DATA_STOP   ', QtNetwork.QHostAddress('192.168.0.251'), 16565)
            self.freqDrSocket.writeDatagram(b'STOP_FREQ       ', QtNetwork.QHostAddress('192.168.0.191'), 32)
            self.freqCorrSocket.writeDatagram(b'STOP_FREQ       ', QtNetwork.QHostAddress('192.168.0.190'), 32)
            
    def onSendAntennasPhases(self):
        antennaPhases = NP.zeros(48)
        nAntennaPhases = NP.zeros(48,dtype='<h')
        meanPairPhases = NP.zeros(46)
        for i in range(46):
            meanPairPhases[i] = NP.mean(NP.unwrap(self.phaVsFreqArray[i,:]))

#        antennaPhases[31] = -self.pha192_64_VsFreqArray.mean()
#        antennaPhases[32] = -self.pha192_65_VsFreqArray.mean()
        for ant in range(15):
            antennaPhases[ant + 1] = -meanPairPhases[ant] + antennaPhases[ant]
            antennaPhases[30 - ant] =  antennaPhases[31 - ant] + meanPairPhases[29 - ant]
            antennaPhases[33 + ant] = -meanPairPhases[31 + ant] + antennaPhases[32 + ant]
            
        antennaPhases[16:32] -= self.pha192_64_VsFreqArray.mean()
        antennaPhases[32:]   -= self.pha192_65_VsFreqArray.mean()

        nAntennaPhases[:] = NP.rad2deg(NP.fmod(antennaPhases[:] + NP.pi*10,2*NP.pi))

        self.dataSocket.writeDatagram(b'SRH_ANT_PHASES  ', QtNetwork.QHostAddress('192.168.0.251'), 16565)
        self.dataSocket.writeDatagram(QtCore.QByteArray(nAntennaPhases.tobytes()),QtNetwork.QHostAddress('192.168.0.251'), 16565)
        print(nAntennaPhases)
        
    def onSendAntennasAmplitudes(self):
        antennaAmps = NP.zeros(48)
        nAntennaAmps = NP.zeros(48,dtype='<h')

        meanPairAmps = NP.zeros(46)
        for i in range(46):
            meanPairAmps[i] = NP.mean(self.ampVsFreqArray[i,:])

        antennaAmps[0] = 1.
        for ant in range(15):
            antennaAmps[ant + 1] = meanPairAmps[ant]**2 / antennaAmps[ant]
            
        antennaAmps[16] = 1.
        for ant in range(31):
            antennaAmps[ant + 17] = meanPairAmps[ant + 16]**2 / antennaAmps[ant + 16]

        nAntennaAmps[:] = int(antennaAmps[:] * 100)

        self.dataSocket.writeDatagram(b'SRH_ANT_AMPS    ', QtNetwork.QHostAddress('192.168.0.251'), 16565)
        self.dataSocket.writeDatagram(QtCore.QByteArray(nAntennaAmps.tobytes()),QtNetwork.QHostAddress('192.168.0.251'), 16565)
        print(nAntennaAmps)
            
    def onSendAntennasDelays(self):
        antennaDelays = NP.zeros(48)
        nAntennaDelays = NP.zeros(48,dtype='<h')

        for ant in range(15):
            antennaDelays[ant + 1] = antennaDelays[ant] + self.phaSlopeVsPair[ant]# + (self.phaSlopePair192_64 + self.phaSlopePair192_65)/2.
            
        for ant in range(31):
            antennaDelays[16 + ant + 1] = antennaDelays[16 + ant] + self.phaSlopeVsPair[15 + ant]
            
        nAntennaDelays[:] = NP.rint(antennaDelays[:]).astype(int)

        self.dataSocket.writeDatagram(b'SRH_ANT_DELAYS  ', QtNetwork.QHostAddress('192.168.0.251'), 16565)
        self.dataSocket.writeDatagram(QtCore.QByteArray(nAntennaDelays.tobytes()),QtNetwork.QHostAddress('192.168.0.251'), 16565)
        print(NP.array_str(nAntennaDelays[0:16], precision=2))
        print(NP.array_str(nAntennaDelays[16:32], precision=2))
        print(NP.array_str(nAntennaDelays[32:], precision=2))

    def onSendAntenna192Delay(self):
        antennaDelays = NP.zeros(48)
        nAntennaDelays = NP.zeros(48,dtype='<h')

        antennaDelays[0:16] -= (self.phaSlopePair192_64 + self.phaSlopePair192_65)/2.
        
        nAntennaDelays[:] = NP.rint(antennaDelays[:]).astype(int)

        self.dataSocket.writeDatagram(b'SRH_ANT_DELAYS  ', QtNetwork.QHostAddress('192.168.0.251'), 16565)
        self.dataSocket.writeDatagram(QtCore.QByteArray(nAntennaDelays.tobytes()),QtNetwork.QHostAddress('192.168.0.251'), 16565)
        print(NP.array_str(nAntennaDelays[0:16], precision=2))
        print(NP.array_str(nAntennaDelays[16:32], precision=2))
        print(NP.array_str(nAntennaDelays[32:], precision=2))

    def onDataRequestTimer(self):
        self.dataSocket.writeDatagram(b'SRH_DATA_REQUEST', QtNetwork.QHostAddress('192.168.0.251'), 16565)

    def onFrequencyIndex(self, index):
        self.frequencyIndex = index
        self.frequencyIndexSpin.setSuffix(', %1.1f GHz'%(1e-9*self.freqList[self.frequencyIndex]))

    def onPairIndex(self, index):
        self.pairIndex = index
        self.pairIndexSpin.setSuffix(self.index2pair[self.pairIndex])
#        self.phaSlopeCanvas.clear()

    def onAntennaIndex(self, index):
        self.antennaIndex = index
        self.antennaIndexSpin.setSuffix(self.index2antenna[self.antennaIndex])

    def onFreqDrSocketReadyToRead(self):
        drMsg = self.freqDrSocket.readDatagram(2*ctypes.sizeof(ctypes.c_double))
        drPacket = struct.unpack('<2d', drMsg[0])
        self.freqDrLabel.setText('%d, MHz'%(drPacket[1]/1e6))
        self.drTimes[1:] = self.drTimes[0:-1]
        self.drFreqs[1:] = self.drFreqs[0:-1]
        self.drTimes[0] = drPacket[0]
        self.drFreqs[0] = drPacket[1]/1e6

    def onFreqCorrSocketReadyToRead(self):
        corrMsg = self.freqCorrSocket.readDatagram(2*ctypes.sizeof(ctypes.c_double))
        corrPacket = struct.unpack('<2d', corrMsg[0])
        self.freqCorrLabel.setText('%d, MHz'%(corrPacket[1]))
        self.corrTimes[1:] = self.corrTimes[0:-1]
        self.corrFreqs[1:] = self.corrFreqs[0:-1]
        self.corrTimes[0] = corrPacket[0]
        self.corrFreqs[0] = corrPacket[1]

    def onDataSocketReadyToRead(self):
        uvSize = 512 + 15 + 31 + 14 + 30# uv, u, v
        uvBytes = self.dataSocket.readDatagram((3 + 2*uvSize) * ctypes.sizeof(ctypes.c_float))
        uvData = struct.unpack('<1207f', uvBytes[0])
        SPairArr = NP.zeros(15, dtype='complex64')
        EWPairArr = NP.zeros(31, dtype='complex64')
        S2PairArr = NP.zeros(14, dtype='complex64')
        EW2PairArr = NP.zeros(30, dtype='complex64')
        visArr = NP.zeros(512, dtype='complex64')
        self.srhControlFreqIndex = int(uvData[0])
        uvFreq = uvData[1]
#        uvTime = uvData[2]
        if (self.srhControlFreqIndex >= 0 and self.srhControlFreqIndex < 32):
            self.freqList[self.srhControlFreqIndex] = uvFreq + 22e6 # + IF1

        if (self.showLcp):
            for i in range(512):
                visArr[i]    = NP.sin(NP.pi/2*uvData[3 + 2*i]) +              1j*NP.sin(NP.pi/2*uvData[3 + 2*i + 1])
            for i in range(15):
                SPairArr[i]  = NP.sin(NP.pi/2*uvData[3 + 2*(i + 512)]) +        1j*NP.sin(NP.pi/2*uvData[3 + 2*(i + 512) + 1])
            for i in range(31):
                EWPairArr[i] = NP.sin(NP.pi/2*uvData[3 + 2*(i + 512 + 15)]) +   1j*NP.sin(NP.pi/2*uvData[3 + 2*(i + 512 + 15) + 1])
            for i in range(14):
                S2PairArr[i]  = NP.sin(NP.pi/2*uvData[3 + 2*(i + 512 + 15 + 31)]) +        1j*NP.sin(NP.pi/2*uvData[3 + 2*(i + 512 + 15 + 31) + 1])
            for i in range(30):
                EW2PairArr[i] = NP.sin(NP.pi/2*uvData[3 + 2*(i + 512 + 15 + 31 + 14)]) +   1j*NP.sin(NP.pi/2*uvData[3 + 2*(i + 512 + 15 + 31 + 14) + 1])
        else:
            for i in range(512):
                visArr[i]    = NP.sin(NP.pi/2*uvData[3 + uvSize + 2*i]) +              1j*NP.sin(NP.pi/2*uvData[3 + uvSize + 2*i + 1])
            for i in range(15):
                SPairArr[i]  = NP.sin(NP.pi/2*uvData[3 + uvSize + 2*(i + 512)]) +        1j*NP.sin(NP.pi/2*uvData[3 + uvSize + 2*(i + 512) + 1])
            for i in range(31):
                EWPairArr[i] = NP.sin(NP.pi/2*uvData[3 + uvSize + 2*(i + 512 + 15)]) +   1j*NP.sin(NP.pi/2*uvData[3 + uvSize + 2*(i + 512 + 15) + 1])
            for i in range(14):
                S2PairArr[i]  = NP.sin(NP.pi/2*uvData[3 + uvSize + 2*(i + 512 + 15 + 31)]) +        1j*NP.sin(NP.pi/2*uvData[3 + uvSize + 2*(i + 512 + 15 + 31) + 1])
            for i in range(30):
                EW2PairArr[i] = NP.sin(NP.pi/2*uvData[3 + uvSize + 2*(i + 512 + 15 + 31 + 14)]) +   1j*NP.sin(NP.pi/2*uvData[3 + uvSize + 2*(i + 512 + 15 + 31 + 14) + 1])
           
        if (self.firstTimePlot == False):
            for i in range(len(self.uvCanvas.subplot.axes.lines)):
                self.uvCanvas.subplot.axes.lines[0].remove()
        
        for i in range(16):
            self.uvCanvas.subplot.plot(0.1*NP.abs(visArr[i*32:(i+1)*32]) + 0.01*i)
        self.uvCanvas.subplot.set_title('frequency %f, tc %s, td %s'%(uvFreq, hhmmssuu_format(self.corrTimes[0]), hhmmssuu_format(self.drTimes[0])))
        self.uvCanvas.redraw()

        self.uvLcp[:,:] = complex(0,0)

        ewPhaLcp = NP.angle(EWPairArr);
        sPhaLcp  = NP.angle(SPairArr);
        ewAmpLcp = NP.abs(EWPairArr);
        sAmpLcp  = NP.abs(SPairArr);
        self.pha192_64_VsFreqArray[self.srhControlFreqIndex] = NP.angle(visArr[15])
        self.pha192_65_VsFreqArray[self.srhControlFreqIndex] = NP.angle(visArr[16])

        self.ewAntPhaLcp, c, d, e = NP.linalg.lstsq(self.ewPhaMatrix,ewPhaLcp)
        self.sAntPhaLcp, c, d, e = NP.linalg.lstsq(self.sPhaMatrix,sPhaLcp)
        
        self.antPhaVsFreqArray[:, self.srhControlFreqIndex] = NP.concatenate((self.sAntPhaLcp[1:], self.ewAntPhaLcp[1:]))
        self.phaVsFreqArray[0:15, self.srhControlFreqIndex] = sPhaLcp
        self.phaVsFreqArray[15:46, self.srhControlFreqIndex] = ewPhaLcp
        self.phaVsFreqArray[46:60, self.srhControlFreqIndex] = NP.angle(S2PairArr)
        self.phaVsFreqArray[60:90, self.srhControlFreqIndex] = NP.angle(EW2PairArr)

        self.ampVsFreqArray[0:15, self.srhControlFreqIndex] = sAmpLcp
        self.ampVsFreqArray[15:46, self.srhControlFreqIndex] = ewAmpLcp
        self.ampVsFreqArray[46:60, self.srhControlFreqIndex] = NP.abs(S2PairArr)
        self.ampVsFreqArray[60:90, self.srhControlFreqIndex] = NP.abs(EW2PairArr)
        if self.srhControlFreqIndex == self.frequencyIndex:
            NP.roll(self.visBuffer,2,1)
            self.visBuffer[:,:,0] = self.ampVsFreqArray * NP.exp(1j*self.phaVsFreqArray)
#        self.ampVsFreqArray[15:29, self.srhControlFreqIndex] = NP.abs(S2PairArr)
#        self.ampVsFreqArray[29:60, self.srhControlFreqIndex] = ewAmpLcp
#        self.ampVsFreqArray[60:90, self.srhControlFreqIndex] = NP.abs(EW2PairArr)

        O = self.lmSize//2 + 1
        for jj in range(16):
            for ii in range(32):
              self.uvLcp[O + jj*2, O + (ii - 16)*2] = visArr[jj*32 +  ii]# * NP.exp(1j*(self.ewAntPhaLcp[ii + 1] - self.sAntPhaLcp[jj + 1]))
              self.uvLcp[O - 2 - jj*2, O + (15 - ii)*2] = NP.conj(self.uvLcp[O + jj*2, O + (ii - 16)*2])

#        if self.frequencyIndex == self.srhControlFreqIndex:
        self.lcp = NP.fft.fft2(NP.roll(NP.roll(self.uvLcp,self.lmSize//2,0),self.lmSize//2,1));
        self.lcp = NP.roll(NP.roll(self.lcp,self.lmSize//2,0),self.lmSize//2,1);

        self.lmCanvas.clear()
        self.lmCanvas.imshowAutoscale(self.lcp.real)
        self.lmCanvas.redraw()
        
        self.visCanvas.clear()
        self.visCanvas.subplot.set_ylim(0., 1.)
        self.visCanvas.plot(self.ampVsFreqArray[:, self.frequencyIndex])
        self.visCanvas.plotXY([14,14],[0,1])
        self.visCanvas.plotXY([45,45],[0,1])
        self.visCanvas.plotXY([59,59],[0,1])
        self.visCanvas.subplot.grid()
        self.visCanvas.redraw()

        self.ampCanvas.clear()
        self.ampCanvas.subplot.set_xlim(4e9,8e9)
        self.ampCanvas.subplot.set_ylim(0,8)
        self.ampCanvas.scatterXY(self.freqList, 10*self.ampVsFreqArray[self.pairIndex, :])
        self.ampCanvas.subplot.set_title('freq %.0f MHz, abs %f'%(self.freqList[self.srhControlFreqIndex]*1e-6, 10*self.ampVsFreqArray[self.pairIndex, self.srhControlFreqIndex]))
        self.ampCanvas.subplot.grid()
        self.ampCanvas.plotXY([uvFreq, uvFreq], [0, 10])
        self.ampCanvas.redraw()

        self.phaCanvas.clear()
        self.phaCanvas.subplot.set_xlim(4e9,8e9)
        self.phaCanvas.subplot.set_ylim(-5,5)
        self.phaCanvas.scatterXY(self.freqList, NP.unwrap(self.pha192_64_VsFreqArray))
        self.phaCanvas.scatterXY(self.freqList, NP.unwrap(self.pha192_65_VsFreqArray))
        self.phaCanvas.scatterXY(self.freqList, NP.unwrap(self.phaVsFreqArray[self.pairIndex, :]))
        self.phaCanvas.subplot.set_title('mean %f'%(NP.rad2deg(self.phaVsFreqArray[self.pairIndex, :].mean())))
        self.phaCanvas.subplot.grid()
        self.phaCanvas.plotXY([uvFreq, uvFreq], [-5, 5])
        self.phaCanvas.redraw()

        if (self.srhControlFreqIndex == 0):
            phaMeanVsPair = []
            ampMeanVsPair = []
            for p in range(46):
                phaSlope, intercept, r_value, p_value, std_err = linregress(self.freqList, NP.unwrap(self.phaVsFreqArray[p, :]))
                self.phaSlopeVsPair[p] = phaSlope/(2*NP.pi)*1e12 #picosecond
                phaMeanVsPair.append(NP.rad2deg(NP.mean(NP.unwrap(self.phaVsFreqArray[p, :]))))
                ampMeanVsPair.append(.1*NP.mean(self.ampVsFreqArray[p, :]))
            
            phaSlope, intercept, r_value, p_value, std_err = linregress(self.freqList, NP.unwrap(self.pha192_64_VsFreqArray))
            self.phaSlopePair192_64 = phaSlope/(2*NP.pi)*1e12 #picosecond
            phaSlope, intercept, r_value, p_value, std_err = linregress(self.freqList, NP.unwrap(self.pha192_65_VsFreqArray))
            self.phaSlopePair192_65 = phaSlope/(2*NP.pi)*1e12 #picosecond
            
            self.phaSlopeCanvas.clear()
            self.phaSlopeCanvas.subplot.set_ylim(-500, 500)
            self.phaSlopeCanvas.plot(self.phaSlopeVsPair)
            self.phaSlopeCanvas.scatter(self.phaSlopeVsPair[0:16])
            self.phaSlopeCanvas.scatterXY(NP.linspace(15,45,31),self.phaSlopeVsPair[15:46])
            self.phaSlopeCanvas.plot(phaMeanVsPair)
            self.phaSlopeCanvas.subplot.grid()
            self.phaSlopeCanvas.subplot.set_title('std slope %f'%(std_err*1e10))
            self.phaSlopeCanvas.redraw()
        
        self.firstTimePlot = False
        uvBytes = 0
        uvData = 0
        SPairArr = 0
        EWPairArr = 0
        
    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self,parent)
        self.index2antenna = {}
        self.index2antenna[0] = ', 192'
        self.index2antenna[1] = ', 191'
        self.index2antenna[2] = ', 190'
        self.index2antenna[3] = ', 189'
        self.index2antenna[4] = ', 188'
        self.index2antenna[5] = ', 187'
        self.index2antenna[6] = ', 186'
        self.index2antenna[7] = ', 185'
        self.index2antenna[8] = ', 184'
        self.index2antenna[9] = ', 183'
        self.index2antenna[10] = ', 182'
        self.index2antenna[11] = ', 181'
        self.index2antenna[12] = ', 180'
        self.index2antenna[13] = ', 179'
        self.index2antenna[14] = ', 178'
        self.index2antenna[15] = ', 177'

        self.index2antenna[16] = ', 49'
        self.index2antenna[17] = ', 50'
        self.index2antenna[18] = ', 51'
        self.index2antenna[19] = ', 52'
        self.index2antenna[20] = ', 53'
        self.index2antenna[21] = ', 54'
        self.index2antenna[22] = ', 55'
        self.index2antenna[23] = ', 56'
        self.index2antenna[24] = ', 57'
        self.index2antenna[25] = ', 58'
        self.index2antenna[26] = ', 59'
        self.index2antenna[27] = ', 60'
        self.index2antenna[28] = ', 61'
        self.index2antenna[29] = ', 62'
        self.index2antenna[30] = ', 63'
        self.index2antenna[31] = ', 64'

        self.index2antenna[32] = ', 65'
        self.index2antenna[33] = ', 66'
        self.index2antenna[34] = ', 67'
        self.index2antenna[35] = ', 68'
        self.index2antenna[36] = ', 69'
        self.index2antenna[37] = ', 70'
        self.index2antenna[38] = ', 71'
        self.index2antenna[39] = ', 72'
        self.index2antenna[40] = ', 73'
        self.index2antenna[41] = ', 74'
        self.index2antenna[42] = ', 75'
        self.index2antenna[43] = ', 76'
        self.index2antenna[44] = ', 77'
        self.index2antenna[45] = ', 78'
        self.index2antenna[46] = ', 79'
        self.index2antenna[47] = ', 80'

        self.index2pair = {}
        self.index2pair[0] = ', 192-191'
        self.index2pair[1] = ', 191-190'
        self.index2pair[2] = ', 190-189'
        self.index2pair[3] = ', 189-188'
        self.index2pair[4] = ', 188-187'
        self.index2pair[5] = ', 187-186'
        self.index2pair[6] = ', 186-185'
        self.index2pair[7] = ', 185-184'
        self.index2pair[8] = ', 184-183'
        self.index2pair[9] = ', 183-182'
        self.index2pair[10]= ', 182-181'
        self.index2pair[11]= ', 181-180'
        self.index2pair[12]= ', 180-179'
        self.index2pair[13]= ', 179-178'
        self.index2pair[14]= ', 178-177'

        self.index2pair[15]= ', 49-50'
        self.index2pair[16]= ', 50-51'
        self.index2pair[17]= ', 51-52'
        self.index2pair[18]= ', 52-53'
        self.index2pair[19]= ', 53-54'
        self.index2pair[20]= ', 54-55'
        self.index2pair[21]= ', 55-56'
        self.index2pair[22]= ', 56-57'
        self.index2pair[23]= ', 57-58'
        self.index2pair[24]= ', 58-59'
        self.index2pair[25]= ', 59-60'
        self.index2pair[26]= ', 60-61'
        self.index2pair[27]= ', 61-62'
        self.index2pair[28]= ', 62-63'
        self.index2pair[29]= ', 63-64'
        
        self.index2pair[30]= ', 64-65'              
        self.index2pair[31]= ', 65-66'
        self.index2pair[32]= ', 66-67'
        self.index2pair[33]= ', 67-68'
        self.index2pair[34]= ', 68-69'
        self.index2pair[35]= ', 69-70'
        self.index2pair[36]= ', 70-71'
        self.index2pair[37]= ', 71-72'
        self.index2pair[38]= ', 72-73'
        self.index2pair[39]= ', 73-74'
        self.index2pair[40]= ', 74-75'
        self.index2pair[41]= ', 75-76'
        self.index2pair[42]= ', 76-77'
        self.index2pair[43]= ', 77-78'
        self.index2pair[44]= ', 78-79'
        self.index2pair[45]= ', 79-80'

        self.index2pair[46] = ', 192-190'
        self.index2pair[47] = ', 191-189'
        self.index2pair[48] = ', 190-188'
        self.index2pair[49] = ', 189-187'
        self.index2pair[50] = ', 188-186'
        self.index2pair[51] = ', 187-185'
        self.index2pair[52] = ', 186-184'
        self.index2pair[53] = ', 185-183'
        self.index2pair[54] = ', 184-182'
        self.index2pair[55] = ', 183-181'
        self.index2pair[56] = ', 182-180'
        self.index2pair[57] = ', 181-179'
        self.index2pair[58] = ', 180-178'
        self.index2pair[59] = ', 179-177'

        self.index2pair[60] = ', 49-51'
        self.index2pair[61] = ', 50-52'
        self.index2pair[62] = ', 51-53'
        self.index2pair[63] = ', 52-54'
        self.index2pair[64] = ', 53-55'
        self.index2pair[65] = ', 54-56'
        self.index2pair[66] = ', 55-57'
        self.index2pair[67] = ', 56-58'
        self.index2pair[68] = ', 57-59'
        self.index2pair[69] = ', 58-60'
        self.index2pair[70] = ', 59-61'
        self.index2pair[71] = ', 60-62'
        self.index2pair[72] = ', 61-63'
        self.index2pair[73] = ', 62-64'
        self.index2pair[74] = ', 63-65'
        self.index2pair[75] = ', 64-66'
        self.index2pair[76] = ', 65-67'
        self.index2pair[77] = ', 66-68'
        self.index2pair[78] = ', 67-69'
        self.index2pair[79] = ', 68-70'
        self.index2pair[80] = ', 69-71'
        self.index2pair[81] = ', 70-72'
        self.index2pair[82] = ', 71-73'
        self.index2pair[83] = ', 72-74'
        self.index2pair[84] = ', 73-75'
        self.index2pair[85] = ', 74-76'
        self.index2pair[86] = ', 75-77'
        self.index2pair[87] = ', 76-78'
        self.index2pair[88] = ', 77-79'
        self.index2pair[89] = ', 78-80'

        self.freqList = NP.zeros(32)
        
        self.firstTimePlot = False
        self.setGeometry(100,100,1600,950)
        self.showLcp = True

        self.dataRequestButton = QtWidgets.QPushButton('Data request', self)
        self.dataRequestButton.setCheckable(True)
        self.dataRequestButton.clicked.connect(self.onDataRequest)
        self.dataRequestButton.setGeometry(0,0,80,30)

        self.sendAntennasPhasesButton = QtWidgets.QPushButton('Send phases', self)
        self.sendAntennasPhasesButton.clicked.connect(self.onSendAntennasPhases)
        self.sendAntennasPhasesButton.setGeometry(600,0,80,30)

        self.sendAntennasDelaysButton = QtWidgets.QPushButton('Send delays', self)
        self.sendAntennasDelaysButton.clicked.connect(self.onSendAntennasDelays)
        self.sendAntennasDelaysButton.setGeometry(690,0,80,30)

        self.sendAntennas192DelayButton = QtWidgets.QPushButton('Send 192 delay', self)
        self.sendAntennas192DelayButton.clicked.connect(self.onSendAntenna192Delay)
        self.sendAntennas192DelayButton.setGeometry(780,0,80,30)

        self.frequencyIndexSpin = QtWidgets.QSpinBox(self)
        self.frequencyIndexSpin.setRange(0,31)
        self.frequencyIndexSpin.setGeometry(90,0,90,30)
        self.frequencyIndexSpin.valueChanged.connect(self.onFrequencyIndex)

        self.pairIndexSpin = QtWidgets.QSpinBox(self)
        self.pairIndexSpin.setSuffix(', 192-191')
        self.pairIndexSpin.setRange(0,89)
        self.pairIndexSpin.setGeometry(180,0,100,30)
        self.pairIndexSpin.valueChanged.connect(self.onPairIndex)

        self.antennaIndex = 0
        self.antennaIndexSpin = QtWidgets.QSpinBox(self)
        self.antennaIndexSpin.setSuffix(', 192')
        self.antennaIndexSpin.setRange(0,47)
        self.antennaIndexSpin.setGeometry(280,0,100,30)
        self.antennaIndexSpin.valueChanged.connect(self.onAntennaIndex)

        self.freqDrLabel = QtWidgets.QLabel(self)
        self.freqDrLabel.setGeometry(290,0,70,30)

        self.freqCorrLabel = QtWidgets.QLabel(self)
        self.freqCorrLabel.setGeometry(290,20,70,30)

        self.uvCanvas = SrhRTVCanvas(self)
        self.uvCanvas.setGeometry(10,30,1000,512)
        self.dataRequestTimer = QtCore.QTimer(self)
        self.dataRequestTimer.timeout.connect(self.onDataRequestTimer)

        self.dataSocket=QtNetwork.QUdpSocket()
        self.dataSocket.bind(QtNetwork.QHostAddress.Any,0)
        self.dataSocket.readyRead.connect(self.onDataSocketReadyToRead)

        self.freqDrSocket=QtNetwork.QUdpSocket()
        self.freqDrSocket.bind(QtNetwork.QHostAddress.Any,0)
        self.freqDrSocket.readyRead.connect(self.onFreqDrSocketReadyToRead)

        self.freqCorrSocket=QtNetwork.QUdpSocket()
        self.freqCorrSocket.bind(QtNetwork.QHostAddress.Any,0)
        self.freqCorrSocket.readyRead.connect(self.onFreqCorrSocketReadyToRead)

        self.lmSize = 256
        self.lmCanvas = SrhRTVCanvas(self)
        self.lmCanvas.setGeometry(1020,30,self.lmSize*2,self.lmSize*2)
        self.uvLcp = NP.zeros((self.lmSize, self.lmSize), dtype='complex')
        self.lcp = NP.zeros((self.lmSize, self.lmSize), dtype='complex')
        self.ewLcpPhaseCorrection = NP.zeros(32)
        self.sLcpPhaseCorrection = NP.zeros(16)
        self.ewPhaMatrix = [ \
            [1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1]];
        self.ewAmpMatrix = [ \
            [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]];
        self.sPhaMatrix = [ \
            [1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1]];
        self.sAmpMatrix = [ \
            [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]];
        
        self.visCanvas = SrhRTVCanvas(self)
        self.visCanvas.setGeometry(0,540,400,400)
        self.ampCanvas = SrhRTVCanvas(self)
        self.ampCanvas.setGeometry(400,540,400,400)
        self.phaCanvas = SrhRTVCanvas(self)
        self.phaCanvas.setGeometry(800,540,400,400)
        self.phaSlopeCanvas = SrhRTVCanvas(self)
        self.phaSlopeCanvas.setGeometry(1200,540,400,400)
#        self.phaSlopeCanvas.clear()
        
        self.visPairArray = NP.zeros((90, 32), dtype='float')
        self.visBuffer = NP.zeros((90, 32, 1000), dtype='complex64')        
    
        self.antPhaVsFreqArray = NP.zeros((48, 32), dtype='float')
        self.phaVsFreqArray = NP.zeros((90, 32), dtype='float')
        self.ampVsFreqArray = NP.zeros((90, 32), dtype='float')
        self.phaSlopeVsPair = NP.zeros(46)
        self.pha192_64_VsFreqArray = NP.zeros((32), dtype='float')
        self.pha192_65_VsFreqArray = NP.zeros((32), dtype='float')
        self.phaSlopePair192_64 = 0.
        self.phaSlopePair192_65 = 0.
        self.srhControlFreqIndex = 0

        self.frequencyIndex = 0
        self.pairIndex = 0
        self.drFreqs = NP.zeros(1000)
        self.corrFreqs = NP.zeros(1000)
        self.drTimes = NP.zeros(1000)
        self.corrTimes = NP.zeros(1000)
        
        self.RAO = BadaryRAO(datetime.date.today())
        self.omegaEarth = coordinates.earth.OMEGA_EARTH.to_value()
        self.delayMatrix = AntennaDelayMatrix(16, self.freqList)

application = QtWidgets.QApplication.instance();
if not application:
    application = QtWidgets.QApplication(sys.argv);

rtView = SrhRealTimeView();
rtView.setWindowTitle('SRH RTV')
rtView.show();
sys.exit(application.exec_());
