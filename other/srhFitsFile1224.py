# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 00:18:47 2016

@author: Sergey
"""

from astropy.io import fits
import numpy as NP
from astropy import coordinates
from astropy import constants
from BadaryRAO import BadaryRAO
from scipy.optimize import least_squares, basinhopping
import sunpy.coordinates
import base2uvw_1224
from skimage.transform import warp, AffineTransform
import scipy.signal
import time
from ZirinTb import ZirinTb
import json

class SrhFitsFile():
    def __init__(self, name, sizeOfUv):
        self.base = 2450
        self.omegaEarth = coordinates.earth.OMEGA_EARTH.to_value()
        self.isOpen = False;
        self.calibIndex = 0;
        self.frequencyChannel = 0;
        self.centerP = 0.;
        self.deltaP = 0.005;
        self.centerQ = 0.;
        self.deltaQ = 0.005;
        self.centerH = 0.;
        self.deltaH = 4.9/3600.*NP.pi;
        self.centerD = 0.;
        self.deltaD = 4.9/3600.*NP.pi;
        self.antNumberEW = 139
        self.antNumberS = 68
        self.averageCalib = False
        self.useNonlinearApproach = True
        self.obsObject = 'Sun'
        self.fringeStopping = False
        
        self.basin = False
        self.centering_ftol = 1e-3
        
        
        self.sizeOfUv = sizeOfUv
        self.baselines = 8
        self.open(name)
        
        self.flagsIndexes = []
        
        self.arcsecPerPixel = 4.91104/2
        
        self.ZirinQSunTb = ZirinTb()
        
        self.convolutionNormCoef = 1.
        
        
        
        
                                    
    def open(self,name):
        try:
            self.hduList = fits.open(name)
            self.isOpen = True
            self.dateObs = self.hduList[0].header['DATE-OBS'] + 'T' + self.hduList[0].header['TIME-OBS']
            self.antennaNumbers = self.hduList[2].data['ant_index']
            self.antennaNumbers = NP.reshape(self.antennaNumbers,self.antennaNumbers.size)
            self.antennaNames = self.hduList[2].data['ant_name']
            self.antennaNames = NP.reshape(self.antennaNames,self.antennaNames.size)
            self.antennaA = self.hduList[4].data['ant_A']
            self.antennaA = NP.reshape(self.antennaA,self.antennaA.size)
            self.antennaB = self.hduList[4].data['ant_B']
            self.antennaB = NP.reshape(self.antennaB,self.antennaB.size)
            self.antX = self.hduList[3].data['ant_X']
            self.antY = self.hduList[3].data['ant_Y']
            self.uvLcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex)
            self.uvRcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex)
            self.freqList = self.hduList[1].data['frequency'];
            self.freqListLength = self.freqList.size;
            self.dataLength = self.hduList[1].data['time'].size // self.freqListLength;
            self.freqTime = self.hduList[1].data['time']
            self.visListLength = self.hduList[1].data['vis_lcp'].size // self.freqListLength // self.dataLength;
            self.visLcp = NP.reshape(self.hduList[1].data['vis_lcp'],(self.freqListLength,self.dataLength,self.visListLength));
            self.visRcp = NP.reshape(self.hduList[1].data['vis_rcp'],(self.freqListLength,self.dataLength,self.visListLength));
            self.visLcp /= float(self.hduList[0].header['VIS_MAX'])
            self.visRcp /= float(self.hduList[0].header['VIS_MAX'])
            self.ampLcp = NP.reshape(self.hduList[1].data['amp_lcp'],(self.freqListLength,self.dataLength,self.antennaNumbers.size));
            self.ampRcp = NP.reshape(self.hduList[1].data['amp_rcp'],(self.freqListLength,self.dataLength,self.antennaNumbers.size));
            ampScale = float(self.hduList[0].header['VIS_MAX']) / 128.
            self.ampLcp = self.ampLcp.astype(float) / ampScale
            self.ampRcp = self.ampRcp.astype(float) / ampScale
            
            self.antZeroRow = []
            for ant in range(self.antNumberEW):
                if ant<69:
                    self.antZeroRow.append(NP.where((self.antennaA==ant) & (self.antennaB==69))[0][0])
                if ant>69:
                    self.antZeroRow.append(NP.where((self.antennaA==69) & (self.antennaB==ant))[0][0])
            # self.antZeroRow = self.hduList[3].data['ant_zero_row'][:97]
            self.RAO = BadaryRAO(self.dateObs.split('T')[0], self.base*1e-3, observedObject = self.obsObject)
            # try:
            #     client = Client('http://ephemeris.rao.istp.ac.ru/?wsdl')
            #     result = client.service.Ephemeride('SSRT','sun',self.dateObs)
            #     self.pAngle = NP.deg2rad(float(result[0]['PAngle']))
            # except:
            self.pAngle = NP.deg2rad(sunpy.coordinates.sun.P(self.dateObs).to_value())
            self.getHourAngle(0)
            
            self.ewAntPhaLcp = NP.zeros((self.freqListLength, self.antNumberEW))
            self.sAntPhaLcp = NP.zeros((self.freqListLength, self.antNumberS))
            self.ewAntPhaRcp = NP.zeros((self.freqListLength, self.antNumberEW))
            self.sAntPhaRcp = NP.zeros((self.freqListLength, self.antNumberS))
            self.ewLcpPhaseCorrection = NP.zeros((self.freqListLength, self.antNumberEW))
            self.ewRcpPhaseCorrection = NP.zeros((self.freqListLength, self.antNumberEW))
            self.sLcpPhaseCorrection = NP.zeros((self.freqListLength, self.antNumberS))
            self.sRcpPhaseCorrection = NP.zeros((self.freqListLength, self.antNumberS))
            self.sSolarPhase = NP.zeros(self.freqListLength)
            self.ewSolarPhase = NP.zeros(self.freqListLength)
            
            self.ewAntAmpLcp = NP.ones((self.freqListLength, self.antNumberEW))
            self.sAntAmpLcp = NP.ones((self.freqListLength, self.antNumberS))
            self.ewAntAmpRcp = NP.ones((self.freqListLength, self.antNumberEW))
            self.sAntAmpRcp = NP.ones((self.freqListLength, self.antNumberS))
            
            self.sLcpStair = NP.zeros(self.freqListLength)
            self.sRcpStair = NP.zeros(self.freqListLength)
            self.ewSlopeLcp = NP.zeros(self.freqListLength)
            self.sSlopeLcp = NP.zeros(self.freqListLength)
            self.ewSlopeRcp = NP.zeros(self.freqListLength)
            self.sSlopeRcp = NP.zeros(self.freqListLength)
            self.diskLevelLcp = NP.ones(self.freqListLength)
            self.diskLevelRcp = NP.ones(self.freqListLength)
            
            x_size = (self.baselines-1)*2 + self.antNumberEW + self.antNumberS
            self.x_ini_lcp = NP.full((self.freqListLength, x_size*2+1), NP.concatenate((NP.ones(x_size+1), NP.zeros(x_size))))
            self.x_ini_rcp = NP.full((self.freqListLength, x_size*2+1), NP.concatenate((NP.ones(x_size+1), NP.zeros(x_size))))
            self.calibrationResultLcp = NP.zeros_like(self.x_ini_lcp)
            self.calibrationResultRcp = NP.zeros_like(self.x_ini_rcp)
            
            self.lcpShift = NP.ones(self.freqListLength) # 0-frequency component in the spectrum
            self.rcpShift = NP.ones(self.freqListLength)
            
        except FileNotFoundError:
            print('File %s  not found'%name);
    
    def append(self,name):
        try:
            hduList = fits.open(name);
            freqTime = hduList[1].data['time']
            dataLength = hduList[1].data['time'].size // self.freqListLength;
            visLcp = NP.reshape(hduList[1].data['vis_lcp'],(self.freqListLength,dataLength,self.visListLength));
            visRcp = NP.reshape(hduList[1].data['vis_rcp'],(self.freqListLength,dataLength,self.visListLength));
            visLcp /= float(self.hduList[0].header['VIS_MAX'])
            visRcp /= float(self.hduList[0].header['VIS_MAX'])
            ampLcp = NP.reshape(hduList[1].data['amp_lcp'],(self.freqListLength,dataLength,self.antennaNumbers.size));
            ampRcp = NP.reshape(hduList[1].data['amp_rcp'],(self.freqListLength,dataLength,self.antennaNumbers.size));
            ampScale = float(self.hduList[0].header['VIS_MAX']) / 128.
            ampLcp = ampLcp.astype(float) / ampScale
            ampRcp = ampRcp.astype(float) / ampScale

            self.freqTime = NP.concatenate((self.freqTime, freqTime), axis = 1)
            self.visLcp = NP.concatenate((self.visLcp, visLcp), axis = 1)
            self.visRcp = NP.concatenate((self.visRcp, visRcp), axis = 1)
            self.ampLcp = NP.concatenate((self.ampLcp, ampLcp), axis = 1)
            self.ampRcp = NP.concatenate((self.ampRcp, ampRcp), axis = 1)
            self.dataLength += dataLength
            hduList.close()

        except FileNotFoundError:
            print('File %s  not found'%name);
            
    def calibrate(self, freq = 'all', phaseCorrect = True, amplitudeCorrect = True, average = 0):
        if freq == 'all':
            self.calculatePhaseCalibration()
            for freq in range(self.freqListLength):
                self.setFrequencyChannel(freq)
                self.vis2uv(scan = self.calibIndex, phaseCorrect=phaseCorrect, amplitudeCorrect=amplitudeCorrect, average=average)
                # self.centerDisk()
        else:
            self.setFrequencyChannel(freq)
            self.solarPhase(freq)
            self.updateAntennaPhase(freq)
            self.vis2uv(scan = self.calibIndex, phaseCorrect=phaseCorrect, amplitudeCorrect=amplitudeCorrect, average=average)
            # self.centerDisk()
            
    def image(self, freq, scan, average = 0, polarization = 'both', phaseCorrect = True, amplitudeCorrect = True, frame = 'heliocentric'):
        self.setFrequencyChannel(freq)    
        self.vis2uv(scan = scan, phaseCorrect=phaseCorrect, amplitudeCorrect=amplitudeCorrect, average=average)
        self.uv2lmImage()
        if frame == 'heliocentric':
            self.lm2Heliocentric()
        if polarization == 'both':
            return NP.flip(self.lcp, 0), NP.flip(self.rcp, 0)
        elif polarization == 'lcp':
            return NP.flip(self.lcp, 0)
        elif polarization == 'rcp':
            return NP.flip(self.rcp, 0)
        
    def saveGains(self, filename):
        currentGainsDict = {}
        currentGainsDict['ewPhaseLcp'] = (self.ewAntPhaLcp + self.ewLcpPhaseCorrection).tolist()
        currentGainsDict['sPhaseLcp'] = (self.sAntPhaLcp + self.sLcpPhaseCorrection).tolist()
        currentGainsDict['ewPhaseRcp'] = (self.ewAntPhaRcp + self.ewRcpPhaseCorrection).tolist()
        currentGainsDict['sPhaseRcp'] = (self.sAntPhaRcp + self.sRcpPhaseCorrection).tolist()
        currentGainsDict['ewAmpLcp'] = self.ewAntAmpLcp.tolist()
        currentGainsDict['sAmpLcp'] = self.sAntAmpLcp.tolist()
        currentGainsDict['ewAmpRcp'] = self.ewAntAmpRcp.tolist()
        currentGainsDict['sAmpRcp'] = self.sAntAmpRcp.tolist()
        with open(filename, 'w') as saveGainFile:
            json.dump(currentGainsDict, saveGainFile)
            
    def loadGains(self, filename):
        with open(filename,'r') as readGainFile:
            currentGains = json.load(readGainFile)
        self.ewAntPhaLcp = NP.array(currentGains['ewPhaseLcp'])
        self.sAntPhaLcp = NP.array(currentGains['sPhaseLcp'])
        self.ewAntPhaRcp = NP.array(currentGains['ewPhaseRcp'])
        self.sAntPhaRcp = NP.array(currentGains['sPhaseRcp'])
        self.ewAntAmpLcp = NP.array(currentGains['ewAmpLcp'])
        self.sAntAmpLcp = NP.array(currentGains['sAmpLcp'])
        self.ewAntAmpRcp = NP.array(currentGains['ewAmpRcp'])
        self.sAntAmpRcp = NP.array(currentGains['sAmpRcp'])

    def changeObject(self, obj):
        self.obsObject = obj
        self.RAO = BadaryRAO(self.dateObs.split('T')[0], observedObject = self.obsObject)
        
    def getHourAngle(self, scan):
        self.hAngle = self.omegaEarth * (self.freqTime[self.frequencyChannel, scan] - self.RAO.culmination)
        if self.hAngle > 1.5*NP.pi:
            self.hAngle -= 2*NP.pi
        return self.hAngle
        
    def setSizeOfUv(self, sizeOfUv):
        self.sizeOfUv = sizeOfUv
        self.uvLcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex);
        self.uvRcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex);
        
    def getDeclination(self):
        return self.RAO.declination

    def getPQScale(self, size, FOV):
        self.cosP = NP.sin(self.hAngle) * NP.cos(self.RAO.declination)
        self.cosQ = NP.cos(self.hAngle) * NP.cos(self.RAO.declination) * NP.sin(self.RAO.observatory.lat) - NP.sin(self.RAO.declination) * NP.cos(self.RAO.observatory.lat)
        FOV_p = 2.*(constants.c / (self.freqList[self.frequencyChannel]*1e3)) / (self.RAO.base*NP.sqrt(1. - self.cosP**2.));
        FOV_q = 2.*(constants.c / (self.freqList[self.frequencyChannel]*1e3)) / (self.RAO.base*NP.sqrt(1. - self.cosQ**2.));
        
        return [int(size*FOV/FOV_p.to_value()), int(size*FOV/FOV_q.to_value())]
        
    def getPQ2HDMatrix(self):
        gP =  NP.arctan(NP.tan(self.hAngle)*NP.sin(self.RAO.declination));
        gQ =  NP.arctan(-(NP.sin(self.RAO.declination) / NP.tan(self.hAngle) + NP.cos(self.RAO.declination) / (NP.sin(self.hAngle)*NP.tan(self.RAO.observatory.lat))));
        
        if self.hAngle > 0:
            gQ = NP.pi + gQ;
        g = gP - gQ;
          
        pqMatrix = NP.zeros((3,3))
        pqMatrix[0, 0] =  NP.cos(gP) - NP.cos(g)*NP.cos(gQ)
        pqMatrix[0, 1] = -NP.cos(g)*NP.cos(gP) + NP.cos(gQ)
        pqMatrix[1, 0] =  NP.sin(gP) - NP.cos(g)*NP.sin(gQ)
        pqMatrix[1, 1] = -NP.cos(g)*NP.sin(gP) + NP.sin(gQ)
        pqMatrix /= NP.sin(g)**2.
        pqMatrix[2, 2] = 1.
        return pqMatrix
        
    def close(self):
        self.hduList.close();
        
    def flag(self, names):
        nameList = names.split(',')
        for i in range(len(nameList)):
            ind = NP.where(self.antennaNames == nameList[i])[0]
            if len(ind):
                self.flagsIndexes.append(int(self.antennaNumbers[ind[0]]))
        self.flagVis = NP.array([], dtype = int)
        for i in range(len(self.flagsIndexes)):
            ind = NP.where(self.antennaA == self.flagsIndexes[i])[0]
            if len(ind):
                self.flagVis = NP.append(self.flagVis, ind)
            ind = NP.where(self.antennaB == self.flagsIndexes[i])[0]
            if len(ind):
                self.flagVis = NP.append(self.flagVis, ind)
        self.visLcp[:,:,self.flagVis] = 0
        self.visRcp[:,:,self.flagVis] = 0
        
    def setBaselinesNumber(self, value):
        self.baselines = value
        x_size = (self.baselines-1)*2 + self.antNumberEW + self.antNumberS
        self.x_ini_lcp = NP.full((self.freqListLength, x_size*2+1), NP.concatenate((NP.ones(x_size+1), NP.zeros(x_size))))
        self.x_ini_rcp = NP.full((self.freqListLength, x_size*2+1), NP.concatenate((NP.ones(x_size+1), NP.zeros(x_size))))
        self.calibrationResultLcp = NP.zeros_like(self.x_ini_lcp)
        self.calibrationResultRcp = NP.zeros_like(self.x_ini_rcp)

    def phaMatrixGenPairsEWN(self, pairs, antNumberEW, antNumberN):
        rowsEW = int(((antNumberEW - 1) + (antNumberEW - pairs))/2 * pairs)
        rowsN = int(((antNumberN) + (antNumberN + 1 - pairs))/2 * pairs)
        colsEW = antNumberEW + pairs
        colsN = antNumberN + pairs
        phaMatrix = NP.zeros((rowsEW + rowsN, colsEW + colsN))
        for pair in range(pairs):
            row0 = int(((antNumberEW - 1) + (antNumberEW - pair))/2 * pair)
            row1 = row0 + (antNumberEW - pair - 1)
            phaMatrix[row0:row1,pair] = 1
            for phaPair in range(antNumberEW - pair - 1):
                phaMatrix[phaPair + row0, phaPair + 2*pairs] = 1
                phaMatrix[phaPair + row0, phaPair + 2*pairs + (pair + 1)] = -1
            row0 = int(((antNumberN) + (antNumberN + 1 - pair))/2 * pair)
            row1 = row0 + (antNumberN - pair)
            phaMatrix[row0 + rowsEW:row1 + rowsEW,pairs + pair] = 1
            for phaPair in range(antNumberN - pair):
                if phaPair == 0:
                    phaMatrix[rowsEW + row0, 2*pairs + 32] = 1
                else:
                    phaMatrix[phaPair + rowsEW + row0, 2*pairs + phaPair + antNumberEW - 1] = 1
                phaMatrix[phaPair + rowsEW + row0, 2*pairs + phaPair + antNumberEW - 1 + (pair + 1)] = -1
        return phaMatrix.copy()
    
    def calculatePhaseCalibration(self, baselinesNumber = 5, lcp = True, rcp = True):
       for freq in range(self.freqListLength):
           self.solarPhase(freq)
           self.updateAntennaPhase(freq, baselinesNumber, lcp = lcp, rcp = rcp)
           
    def calculateAmpCalibration(self, baselinesNumber = 5):
       for freq in range(self.freqListLength):
           self.calculateAmplitude_linear(freq, baselinesNumber)

    def updateAntennaPhase(self, freqChannel, lcp = True, rcp = True):
        if self.useNonlinearApproach:
            if lcp:
                self.calculatePhaseLcp_nonlinear(freqChannel)
            if rcp:
                self.calculatePhaseRcp_nonlinear(freqChannel)
        else:
            self.calculatePhase_linear(freqChannel)
            
    def solarPhase(self, freq):
        u,v,w = base2uvw_1224.base2uvw(self.hAngle, self.RAO.declination, 140, 141)
        baseWave = NP.sqrt(u**2+v**2)*self.freqList[freq]*1e3/constants.c.to_value()
        if baseWave > 120:
            self.nSolarPhase[freq] = NP.pi
        else:
            self.nSolarPhase[freq] = 0
        u,v,w = base2uvw_1224.base2uvw(self.hAngle, self.RAO.declination, 1, 2)
        baseWave = NP.sqrt(u**2+v**2)*self.freqList[freq]*1e3/constants.c.to_value()
        if baseWave > 120:
            self.ewSolarPhase[freq] = NP.pi
        else:
            self.ewSolarPhase[freq] = 0
            
    def calculatePhase_linear(self, freqChannel, baselinesNumber = 1):
        antNumberN = 31
        antNumberEW = 97
        redIndexesN = []
        for baseline in range(1, baselinesNumber+1):
            redIndexesN.append(NP.where((self.antennaA==98-1+baseline) & (self.antennaB==33))[0][0])
            for i in range(antNumberN - baseline):
                redIndexesN.append(NP.where((self.antennaB==98+i) & (self.antennaA==98+i+baseline))[0][0])
    
        redIndexesEW = []
        for baseline in range(1, baselinesNumber+1):
            for i in range(antNumberEW - baseline):
                redIndexesEW.append(NP.where((self.antennaA==1+i) & (self.antennaB==1+i+baseline))[0][0])
                
        phaMatrix = self.phaMatrixGenPairsEWN(baselinesNumber, antNumberEW, antNumberN)
        redundantVisLcp = self.visLcp[freqChannel, self.calibIndex, NP.append(redIndexesEW, redIndexesN)]
        sunVisPhases = NP.zeros(2)
#        if self.freqList[freqChannel] > 4e6:
#            sunVisPhases = NP.array((NP.pi, NP.pi))
        phasesLcp = NP.concatenate((NP.angle(redundantVisLcp), sunVisPhases))
        antPhaLcp, c, d, e = NP.linalg.lstsq(phaMatrix, phasesLcp, rcond=None)
        self.ewAntPhaLcp[freqChannel] = antPhaLcp[baselinesNumber*2:baselinesNumber*2+antNumberEW]
        self.nAntPhaLcp[freqChannel] = antPhaLcp[baselinesNumber*2+antNumberEW:]
        
        redundantVisRcp = self.visRcp[freqChannel, self.calibIndex, NP.append(redIndexesEW, redIndexesN)]
        phasesRcp = NP.concatenate((NP.angle(redundantVisRcp), NP.array((0,0))))
        antPhaRcp, c, d, e = NP.linalg.lstsq(phaMatrix, phasesRcp, rcond=None)
        self.ewAntPhaRcp[freqChannel] = antPhaRcp[baselinesNumber*2:baselinesNumber*2+antNumberEW]
        self.nAntPhaRcp[freqChannel] = antPhaRcp[baselinesNumber*2+antNumberEW:]
        
    def calculatePhaseLcp_nonlinear(self, freqChannel):
        antY_diff = (self.antY[self.antennaB] - self.antY[self.antennaA])/self.base
        antX_diff = (self.antX[self.antennaB] - self.antX[self.antennaA])/self.base
        self.redIndexesS = NP.array(())
        self.redIndexesEW = NP.array(())
        self.redIndexesS_len = []
        self.redIndexesEW_len = []
        for baseline in range(1, self.baselines+1):
            ind = NP.intersect1d(NP.where(NP.abs(antX_diff)==baseline)[0], NP.where(antY_diff == 0)[0])
            self.redIndexesS = NP.append(self.redIndexesS, ind)
            self.redIndexesS_len.append(len(ind))
            ind = NP.intersect1d(NP.where(NP.abs(antY_diff)==baseline)[0], NP.where(antX_diff == 0)[0])
            self.redIndexesEW = NP.append(self.redIndexesEW, ind)
            self.redIndexesEW_len.append(len(ind))

        if self.averageCalib:
            redundantVisS = NP.mean(self.visLcp[freqChannel, :20, self.redIndexesS.astype(int)], axis = 1)
            redundantVisEW = NP.mean(self.visLcp[freqChannel, :20, self.redIndexesEW.astype(int)], axis = 1)
            redundantVisAll = NP.append(redundantVisEW, redundantVisS)
        else:
            redundantVisS = self.visLcp[freqChannel, self.calibIndex, self.redIndexesS.astype(int)]
            redundantVisEW = self.visLcp[freqChannel, self.calibIndex, self.redIndexesEW.astype(int)]
            redundantVisAll = NP.append(redundantVisEW, redundantVisS)

        ls_res = least_squares(self.allGainsFunc_constrained, self.x_ini_lcp[freqChannel], args = (redundantVisAll, self.antNumberEW, self.antNumberS, self.baselines, freqChannel), max_nfev = 400)
        self.calibrationResultLcp[freqChannel] = ls_res['x']
        gains = self.real_to_complex(ls_res['x'][1:])[(self.baselines-1)*2:]
        self.ew_gains_lcp = gains[:self.antNumberEW]
        self.ewAntPhaLcp[freqChannel] = NP.angle(self.ew_gains_lcp)
        self.s_gains_lcp = gains[self.antNumberEW:]
        self.sAntPhaLcp[freqChannel] = NP.angle(self.s_gains_lcp)
        
        norm = NP.mean(NP.abs(gains))#[NP.abs(gains)<NP.median(NP.abs(gains))*0.6]))
        self.ewAntAmpLcp[freqChannel] = NP.abs(self.ew_gains_lcp)/norm
        self.ewAntAmpLcp[freqChannel][self.ewAntAmpLcp[freqChannel]<NP.median(self.ewAntAmpLcp[freqChannel])*0.6] = 1e6
        self.sAntAmpLcp[freqChannel] = NP.abs(self.s_gains_lcp)/norm
        self.sAntAmpLcp[freqChannel][self.sAntAmpLcp[freqChannel]<NP.median(self.sAntAmpLcp[freqChannel])*0.6] = 1e6
        
    def calculatePhaseRcp_nonlinear(self, freqChannel):
        antY_diff = (self.antY[self.antennaB] - self.antY[self.antennaA])/self.base
        antX_diff = (self.antX[self.antennaB] - self.antX[self.antennaA])/self.base
        self.redIndexesS = NP.array(())
        self.redIndexesEW = NP.array(())
        self.redIndexesS_len = []
        self.redIndexesEW_len = []
        for baseline in range(1, self.baselines+1):
            ind = NP.intersect1d(NP.where(NP.abs(antX_diff)==baseline)[0], NP.where(antY_diff == 0)[0])
            self.redIndexesS = NP.append(self.redIndexesS, ind)
            self.redIndexesS_len.append(len(ind))
            ind = NP.intersect1d(NP.where(NP.abs(antY_diff)==baseline)[0], NP.where(antX_diff == 0)[0])
            self.redIndexesEW = NP.append(self.redIndexesEW, ind)
            self.redIndexesEW_len.append(len(ind))

        if self.averageCalib:
            redundantVisS = NP.mean(self.visRcp[freqChannel, :20, self.redIndexesS.astype(int)], axis = 1)
            redundantVisEW = NP.mean(self.visRcp[freqChannel, :20, self.redIndexesEW.astype(int)], axis = 1)
            redundantVisAll = NP.append(redundantVisEW, redundantVisS)
        else:
            redundantVisS = self.visRcp[freqChannel, self.calibIndex, self.redIndexesS.astype(int)]
            redundantVisEW = self.visRcp[freqChannel, self.calibIndex, self.redIndexesEW.astype(int)]
            redundantVisAll = NP.append(redundantVisEW, redundantVisS)

        ls_res = least_squares(self.allGainsFunc_constrained, self.x_ini_rcp[freqChannel], args = (redundantVisAll, self.antNumberEW, self.antNumberS, self.baselines, freqChannel), max_nfev = 400)
        self.calibrationResultRcp[freqChannel] = ls_res['x']
        gains = self.real_to_complex(ls_res['x'][1:])[(self.baselines-1)*2:]
        self.ew_gains_rcp = gains[:self.antNumberEW]
        self.ewAntPhaRcp[freqChannel] = NP.angle(self.ew_gains_rcp)
        self.s_gains_rcp = gains[self.antNumberEW:]
        self.sAntPhaRcp[freqChannel] = NP.angle(self.s_gains_rcp)
        
        norm = NP.mean(NP.abs(gains))#[NP.abs(gains)<NP.median(NP.abs(gains))*0.6]))
        self.ewAntAmpRcp[freqChannel] = NP.abs(self.ew_gains_rcp)/norm
        self.ewAntAmpRcp[freqChannel][self.ewAntAmpRcp[freqChannel]<NP.median(self.ewAntAmpRcp[freqChannel])*0.6] = 1e6
        self.sAntAmpRcp[freqChannel] = NP.abs(self.s_gains_rcp)/norm
        self.sAntAmpRcp[freqChannel][self.sAntAmpRcp[freqChannel]<NP.median(self.sAntAmpRcp[freqChannel])*0.6] = 1e6
     
    def calculateAmplitude_linear(self, freqChannel, baselinesNumber = 3):    
        antNumberN = 31
        antNumberEW = self.antNumberEW
        redIndexesN = []
        for baseline in range(1, baselinesNumber+1):
            redIndexesN.append(NP.where((self.antennaA==98-1+baseline) & (self.antennaB==33))[0][0])
            for i in range(antNumberN - baseline):
                redIndexesN.append(NP.where((self.antennaB==98+i) & (self.antennaA==98+i+baseline))[0][0])
    
        redIndexesEW = []
        for baseline in range(1, baselinesNumber+1):
            for i in range(antNumberEW - baseline):
                redIndexesEW.append(NP.where((self.antennaA==1+i) & (self.antennaB==1+i+baseline))[0][0])
             
        if self.averageCalib:
            redundantVisN = NP.mean(self.visLcp[freqChannel, :20, redIndexesN], axis = 1)
            redundantVisEW = NP.mean(self.visLcp[freqChannel, :20, redIndexesEW], axis = 1)
            redundantVisAllLcp = NP.append(redundantVisEW, redundantVisN)
            redundantVisN = NP.mean(self.visRcp[freqChannel, :20, redIndexesN], axis = 1)
            redundantVisEW = NP.mean(self.visRcp[freqChannel, :20, redIndexesEW], axis = 1)
            redundantVisAllRcp = NP.append(redundantVisEW, redundantVisN)
        else:
            redundantVisN = self.visLcp[freqChannel, self.calibIndex, redIndexesN]
            redundantVisEW = self.visLcp[freqChannel, self.calibIndex, redIndexesEW]
            redundantVisAllLcp = NP.append(redundantVisEW, redundantVisN)
            redundantVisN = self.visRcp[freqChannel, self.calibIndex, redIndexesN]
            redundantVisEW = self.visRcp[freqChannel, self.calibIndex, redIndexesEW]
            redundantVisAllRcp = NP.append(redundantVisEW, redundantVisN)
            
        ampMatrix = NP.abs(self.phaMatrixGenPairsEWN(baselinesNumber, antNumberEW, antNumberN))
        
        allAmp = NP.abs(redundantVisAllLcp)
        antAmp, c, d, e = NP.linalg.lstsq(ampMatrix,NP.log(allAmp), rcond=None)
        antAmp= NP.exp(antAmp[baselinesNumber*2:])
        self.ewAntAmpLcp[freqChannel] = antAmp[:antNumberEW]
        self.nAntAmpLcp[freqChannel] = antAmp[antNumberEW:]
        self.ewAntAmpLcp[freqChannel][self.ewAntAmpLcp[freqChannel]<NP.median(self.ewAntAmpLcp[freqChannel])*0.6] = 1e6
        self.nAntAmpLcp[freqChannel][self.nAntAmpLcp[freqChannel]<NP.median(self.nAntAmpLcp[freqChannel])*0.6] = 1e6
        
        allAmp = NP.abs(redundantVisAllRcp)
        antAmp, c, d, e = NP.linalg.lstsq(ampMatrix,NP.log(allAmp), rcond=None)
        antAmp= NP.exp(antAmp[baselinesNumber*2:])
        self.ewAntAmpRcp[freqChannel] = antAmp[:antNumberEW]
        self.nAntAmpRcp[freqChannel] = antAmp[antNumberEW:]
        self.ewAntAmpRcp[freqChannel][self.ewAntAmpRcp[freqChannel]<NP.median(self.ewAntAmpRcp[freqChannel])*0.6] = 1e6
        self.nAntAmpRcp[freqChannel][self.nAntAmpRcp[freqChannel]<NP.median(self.nAntAmpRcp[freqChannel])*0.6] = 1e6
    
    def eastWestGainsFunc_constrained(self, x, obsVis, antNumber, baselineNumber, freq):
        res = NP.zeros_like(obsVis, dtype = complex)
        x_complex = self.real_to_complex(x)
        solVis = NP.append((NP.exp(1j*self.ewSolarPhase[freq])), x_complex[:baselineNumber-1])
        gains = NP.insert(x_complex[baselineNumber-1:], 32, (1+0j))
        
        solVisArray = NP.array(())
        antAGains = NP.array(())
        antBGains = NP.array(())
        for baseline in range(1, baselineNumber+1):
            solVisArray = NP.append(solVisArray, NP.full(antNumber-baseline, solVis[baseline-1]))
            antAGains = NP.append(antAGains, gains[:antNumber-baseline])
            antBGains = NP.append(antBGains, gains[baseline:])
        res = solVisArray * antAGains * NP.conj(antBGains) - obsVis
        return self.complex_to_real(res)

    def northGainsFunc_constrained(self, x, obsVis, antNumber, baselineNumber, freq):
        antNumber_c = antNumber + 1
        res = NP.zeros_like(obsVis, dtype = complex)
        x_complex = self.real_to_complex(x)
        solVis = NP.append((NP.exp(1j*self.nSolarPhase[freq])), x_complex[:baselineNumber-1])
        gains = NP.append((1+0j), x_complex[baselineNumber-1:])
        
        solVisArray = NP.array(())
        antAGains = NP.array(())
        antBGains = NP.array(())
        for baseline in range(1, baselineNumber+1):
            solVisArray = NP.append(solVisArray, NP.full(antNumber_c-baseline, solVis[baseline-1]))
            antAGains = NP.append(antAGains, gains[:antNumber_c-baseline])
            antBGains = NP.append(antBGains, gains[baseline:])
        res = solVisArray * antAGains * NP.conj(antBGains) - obsVis
        return self.complex_to_real(res)
    
    def allGainsFunc_constrained(self, x, obsVis, ewAntNumber, sAntNumber, baselineNumber, freq):
        res = NP.zeros_like(obsVis, dtype = complex)
        ewSolarAmp = 1
        sSolarAmp = NP.abs(x[0])
        x_complex = self.real_to_complex(x[1:])
        
        sAntNumber_c = sAntNumber + 1
        
        sGainsNumber = sAntNumber
        ewGainsNumber = ewAntNumber
        sSolVisNumber = baselineNumber - 1
        ewSolVisNumber = baselineNumber - 1
        ewSolVis = NP.append((ewSolarAmp * NP.exp(1j*self.ewSolarPhase[freq])), x_complex[: ewSolVisNumber])
        sSolVis = NP.append((sSolarAmp * NP.exp(1j*self.sSolarPhase[freq])), x_complex[ewSolVisNumber : ewSolVisNumber+sSolVisNumber])
        ewGains = x_complex[ewSolVisNumber+sSolVisNumber : ewSolVisNumber+sSolVisNumber+ewGainsNumber]
        sGains = NP.append(ewGains[69], x_complex[ewSolVisNumber+sSolVisNumber+ewGainsNumber :])
        
        solVisArrayS = NP.array(())
        antAGainsS = NP.array(())
        antBGainsS = NP.array(())
        solVisArrayEW = NP.array(())
        antAGainsEW = NP.array(())
        antBGainsEW = NP.array(())
        for baseline in range(1, baselineNumber+1):
            solVisArrayS = NP.append(solVisArrayS, NP.full(self.redIndexesS_len[baseline-1], sSolVis[baseline-1]))
            solVisArrayEW = NP.append(solVisArrayEW, NP.full(self.redIndexesEW_len[baseline-1], ewSolVis[baseline-1]))
            
        antA = self.antennaA[self.redIndexesS.astype(int)] - 138
        antB = self.antennaB[self.redIndexesS.astype(int)] - 138
        ind2swap = NP.where(antB < 0)[0]
        antB[ind2swap] = 0
        antA[ind2swap], antB[ind2swap] = antB[ind2swap], antA[ind2swap]
        antAGainsS = NP.append(antAGainsS, sGains[antA])
        antBGainsS = NP.append(antBGainsS, sGains[antB])
        
        antAGainsEW = NP.append(antAGainsEW, ewGains[self.antennaA[self.redIndexesEW.astype(int)]])
        antBGainsEW = NP.append(antBGainsEW, ewGains[self.antennaB[self.redIndexesEW.astype(int)]])
            
        res = NP.append(solVisArrayEW, solVisArrayS) * NP.append(antAGainsEW, antAGainsS) * NP.conj(NP.append(antBGainsEW, antBGainsS)) - obsVis
        return self.complex_to_real(res)  
    
    
    def buildEwPhase(self):
        newLcpPhaseCorrection = NP.zeros(self.antNumberEW)
        newRcpPhaseCorrection = NP.zeros(self.antNumberEW)
        for j in range(self.antNumberEW):
                newLcpPhaseCorrection[j] += NP.deg2rad(self.ewSlopeLcp[self.frequencyChannel] * (j - 32)) 
                newRcpPhaseCorrection[j] += NP.deg2rad(self.ewSlopeRcp[self.frequencyChannel] * (j - 32))
        self.ewLcpPhaseCorrection[self.frequencyChannel, :] = newLcpPhaseCorrection[:]
        self.ewRcpPhaseCorrection[self.frequencyChannel, :] = newRcpPhaseCorrection[:]
        
    def buildSPhase(self):
        newLcpPhaseCorrection = NP.zeros(self.antNumberS)
        newRcpPhaseCorrection = NP.zeros(self.antNumberS)
        for j in range(self.antNumberS):
                newLcpPhaseCorrection[j] += NP.deg2rad(-self.sSlopeLcp[self.frequencyChannel] * (j + 1)) 
                newRcpPhaseCorrection[j] += NP.deg2rad(-self.sSlopeRcp[self.frequencyChannel] * (j + 1))
        self.sLcpPhaseCorrection[self.frequencyChannel, :] = newLcpPhaseCorrection[:]
        self.sRcpPhaseCorrection[self.frequencyChannel, :] = newRcpPhaseCorrection[:]
    
    def real_to_complex(self, z):
        return z[:len(z)//2] + 1j * z[len(z)//2:]
    
    def complex_to_real(self, z):
        return NP.concatenate((NP.real(z), NP.imag(z)))
    
    def setCalibIndex(self, calibIndex):
        self.calibIndex = calibIndex;

    def setFrequencyChannel(self, channel):
        self.frequencyChannel = channel
        
    def vis2uv(self, scan, phaseCorrect = True, amplitudeCorrect = False, PSF=False, average = 0):

        antX = (self.antX/self.base).astype(int)
        antY = (self.antY/self.base).astype(int)
        
        self.uvLcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex)
        self.uvRcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex)
        
        if average:
            firstScan = scan
            if  self.visLcp.shape[1] < (scan + average):
                lastScan = self.dataLength
            else:
                lastScan = scan + average
        
        O = self.sizeOfUv//2
        for i in range(self.antNumberS):
            for j in range(self.antNumberEW):
                vis = i*self.antNumberEW + j
                antA = min(self.antennaA[vis], self.antennaB[vis])
                antB = max(self.antennaA[vis], self.antennaB[vis])
                antAInd = NP.where(self.antennaNumbers==str(antA))[0][0]
                antBInd = NP.where(self.antennaNumbers==str(antB))[0][0]
                antAX, antAY = antX[antAInd], antY[antAInd]
                antBX, antBY = antX[antBInd], antY[antBInd]
                u = antAY - antBY
                v = antBX - antAX
                                
                if average:
                    self.uvLcp[O + v*2, O + u*2] = NP.mean(self.visLcp[self.frequencyChannel, firstScan:lastScan, vis])
                    self.uvRcp[O + v*2, O + u*2] = NP.mean(self.visRcp[self.frequencyChannel, firstScan:lastScan, vis])
                else:
                    self.uvLcp[O + v*2, O + u*2] = self.visLcp[self.frequencyChannel, scan, vis]
                    self.uvRcp[O + v*2, O + u*2] = self.visRcp[self.frequencyChannel, scan, vis]
                
                if (phaseCorrect):
                    ewPh = self.ewAntPhaLcp[self.frequencyChannel, j]+self.ewLcpPhaseCorrection[self.frequencyChannel, j]
                    sPh = self.sAntPhaLcp[self.frequencyChannel, i]+self.sLcpPhaseCorrection[self.frequencyChannel, i]
                    self.uvLcp[O + v*2, O + u*2] *= NP.exp(1j * (-ewPh + sPh))
                    ewPh = self.ewAntPhaRcp[self.frequencyChannel, j]+self.ewRcpPhaseCorrection[self.frequencyChannel, j]
                    sPh = self.sAntPhaRcp[self.frequencyChannel, i]+self.sRcpPhaseCorrection[self.frequencyChannel, i]
                    self.uvRcp[O + v*2, O + u*2] *= NP.exp(1j * (-ewPh + sPh))
                if (amplitudeCorrect):
                    self.uvLcp[O + v*2, O + u*2] /= (self.ewAntAmpLcp[self.frequencyChannel, j] * self.sAntAmpLcp[self.frequencyChannel, i])
                    self.uvRcp[O + v*2, O + u*2] /= (self.ewAntAmpRcp[self.frequencyChannel, j] * self.sAntAmpRcp[self.frequencyChannel, i])
                
                self.uvLcp[O - v*2, O - u*2] = NP.conj(self.uvLcp[O + v*2, O + u*2])
                self.uvRcp[O - v*2, O - u*2] = NP.conj(self.uvRcp[O + v*2, O + u*2])
                
                
        for i in range(len(self.antZeroRow)):
            vis = self.antZeroRow[i]
            if i<69:
                antA = self.antennaA[vis]
                antB = self.antennaB[vis]
                antAInd = NP.where(self.antennaNumbers==str(antA))[0][0]
                antBInd = NP.where(self.antennaNumbers==str(antB))[0][0]
                antAX, antAY = antX[antAInd], antY[antAInd]
                antBX, antBY = antX[antBInd], antY[antBInd]
                u = antAY - antBY
                v = antBX - antAX
                
                if average:
                    self.uvLcp[O + v*2, O + u*2] = NP.mean(self.visLcp[self.frequencyChannel, firstScan:lastScan, vis])
                    self.uvRcp[O + v*2, O + u*2] = NP.mean(self.visRcp[self.frequencyChannel, firstScan:lastScan, vis])
                else:
                    self.uvLcp[O + v*2, O + u*2] = self.visLcp[self.frequencyChannel, scan, vis]
                    self.uvRcp[O + v*2, O + u*2] = self.visRcp[self.frequencyChannel, scan, vis]
                
                if (phaseCorrect):
                    ewPh1 = self.ewAntPhaLcp[self.frequencyChannel, i]+self.ewLcpPhaseCorrection[self.frequencyChannel, i]
                    ewPh2 = self.ewAntPhaLcp[self.frequencyChannel, 69]+self.ewLcpPhaseCorrection[self.frequencyChannel, 69]
                    self.uvLcp[O + v*2, O + u*2] *= NP.exp(1j * (-ewPh1 + ewPh2))
                    ewPh1 = self.ewAntPhaRcp[self.frequencyChannel, i]+self.ewRcpPhaseCorrection[self.frequencyChannel, i]
                    ewPh2 = self.ewAntPhaRcp[self.frequencyChannel, 69]+self.ewRcpPhaseCorrection[self.frequencyChannel, 69]
                    self.uvRcp[O + v*2, O + u*2] *= NP.exp(1j * (-ewPh1 + ewPh2))
                if (amplitudeCorrect):
                    self.uvLcp[O + v*2, O + u*2] /= (self.ewAntAmpLcp[self.frequencyChannel, i] * self.ewAntAmpLcp[self.frequencyChannel, 69])
                    self.uvRcp[O + v*2, O + u*2] /= (self.ewAntAmpRcp[self.frequencyChannel, i] * self.ewAntAmpRcp[self.frequencyChannel, 69])
                
            else:
                antA = self.antennaB[vis]
                antB = self.antennaA[vis]
                antAInd = NP.where(self.antennaNumbers==str(antA))[0][0]
                antBInd = NP.where(self.antennaNumbers==str(antB))[0][0]
                antAX, antAY = antX[antAInd], antY[antAInd]
                antBX, antBY = antX[antBInd], antY[antBInd]
                u = antAY - antBY
                v = antBX - antAX
                
                if average:
                    self.uvLcp[O + v*2, O + u*2] =  NP.conj(NP.mean(self.visLcp[self.frequencyChannel, firstScan:lastScan, vis]))
                    self.uvRcp[O + v*2, O + u*2] =  NP.conj(NP.mean(self.visRcp[self.frequencyChannel, firstScan:lastScan, vis]))
                else:
                    self.uvLcp[O + v*2, O + u*2] = NP.conj(self.visLcp[self.frequencyChannel, scan, vis])
                    self.uvRcp[O + v*2, O + u*2] = NP.conj(self.visRcp[self.frequencyChannel, scan, vis])
                
                if (phaseCorrect):
                    ewPh1 = self.ewAntPhaLcp[self.frequencyChannel, i]+self.ewLcpPhaseCorrection[self.frequencyChannel, i]
                    ewPh2 = self.ewAntPhaLcp[self.frequencyChannel, 69]+self.ewLcpPhaseCorrection[self.frequencyChannel, 69]
                    self.uvLcp[O + v*2, O + u*2] *= NP.exp(1j * (-ewPh1 + ewPh2))
                    ewPh1 = self.ewAntPhaRcp[self.frequencyChannel, i]+self.ewRcpPhaseCorrection[self.frequencyChannel, i]
                    ewPh2 = self.ewAntPhaRcp[self.frequencyChannel, 69]+self.ewRcpPhaseCorrection[self.frequencyChannel, 69]
                    self.uvRcp[O + v*2, O + u*2] *= NP.exp(1j * (-ewPh1 + ewPh2))
                if (amplitudeCorrect):
                    self.uvLcp[O + v*2, O + u*2] /= (self.ewAntAmpLcp[self.frequencyChannel, i] * self.ewAntAmpLcp[self.frequencyChannel, 69])
                    self.uvRcp[O + v*2, O + u*2] /= (self.ewAntAmpRcp[self.frequencyChannel, i] * self.ewAntAmpRcp[self.frequencyChannel, 69])
                
       
        self.uvLcp[O,O] = self.lcpShift[self.frequencyChannel]
        self.uvRcp[O,O] = self.rcpShift[self.frequencyChannel]
        
        if PSF:
            self.uvLcp[NP.abs(self.uvLcp)>1e-8] = 1
            self.uvRcp[NP.abs(self.uvRcp)>1e-8] = 1
            
        self.uvLcp[NP.abs(self.uvLcp)<1e-6] = 0.
        self.uvRcp[NP.abs(self.uvRcp)<1e-6] = 0.
        self.uvLcp /= NP.count_nonzero(self.uvLcp)
        self.uvRcp /= NP.count_nonzero(self.uvRcp)

    def uv2lmImage(self):
        self.lcp = NP.fft.fft2(NP.roll(NP.roll(self.uvLcp,self.sizeOfUv//2+1,0),self.sizeOfUv//2+1,1));
        self.lcp = NP.roll(NP.roll(self.lcp,self.sizeOfUv//2-1,0),self.sizeOfUv//2-1,1);
        self.lcp = NP.flip(self.lcp, 1)
        self.rcp = NP.fft.fft2(NP.roll(NP.roll(self.uvRcp,self.sizeOfUv//2+1,0),self.sizeOfUv//2+1,1));
        self.rcp = NP.roll(NP.roll(self.rcp,self.sizeOfUv//2-1,0),self.sizeOfUv//2-1,1);
        self.rcp = NP.flip(self.rcp, 1)
        
    def lm2Heliocentric(self):
        scaling = self.getPQScale(self.sizeOfUv, NP.deg2rad(self.arcsecPerPixel * (self.sizeOfUv - 1)/3600.)*2)
        scale = AffineTransform(scale=(self.sizeOfUv/scaling[0], self.sizeOfUv/scaling[1]))
        shift = AffineTransform(translation=(-self.sizeOfUv/2,-self.sizeOfUv/2))
        rotate = AffineTransform(rotation = self.pAngle)
        matrix = AffineTransform(matrix = self.getPQ2HDMatrix())
        back_shift = AffineTransform(translation=(self.sizeOfUv/2,self.sizeOfUv/2))

        O = self.sizeOfUv//2
        Q = self.sizeOfUv//4
        dataResult0 = warp(self.lcp.real,(shift + (scale + back_shift)).inverse)
        self.lcp = warp(dataResult0,(shift + (matrix + back_shift)).inverse)
        dataResult0 = warp(self.rcp.real,(shift + (scale + back_shift)).inverse)
        self.rcp = warp(dataResult0,(shift + (matrix + back_shift)).inverse)
        dataResult0 = 0
        self.lcp = warp(self.lcp,(shift + (rotate + back_shift)).inverse)[O-Q:O+Q,O-Q:O+Q]
        self.rcp = warp(self.rcp,(shift + (rotate + back_shift)).inverse)[O-Q:O+Q,O-Q:O+Q]

    def calculateDelays(self):
        self.nDelays = NP.zeros((self.antNumberN, self.dataLength))
        self.ewDelays = NP.zeros((self.antNumberEW, self.dataLength))
        base = 9.8
        ew_center = 32
        ns_center = 0
        _PHI = 0.903338787600965
        for scan in range(self.dataLength):
            hourAngle = self.omegaEarth * (self.freqTime[self.frequencyChannel, scan] - self.RAO.culmination)
            dec = self.getDeclination()
            for ant in range(self.antNumberN):
                cosQ = NP.cos(hourAngle) * NP.cos(dec)*NP.sin(_PHI) - NP.sin(dec)*NP.cos(_PHI);
                M = ns_center - ant - 1;
                self.nDelays[ant, scan] = base * M * cosQ / constants.c.to_value();
            for ant in range(self.antNumberEW):
                cosP = NP.sin(hourAngle) * NP.cos(dec);
                N = ew_center - ant;
                self.ewDelays[ant, scan] = base * N * cosP / constants.c.to_value();
                
    def resetDelays(self):
        self.nDelays = NP.zeros((self.antNumberN, self.dataLength))
        self.ewDelays = NP.zeros((self.antNumberEW, self.dataLength))
        
    def wrap(self, value):
        while value<-180:
            value+=360
        while value>180:
            value-=360
        return value
    
    def createDisk(self, radius, arcsecPerPixel = 2.45552):
        qSun = NP.zeros((self.sizeOfUv, self.sizeOfUv))
        sunRadius = radius / (arcsecPerPixel*2)
        for i in range(self.sizeOfUv):
            x = i - self.sizeOfUv//2 - 1
            for j in range(self.sizeOfUv):
                y = j - self.sizeOfUv//2 - 1
                if (NP.sqrt(x*x + y*y) < sunRadius):
                    qSun[i , j] = 1
                    
        dL = 2*( 12//2) + 1
        arg_x = NP.linspace(-1.,1,dL)
        arg_y = NP.linspace(-1.,1,dL)
        xx, yy = NP.meshgrid(arg_x, arg_y)
        
        scaling = self.getPQScale(self.sizeOfUv, NP.deg2rad(arcsecPerPixel*(self.sizeOfUv-1)/3600.)*2)
        scale = AffineTransform(scale=(scaling[0]/self.sizeOfUv, scaling[1]/self.sizeOfUv))
        back_shift = AffineTransform(translation=(self.sizeOfUv/2, self.sizeOfUv/2))
        shift = AffineTransform(translation=(-self.sizeOfUv/2, -self.sizeOfUv/2))
        matrix = AffineTransform(matrix = NP.linalg.inv(self.getPQ2HDMatrix()))
        rotate = AffineTransform(rotation = -self.pAngle)
        
        gKern =   NP.exp(-0.5*(xx**2 + yy**2))
        qSmoothSun = scipy.signal.fftconvolve(qSun,gKern) / dL**2
        qSmoothSun = qSmoothSun[dL//2:dL//2+self.sizeOfUv,dL//2:dL//2+self.sizeOfUv]
        smoothCoef = qSmoothSun[512, 512]
        qSmoothSun /= smoothCoef
        qSun_el_hd = warp(qSmoothSun,(shift + (rotate + back_shift)).inverse)
        
        res = warp(qSun_el_hd, (shift + (matrix + back_shift)).inverse)
        qSun_lm = warp(res,(shift + (scale + back_shift)).inverse)
        qSun_lm_fft = NP.fft.fft2(NP.roll(NP.roll(qSun_lm,self.sizeOfUv//2,0),self.sizeOfUv//2,1));
        qSun_lm_fft = NP.roll(NP.roll(qSun_lm_fft,self.sizeOfUv//2,0),self.sizeOfUv//2,1) / self.sizeOfUv;
        qSun_lm_fft = NP.flip(qSun_lm_fft, 0)
#        qSun_lm_uv = qSun_lm_fft * uvPsf
#        qSun_lm_conv = NP.fft.fft2(NP.roll(NP.roll(qSun_lm_uv,self.sizeOfUv//2+1,0),self.sizeOfUv//2+1,1));
#        qSun_lm_conv = NP.roll(NP.roll(qSun_lm_conv,self.sizeOfUv//2-1,0),self.sizeOfUv//2-1,1);
#        qSun_lm_conv = NP.flip(NP.flip(qSun_lm_conv, 1), 0)
        self.fftDisk = qSun_lm_fft #qSun_lm_conv, 
    
    def createUvUniform(self):
        self.uvUniform = NP.zeros((self.sizeOfUv, self.sizeOfUv), dtype = complex)
        flags_ew = NP.where(self.ewAntAmpLcp[self.frequencyChannel]==1e6)[0]
        flags_n = NP.where(self.nAntAmpLcp[self.frequencyChannel]==1e6)[0]
        O = self.sizeOfUv//2
        for i in range(self.antNumberN):
            for j in range(self.antNumberEW):
                if not (NP.any(flags_ew == j) or NP.any(flags_n == i)):
                    self.uvUniform[O + (i+1)*2, O + (j-32)*2] = 1
                    self.uvUniform[O - (i+1)*2, O - (j-32)*2] = 1
        for i in range(self.antNumberEW):
            if i != 32:
                if not (NP.any(flags_ew == i) or NP.any(flags_ew == 32)):
                    self.uvUniform[O, O + (i-32)*2] = 1
        self.uvUniform[O, O] = 1
                    
    def createUvPsf(self, T, ewSlope, nSlope, shift):
        self.uvPsf = self.uvUniform.copy()
        O = self.sizeOfUv//2
        ewSlope = NP.deg2rad(ewSlope)
        nSlope = NP.deg2rad(nSlope)
        ewSlopeUv = NP.linspace(-O * ewSlope/2., O * ewSlope/2., self.sizeOfUv)
        nSlopeUv = NP.linspace(-O * nSlope/2., O * nSlope/2., self.sizeOfUv)
        ewGrid,nGrid = NP.meshgrid(ewSlopeUv, nSlopeUv)
        slopeGrid = ewGrid + nGrid
        slopeGrid[self.uvUniform == 0] = 0
        self.uvPsf *= T * NP.exp(1j * slopeGrid)
        self.uvPsf[O,O] = shift
    
    def diskDiff(self, x, pol):
        self.createUvPsf(x[0], x[1], x[2], x[3])
        uvDisk = self.fftDisk * self.uvPsf
        if pol == 0:
            diff = self.uvLcp - uvDisk
        if pol == 1:
            diff = self.uvRcp - uvDisk
        return self.complex_to_real(diff[self.uvUniform!=0])
#        qSun_lm_conv = NP.fft.fft2(NP.roll(NP.roll(diff,uvSize//2+1,0),uvSize//2+1,1));
#        return NP.abs(NP.reshape(qSun_lm_conv, uvSize**2))
    
    def findDisk(self):
        self.createDisk(980)
        self.createUvUniform()
        self.x_ini = [1,0,0,1]
        # x_ini = [1,0,0]
        self.center_ls_res_lcp = least_squares(self.diskDiff, self.x_ini, args = (0,))
        _diskLevelLcp, _ewSlopeLcp, _nSlopeLcp, _shiftLcp = self.center_ls_res_lcp['x']
        self.center_ls_res_rcp = least_squares(self.diskDiff, self.x_ini, args = (1,))
        _diskLevelRcp, _ewSlopeRcp, _nSlopeRcp, _shiftRcp = self.center_ls_res_rcp['x']
        
        self.diskLevelLcp[self.frequencyChannel] = _diskLevelLcp
        self.diskLevelRcp[self.frequencyChannel] = _diskLevelRcp
        
        Tb = self.ZirinQSunTb.getTbAtFrequency(self.freqList[self.frequencyChannel]*1e-6) * 1e3
        
        self.lcpShift[self.frequencyChannel] = self.lcpShift[self.frequencyChannel]/(_shiftLcp * self.convolutionNormCoef / Tb)
        self.rcpShift[self.frequencyChannel] = self.rcpShift[self.frequencyChannel]/(_shiftRcp * self.convolutionNormCoef / Tb)
        
        self.ewAntAmpLcp[self.frequencyChannel][self.ewAntAmpLcp[self.frequencyChannel]!=1e6] *= NP.sqrt(_diskLevelLcp*self.convolutionNormCoef / Tb)
        self.nAntAmpLcp[self.frequencyChannel][self.nAntAmpLcp[self.frequencyChannel]!=1e6] *= NP.sqrt(_diskLevelLcp*self.convolutionNormCoef / Tb)
        self.ewAntAmpRcp[self.frequencyChannel][self.ewAntAmpRcp[self.frequencyChannel]!=1e6] *= NP.sqrt(_diskLevelRcp*self.convolutionNormCoef / Tb)
        self.nAntAmpRcp[self.frequencyChannel][self.nAntAmpRcp[self.frequencyChannel]!=1e6] *= NP.sqrt(_diskLevelRcp*self.convolutionNormCoef / Tb)
        
        self.ewSlopeLcp[self.frequencyChannel] = self.wrap(self.ewSlopeLcp[self.frequencyChannel] + _ewSlopeLcp)
        self.nSlopeLcp[self.frequencyChannel] = self.wrap(self.nSlopeLcp[self.frequencyChannel] + _nSlopeLcp)
        self.ewSlopeRcp[self.frequencyChannel] = self.wrap(self.ewSlopeRcp[self.frequencyChannel] + _ewSlopeRcp)
        self.nSlopeRcp[self.frequencyChannel] = self.wrap(self.nSlopeRcp[self.frequencyChannel] + _nSlopeRcp)
          
        
    def diskDiff_2(self, x, pol):
        self.createUvPsf(x[0], x[1], x[2])
        uvDisk = self.fftDisk * self.uvPsf
        if pol == 0:
            diff = self.uvLcp - uvDisk
        if pol == 1:
            diff = self.uvRcp - uvDisk
        return NP.sum(self.complex_to_real(diff[self.uvUniform!=0])**2)
#        qSun_lm_conv = NP.fft.fft2(NP.roll(NP.roll(diff,uvSize//2+1,0),uvSize//2+1,1));
#        return NP.abs(NP.reshape(qSun_lm_conv, uvSize**2))
    
    def findDisk_2(self):
        self.createDisk(980)
        self.createUvUniform()
        self.x_ini = [1,0,0]
        ls_res = basinhopping(self.diskDiff_2, self.x_ini, stepsize=180, minimizer_kwargs = {'args':(0,)})
        self.diskLevelLcp, self.ewSlopeLcp, self.nSlopeLcp = ls_res['x']
        ls_res = basinhopping(self.diskDiff_2, self.x_ini, stepsize=180, minimizer_kwargs = {'args':(1,)})
        self.diskLevelRcp, self.ewSlopeRcp, self.nSlopeRcp = ls_res['x']
        
    def findDisk_3(self):
        self.createDisk(980)
        self.createUvUniform()
        fun_lcp = 10
        fun_rcp = 10
        for i in range(3):
            for j in range(3):
                self.x_ini = [1, -90+i*90, -90+j*90]
                ls_res = least_squares(self.diskDiff, self.x_ini, args = (0,))
                print(NP.sum(ls_res['fun']**2))
                if NP.sum(ls_res['fun']**2)<fun_lcp:
                    print('min updated')
                    self.diskLevelLcp, self.ewSlopeLcp, self.nSlopeLcp = ls_res['x']
                    fun_lcp = NP.sum(ls_res['fun']**2)
                ls_res = least_squares(self.diskDiff, self.x_ini, args = (1,))
                if NP.sum(ls_res['fun']**2)<fun_rcp:
                    self.diskLevelRcp, self.ewSlopeRcp, self.nSlopeRcp = ls_res['x']
                    fun_rcp = NP.sum(ls_res['fun']**2)
        while self.ewSlopeLcp<-180:
            self.ewSlopeLcp+=360
        while self.ewSlopeLcp>180:
            self.ewSlopeLcp-=360
        while self.nSlopeLcp<-180:
            self.nSlopeLcp+=360
        while self.nSlopeLcp>180:
            self.nSlopeLcp-=360
        while self.ewSlopeRcp<-180:
            self.ewSlopeRcp+=360
        while self.ewSlopeRcp>180:
            self.ewSlopeRcp-=360
        while self.nSlopeRcp<-180:
            self.nSlopeRcp+=360
        while self.nSlopeRcp>180:
            self.nSlopeRcp-=360
        
    def findDisk_4(self):
        start_time_all = time.time()
        self.createDisk(980)
        self.createUvUniform()
        fun_lcp = 10
        fun_rcp = 10
        _diskLevelLcp = self.diskLevelLcp[self.frequencyChannel].copy()
        _diskLevelRcp = self.diskLevelRcp[self.frequencyChannel].copy()
                
        for i in range(3):
            for j in range(3):
                start_time = time.time()
                
                self.x_ini = [_diskLevelLcp, -90+i*90, -90+j*90]
                
                ls_res = least_squares(self.diskDiff, self.x_ini, args = (0,), ftol=self.centering_ftol)
                print(NP.sum(ls_res['fun']**2))
                if i==0 and j==0:
                    _diskLevelLcp, _ewSlopeLcp, _nSlopeLcp = ls_res['x']
                    fun_lcp = NP.sum(ls_res['fun']**2)
                    _diskLevelRcp, _ewSlopeRcp, _nSlopeRcp = ls_res['x']
                    fun_rcp = NP.sum(ls_res['fun']**2)
                    
                else:
                    if NP.sum(ls_res['fun']**2)<fun_lcp and ls_res['x'][0]>0:
                        print('min updated')
                        _diskLevelLcp, _ewSlopeLcp, _nSlopeLcp = ls_res['x']
                        fun_lcp = NP.sum(ls_res['fun']**2)
                        print((_diskLevelLcp, _ewSlopeLcp, _nSlopeLcp))
 
                    self.x_ini = [_diskLevelRcp, -90+i*90, -90+j*90]
                    ls_res = least_squares(self.diskDiff, self.x_ini, args = (1,), ftol=self.centering_ftol)
                    if NP.sum(ls_res['fun']**2)<fun_rcp and ls_res['x'][0]>0:
                        _diskLevelRcp, _ewSlopeRcp, _nSlopeRcp = ls_res['x']
                        fun_rcp = NP.sum(ls_res['fun']**2)
                    
                print("ITER " + str(i*3+j) + " --- %s seconds ---" % (time.time() - start_time))
                
        self.x_ini = [_diskLevelLcp, _ewSlopeLcp, _nSlopeLcp]               
        ls_res = least_squares(self.diskDiff, self.x_ini, args = (0,), ftol=1e-10)
        _diskLevelLcp, _ewSlopeLcp, _nSlopeLcp = ls_res['x']
        
        self.x_ini = [_diskLevelRcp, _ewSlopeRcp, _nSlopeRcp]               
        ls_res = least_squares(self.diskDiff, self.x_ini, args = (0,), ftol=1e-10)
        _diskLevelRcp, _ewSlopeRcp, _nSlopeRcp = ls_res['x']
        
        self.diskLevelLcp[self.frequencyChannel] = _diskLevelLcp
        self.diskLevelRcp[self.frequencyChannel] = _diskLevelRcp
        
        Tb = self.ZirinQSunTb.getTbAtFrequency(self.freqList[self.frequencyChannel]*1e-6) * 1e3
        
        self.ewAntAmpLcp[self.frequencyChannel][self.ewAntAmpLcp[self.frequencyChannel]!=1e6] *= NP.sqrt(_diskLevelLcp*self.convolutionNormCoef / Tb)
        self.nAntAmpLcp[self.frequencyChannel][self.nAntAmpLcp[self.frequencyChannel]!=1e6] *= NP.sqrt(_diskLevelLcp*self.convolutionNormCoef / Tb)
        self.ewAntAmpRcp[self.frequencyChannel][self.ewAntAmpRcp[self.frequencyChannel]!=1e6] *= NP.sqrt(_diskLevelRcp*self.convolutionNormCoef / Tb)
        self.nAntAmpRcp[self.frequencyChannel][self.nAntAmpRcp[self.frequencyChannel]!=1e6] *= NP.sqrt(_diskLevelRcp*self.convolutionNormCoef / Tb)
        
        self.ewSlopeLcp[self.frequencyChannel] = self.wrap(self.ewSlopeLcp[self.frequencyChannel] + _ewSlopeLcp)
        self.nSlopeLcp[self.frequencyChannel] = self.wrap(self.nSlopeLcp[self.frequencyChannel] + _nSlopeLcp)
        self.ewSlopeRcp[self.frequencyChannel] = self.wrap(self.ewSlopeRcp[self.frequencyChannel] + _ewSlopeRcp)
        self.nSlopeRcp[self.frequencyChannel] = self.wrap(self.nSlopeRcp[self.frequencyChannel] + _nSlopeRcp)
           
        print("FUNCTION RUNTIME  --- %s seconds ---" % (time.time() - start_time_all))
    
    def centerDisk(self):
        self.findDisk()
        self.buildEwPhase()
        self.buildNPhase()
        
    def modelDiskConv(self):
        # self.createUvPsf(self.diskLevelLcp,0,0,0)
        currentDiskTb = self.ZirinQSunTb.getTbAtFrequency(self.freqList[self.frequencyChannel]*1e-6)*1e3
        self.createUvPsf(currentDiskTb/self.convolutionNormCoef,0,0,currentDiskTb/self.convolutionNormCoef)
        self.uvDiskConv = self.fftDisk * self.uvPsf# - self.uvLcp
        qSun_lm = NP.fft.fft2(NP.roll(NP.roll(self.uvDiskConv,self.sizeOfUv//2+1,0),self.sizeOfUv//2+1,1));
        qSun_lm = NP.roll(NP.roll(qSun_lm,self.sizeOfUv//2-1,0),self.sizeOfUv//2-1,1)# / self.sizeOfUv;
        qSun_lm = NP.flip(qSun_lm, 0)
        self.modelDisk = qSun_lm
        
    def modelDiskConv_unity(self):
        self.createDisk(980)
        self.createUvUniform()
        self.createUvPsf(1,0,0,1)
        self.uvDiskConv = self.fftDisk * self.uvPsf# - self.uvLcp
        qSun_lm = NP.fft.fft2(NP.roll(NP.roll(self.uvDiskConv,self.sizeOfUv//2+1,0),self.sizeOfUv//2+1,1));
        qSun_lm = NP.roll(NP.roll(qSun_lm,self.sizeOfUv//2-1,0),self.sizeOfUv//2-1,1)# / self.sizeOfUv;
        qSun_lm = NP.flip(qSun_lm, 0)
        self.modelDisk = qSun_lm