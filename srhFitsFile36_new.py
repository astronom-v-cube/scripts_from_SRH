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
from scipy.optimize import least_squares
import base2uvw_36

class SrhFitsFile():
    def __init__(self, name, sizeOfUv):
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
        self.antNumberEW = 65
        self.antNumberN = 31
        self.averageCalib = False
        self.useNonlinearApproach = True
        
        
        self.badAntsLcp = NP.zeros(128)
        self.badAntsRcp = NP.zeros(128)
        self.sizeOfUv = sizeOfUv
        self.baselines = 3
        self.open(name)
        
        self.flagsIndexes = []
        
                                    
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
            
            self.antZeroRow = self.hduList[3].data['ant_zero_row'][:97]
            self.RAO = BadaryRAO(self.dateObs.split('T')[0])
            
            self.ewAntPhaLcp = NP.zeros((self.freqListLength, self.antNumberEW))
            self.nAntPhaLcp = NP.zeros((self.freqListLength, self.antNumberN))
            self.ewAntPhaRcp = NP.zeros((self.freqListLength, self.antNumberEW))
            self.nAntPhaRcp = NP.zeros((self.freqListLength, self.antNumberN))
            self.ewLcpPhaseCorrection = NP.zeros((self.freqListLength, self.antNumberEW))
            self.ewRcpPhaseCorrection = NP.zeros((self.freqListLength, self.antNumberEW))
            self.nLcpPhaseCorrection = NP.zeros((self.freqListLength, self.antNumberN))
            self.nRcpPhaseCorrection = NP.zeros((self.freqListLength, self.antNumberN))
            self.nSolarPhase = NP.zeros(self.freqListLength)
            self.ewSolarPhase = NP.zeros(self.freqListLength)

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

    def getHourAngle(self, scan):
        self.hAngle = self.omegaEarth * (self.freqTime[self.frequencyChannel, scan] - self.RAO.culmination)
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

    def phaMatrixGenPairsEWN(self, pairs, antNumberEW, antNumberN):
        rowsEW = int(((antNumberEW - 1) + (antNumberEW - pairs))/2 * pairs)
        rowsN = int(((antNumberN) + (antNumberN + 1 - pairs))/2 * pairs)
        colsEW = antNumberEW + pairs
        colsN = antNumberN + pairs
        phaMatrix = NP.zeros((rowsEW + rowsN + 2, colsEW + colsN))
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
                    phaMatrix[rowsEW + row0, 2*pairs + antNumberEW//2] = 1
                else:
                    phaMatrix[phaPair + rowsEW + row0, 2*pairs + phaPair + antNumberEW - 1] = 1
                phaMatrix[phaPair + rowsEW + row0, 2*pairs + phaPair + antNumberEW - 1 + (pair + 1)] = -1
        phaMatrix[rowsEW + rowsN, 0] = 1
        phaMatrix[rowsEW + rowsN+1, pairs] = 1
        return phaMatrix.copy()
    
    def calculatePhaseCalibration(self, baselinesNumber = 1):
       for freq in range(self.freqListLength):
           self.solarPhase(freq)
           self.updateAntennaPhase(freq, baselinesNumber)

    def updateAntennaPhase(self, freqChannel, baselinesNumber = 1):
        if self.useNonlinearApproach:
            self.calculatePhaseLcp_nonlinear(freqChannel, baselinesNumber = baselinesNumber)
#            self.calculatePhaseRcp_nonlinear(freqChannel, baselinesNumber = baselinesNumber)
        else:
            self.calculatePhase_linear(freqChannel, baselinesNumber = baselinesNumber)
            
    def solarPhase(self, freq):
        u,v,w = base2uvw_36.base2uvw(self.hAngle, self.RAO.declination, 98, 99)
        baseWave = NP.sqrt(u**2+v**2)*self.freqList[freq]*1e3/constants.c.to_value()
        if baseWave > 120:
            self.nSolarPhase[freq] = NP.pi
        else:
            self.nSolarPhase[freq] = 0
        u,v,w = base2uvw_36.base2uvw(self.hAngle, self.RAO.declination, 1, 2)
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
        
    def calculatePhaseLcp_nonlinear(self, freqChannel, baselinesNumber = 2):
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
             
        redundantVisN = self.visLcp[freqChannel, self.calibIndex, redIndexesN]
        redundantVisEW = self.visLcp[freqChannel, self.calibIndex, redIndexesEW]
        redundantVisAll = NP.append(redundantVisEW, redundantVisN)
        
        x_size = (baselinesNumber-1)*2 + antNumberEW + antNumberN
        x_ini = NP.concatenate((NP.ones(x_size), NP.zeros(x_size)))
        ls_res = least_squares(self.allGainsFunc_constrained, x_ini, args = (redundantVisAll, antNumberEW, antNumberN, baselinesNumber, freqChannel), max_nfev = 150)
        
        gains = self.real_to_complex(ls_res['x'])[(baselinesNumber-1)*2:]
        self.ew_gains_lcp = gains[:antNumberEW]
        self.ewAntPhaLcp[freqChannel] = NP.angle(self.ew_gains_lcp)
        self.n_gains_lcp = gains[antNumberEW:]
        self.nAntPhaLcp[freqChannel] = NP.angle(self.n_gains_lcp)
        
    def calculatePhaseRcp_nonlinear(self, freqChannel, baselinesNumber = 2):
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
             
        redundantVisN = self.visRcp[freqChannel, self.calibIndex, redIndexesN]
        redundantVisEW = self.visRcp[freqChannel, self.calibIndex, redIndexesEW]
        redundantVisAll = NP.append(redundantVisEW, redundantVisN)
        
        x_size = (baselinesNumber-1)*2 + antNumberEW + antNumberN
        x_ini = NP.concatenate((NP.ones(x_size), NP.zeros(x_size)))
        ls_res = least_squares(self.allGainsFunc_constrained, x_ini, args = (redundantVisAll, antNumberEW, antNumberN, baselinesNumber, freqChannel), max_nfev = 200)
        
        gains = self.real_to_complex(ls_res['x'])[(baselinesNumber-1)*2:]
        self.ew_gains_rcp = gains[:antNumberEW]
        self.ewAntPhaRcp[freqChannel] = NP.angle(self.ew_gains_rcp)
        self.n_gains_rcp = gains[antNumberEW:]
        self.nAntPhaRcp[freqChannel] = NP.angle(self.n_gains_rcp)
                
    def eastWestGainsFunc_constrained(self, x, obsVis, antNumber, baselineNumber):
        res = NP.zeros_like(obsVis, dtype = complex)
        x_complex = self.real_to_complex(x)
        solVis = NP.append((1+0j), x_complex[:baselineNumber-1])
        gains = NP.insert(x_complex[baselineNumber-1:], 32, (1+0j))
        i = 0
        for baseline in range(1, baselineNumber+1):
            for antA in range(antNumber-baseline):
                antB = antA + baseline
                res[i] = solVis[baseline-1] * gains[antA] * NP.conj(gains[antB]) - obsVis[i]
                i+=1
        return self.complex_to_real(res)

    def northGainsFunc_constrained(self, x, obsVis, antNumber, baselineNumber):
        res = NP.zeros_like(obsVis, dtype = complex)
        x_complex = self.real_to_complex(x)
        solVis = NP.append((1+0j), x_complex[:baselineNumber-1])
        gains = NP.append((1+0j), x_complex[baselineNumber-1:])
        i = 0
        for baseline in range(1, baselineNumber+1):
            for antA in range(antNumber-baseline):
                antB = antA + baseline
                res[i] = solVis[baseline-1] * gains[antA] * NP.conj(gains[antB]) - obsVis[i]
                i+=1
        return self.complex_to_real(res)    
    
    def allGainsFunc_constrained(self, x, obsVis, ewAntNumber, nAntNumber, baselineNumber, freq):
        res = NP.zeros_like(obsVis, dtype = complex)
        x_complex = self.real_to_complex(x)
        
        nGainsNumber = nAntNumber
        ewGainsNumber = ewAntNumber
        nSolVisNumber = baselineNumber - 1
        ewSolVisNumber = baselineNumber - 1
        ewSolVis = NP.append((NP.exp(1j*self.ewSolarPhase[freq])), x_complex[: ewSolVisNumber])
        nSolVis = NP.append((NP.exp(1j*self.nSolarPhase[freq])), x_complex[ewSolVisNumber : ewSolVisNumber+nSolVisNumber])
        ewGains = x_complex[ewSolVisNumber+nSolVisNumber : ewSolVisNumber+nSolVisNumber+ewGainsNumber]
        nGains = NP.append(ewGains[32], x_complex[ewSolVisNumber+nSolVisNumber+ewGainsNumber :])
        i = 0
        for baseline in range(1, baselineNumber+1):
            for antA in range(ewAntNumber-baseline):
                antB = antA + baseline
                res[i] = ewSolVis[baseline-1] * ewGains[antA] * NP.conj(ewGains[antB]) - obsVis[i]
                i+=1
        for baseline in range(1, baselineNumber+1):
            for antA in range(nAntNumber-baseline):
                antB = antA + baseline
                res[i] = nSolVis[baseline-1] * nGains[antA] * NP.conj(nGains[antB]) - obsVis[i]
                i+=1
        return self.complex_to_real(res)  
    
    def changeEastWestPhase(self, newLcpPhaseCorrection, newRcpPhaseCorrection):
        self.ewLcpPhaseCorrection[self.frequencyChannel, :] = newLcpPhaseCorrection[:]
        self.ewRcpPhaseCorrection[self.frequencyChannel, :] = newRcpPhaseCorrection[:]
        
    def changeNorthPhase(self, newLcpPhaseCorrection, newRcpPhaseCorrection):
        self.nLcpPhaseCorrection[self.frequencyChannel, :] = newLcpPhaseCorrection[:]
        self.nRcpPhaseCorrection[self.frequencyChannel, :] = newRcpPhaseCorrection[:]
    
    def real_to_complex(self, z):
        return z[:len(z)//2] + 1j * z[len(z)//2:]
    
    def complex_to_real(self, z):
        return NP.concatenate((NP.real(z), NP.imag(z)))
    
    def setCalibIndex(self, calibIndex):
        self.calibIndex = calibIndex;

    def setFrequencyChannel(self, channel):
        self.frequencyChannel = channel
        
    def vis2uv(self, scan, phaseCorrect = True, amplitudeCorrect = False, PSF=False, average = 0):
        self.uvLcp[:,:] = complex(0,0)
        self.uvRcp[:,:] = complex(0,0)
        O = self.sizeOfUv//2
        if average:
            for i in range(31):
                for j in range(self.antNumberEW):
                    self.uvLcp[O + (i+1)*2, O + (j-32)*2] = NP.mean(self.visLcp[self.frequencyChannel, :, i*97+j])
                    self.uvRcp[O + (i+1)*2, O + (j-32)*2] = NP.mean(self.visRcp[self.frequencyChannel, :, i*97+j])
                    if (phaseCorrect):
                        ewPh = self.ewAntPhaLcp[self.frequencyChannel, j]+self.ewLcpPhaseCorrection[self.frequencyChannel, j]
                        nPh = self.nAntPhaLcp[self.frequencyChannel, i]+self.nLcpPhaseCorrection[self.frequencyChannel, i]
                        self.uvLcp[O + (i+1)*2, O + (j-32)*2] *= NP.exp(1j * (-ewPh + nPh))
                        ewPh = self.ewAntPhaRcp[self.frequencyChannel, j]+self.ewRcpPhaseCorrection[self.frequencyChannel, j]
                        nPh = self.nAntPhaRcp[self.frequencyChannel, i]+self.nRcpPhaseCorrection[self.frequencyChannel, i]
                        self.uvRcp[O + (i+1)*2, O + (j-32)*2] *= NP.exp(1j * (-ewPh + nPh))
                    self.uvLcp[O - (i+1)*2, O - (j-32)*2] = NP.conj(self.uvLcp[O + (i+1)*2, O + (j-32)*2])
                    self.uvRcp[O - (i+1)*2, O - (j-32)*2] = NP.conj(self.uvRcp[O + (i+1)*2, O + (j-32)*2])
            for i in range(self.antNumberEW):
                if i<32:
                    self.uvLcp[O, O + (i-32)*2] = NP.mean(self.visLcp[self.frequencyChannel, :, self.antZeroRow[i]])
                    self.uvRcp[O, O + (i-32)*2] = NP.mean(self.visRcp[self.frequencyChannel, :, self.antZeroRow[i]])
                    if (phaseCorrect):
                        ewPh1 = self.ewAntPhaLcp[self.frequencyChannel, i]+self.ewLcpPhaseCorrection[self.frequencyChannel, i]
                        ewPh2 = self.ewAntPhaLcp[self.frequencyChannel, 32]+self.ewLcpPhaseCorrection[self.frequencyChannel, 32]
                        self.uvLcp[O, O + (i-32)*2] *= NP.exp(1j * (-ewPh1 + ewPh2))
                        ewPh1 = self.ewAntPhaRcp[self.frequencyChannel, i]+self.ewRcpPhaseCorrection[self.frequencyChannel, i]
                        ewPh2 = self.ewAntPhaRcp[self.frequencyChannel, 32]+self.ewRcpPhaseCorrection[self.frequencyChannel, 32]
                        self.uvRcp[O, O + (i-32)*2] *= NP.exp(1j * (-ewPh1 + ewPh2))
                    self.uvLcp[O, O + (32-i)*2] = NP.conj(self.uvLcp[O, O + (i-32)*2])
                    self.uvRcp[O, O + (32-i)*2] = NP.conj(self.uvRcp[O, O + (i-32)*2])
#                if i>32:
#                    self.uvLcp[O, O + (i-32)*2] = NP.conj(NP.mean(self.visLcp[self.frequencyChannel, :, self.antZeroRow[i]]))
#                    self.uvRcp[O, O + (i-32)*2] = NP.conj(NP.mean(self.visRcp[self.frequencyChannel, :, self.antZeroRow[i]]))
#                    if (phaseCorrect):
#                        ewPh1 = self.ewAntPhaLcp[self.frequencyChannel, i]+self.ewLcpPhaseCorrection[self.frequencyChannel, i]
#                        ewPh2 = self.ewAntPhaLcp[self.frequencyChannel, 32]+self.ewLcpPhaseCorrection[self.frequencyChannel, 32]
#                        self.uvLcp[O, O + (i-32)*2] *= NP.exp(1j * (ewPh1 - ewPh2))
#                        ewPh1 = self.ewAntPhaRcp[self.frequencyChannel, i]+self.ewRcpPhaseCorrection[self.frequencyChannel, i]
#                        ewPh2 = self.ewAntPhaRcp[self.frequencyChannel, 32]+self.ewRcpPhaseCorrection[self.frequencyChannel, 32]
#                        self.uvRcp[O, O + (i-32)*2] *= NP.exp(1j * (ewPh1 - ewPh2))
#                    self.uvLcp[O, O + (32-i)*2] = NP.conj(self.uvLcp[O, O + (i-32)*2])
#                    self.uvRcp[O, O + (32-i)*2] = NP.conj(self.uvRcp[O, O + (i-32)*2])
                
        else:
            for i in range(31):
                for j in range(self.antNumberEW):
                    self.uvLcp[O + (i+1)*2, O + (j-32)*2] = self.visLcp[self.frequencyChannel, scan, i*97+j]
                    self.uvRcp[O + (i+1)*2, O + (j-32)*2] = self.visRcp[self.frequencyChannel, scan, i*97+j]
                    if (phaseCorrect):
                        ewPh = self.ewAntPhaLcp[self.frequencyChannel, j]+self.ewLcpPhaseCorrection[self.frequencyChannel, j]
                        nPh = self.nAntPhaLcp[self.frequencyChannel, i]+self.nLcpPhaseCorrection[self.frequencyChannel, i]
                        self.uvLcp[O + (i+1)*2, O + (j-32)*2] *= NP.exp(1j * (-ewPh + nPh))
                        ewPh = self.ewAntPhaRcp[self.frequencyChannel, j]+self.ewRcpPhaseCorrection[self.frequencyChannel, j]
                        nPh = self.nAntPhaRcp[self.frequencyChannel, i]+self.nRcpPhaseCorrection[self.frequencyChannel, i]
                        self.uvRcp[O + (i+1)*2, O + (j-32)*2] *= NP.exp(1j * (-ewPh + nPh))
                    self.uvLcp[O - (i+1)*2, O - (j-32)*2] = NP.conj(self.uvLcp[O + (i+1)*2, O + (j-32)*2])
                    self.uvRcp[O - (i+1)*2, O - (j-32)*2] = NP.conj(self.uvRcp[O + (i+1)*2, O + (j-32)*2])
            for i in range(self.antNumberEW):
                if i<32:
                    self.uvLcp[O, O + (i-32)*2] = self.visLcp[self.frequencyChannel, scan, self.antZeroRow[i]]
                    self.uvRcp[O, O + (i-32)*2] = self.visRcp[self.frequencyChannel, scan, self.antZeroRow[i]]
                    if (phaseCorrect):
                        ewPh1 = self.ewAntPhaLcp[self.frequencyChannel, i]+self.ewLcpPhaseCorrection[self.frequencyChannel, i]
                        ewPh2 = self.ewAntPhaLcp[self.frequencyChannel, 32]+self.ewLcpPhaseCorrection[self.frequencyChannel, 32]
                        self.uvLcp[O, O + (i-32)*2] *= NP.exp(1j * (-ewPh1 + ewPh2))
                        ewPh1 = self.ewAntPhaRcp[self.frequencyChannel, i]+self.ewRcpPhaseCorrection[self.frequencyChannel, i]
                        ewPh2 = self.ewAntPhaRcp[self.frequencyChannel, 32]+self.ewRcpPhaseCorrection[self.frequencyChannel, 32]
                        self.uvRcp[O, O + (i-32)*2] *= NP.exp(1j * (-ewPh1 + ewPh2))
                    self.uvLcp[O, O + (32-i)*2] = NP.conj(self.uvLcp[O, O + (i-32)*2])
                    self.uvRcp[O, O + (32-i)*2] = NP.conj(self.uvRcp[O, O + (i-32)*2])
                # if i>32:
                #     self.uvLcp[O, O + (i-32)*2] = NP.conj(self.visLcp[self.frequencyChannel, scan, self.antZeroRow[i]])
                #     self.uvRcp[O, O + (i-32)*2] = NP.conj(self.visRcp[self.frequencyChannel, scan, self.antZeroRow[i]])
                #     if (phaseCorrect):
                #         ewPh1 = self.ewAntPhaLcp[self.frequencyChannel, i]+self.ewLcpPhaseCorrection[self.frequencyChannel, i]
                #         ewPh2 = self.ewAntPhaLcp[self.frequencyChannel, 32]+self.ewLcpPhaseCorrection[self.frequencyChannel, 32]
                #         self.uvLcp[O, O + (i-32)*2] *= NP.exp(1j * (ewPh1 - ewPh2))
                #         ewPh1 = self.ewAntPhaRcp[self.frequencyChannel, i]+self.ewRcpPhaseCorrection[self.frequencyChannel, i]
                #         ewPh2 = self.ewAntPhaRcp[self.frequencyChannel, 32]+self.ewRcpPhaseCorrection[self.frequencyChannel, 32]
                #         self.uvRcp[O, O + (i-32)*2] *= NP.exp(1j * (ewPh1 - ewPh2))
                #     # self.uvLcp[O, O + (32-i)*2] = NP.conj(self.uvLcp[O, O + (i-32)*2])
                #     # self.uvRcp[O, O + (32-i)*2] = NP.conj(self.uvRcp[O, O + (i-32)*2])

 
    def uv2lmImage(self):
        self.lcp = NP.fft.fft2(NP.roll(NP.roll(self.uvLcp,self.sizeOfUv//2+1,0),self.sizeOfUv//2+1,1));
        self.lcp = NP.roll(NP.roll(self.lcp,self.sizeOfUv//2-1,0),self.sizeOfUv//2-1,1);
        self.lcp = NP.flip(self.lcp, 1)
        self.rcp = NP.fft.fft2(NP.roll(NP.roll(self.uvRcp,self.sizeOfUv//2+1,0),self.sizeOfUv//2+1,1));
        self.rcp = NP.roll(NP.roll(self.rcp,self.sizeOfUv//2-1,0),self.sizeOfUv//2-1,1);
        self.rcp = NP.flip(self.rcp, 1)
        