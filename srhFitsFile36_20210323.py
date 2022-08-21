# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 00:18:47 2016

@author: Sergey
"""

from astropy.io import fits
import numpy as NP
from scipy import ndimage
from astropy import coordinates
from astropy import constants
from BadaryRAO import BadaryRAO
from scipy.optimize import least_squares

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
        self.averageCalib = False
        self.useNonlinearApproach = False
        self.ewLcpPhaseCorrection = NP.zeros(32)
        self.ewRcpPhaseCorrection = NP.zeros(32)
        self.sLcpPhaseCorrection = NP.zeros(16)
        self.sRcpPhaseCorrection = NP.zeros(16)
        self.badAntsLcp = NP.zeros((32, 48))
        self.badAntsRcp = NP.zeros((32, 48))
        self.fringeStopping = False
        self.phaseCalibFull = True
        self.updateAmpCalibration = True
        self.sizeOfUv = sizeOfUv
        self.baselines = 1
        self.open(name)
        self.ewAntPhaLcp = NP.zeros(65)
        self.nAntPhaLcp = NP.zeros(31)
        self.ewAntPhaRcp = NP.zeros(65)
        self.nAntPhaRcp = NP.zeros(31)
                                    
    def open(self,name):
        try:
            self.hduList = fits.open(name)
            self.isOpen = True
            self.dateObs = self.hduList[0].header['DATE-OBS'] + 'T' + self.hduList[0].header['TIME-OBS']
            self.antennaNumbers = self.hduList[2].data['ant_index']
            self.antennaNumbers = NP.reshape(self.antennaNumbers,self.antennaNumbers.size)
            self.antennaA = self.hduList[4].data['ant_B']
            self.antennaA = NP.reshape(self.antennaA,self.antennaA.size)
            self.antennaB = self.hduList[4].data['ant_A']
            self.antennaB = NP.reshape(self.antennaB,self.antennaB.size)
            self.uvLcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex)
            self.uvRcp = NP.zeros((self.sizeOfUv,self.sizeOfUv),dtype=complex)
            self.freqList = self.hduList[1].data['frequency'];
            self.freqListLength = self.freqList.size;
            self.dataLength = self.hduList[1].data['time'].size // self.freqListLength;
#            self.freqTime = NP.reshape(self.hduList[1].data['TIME'],(self.freqListLength,self.dataLength));
            self.freqTime = self.hduList[1].data['time']
            self.visListLength = self.hduList[1].data['vis_lcp'].size // self.freqListLength // self.dataLength;
            self.visLcp = NP.reshape(self.hduList[1].data['vis_lcp'],(self.freqListLength,self.dataLength,self.visListLength));
            self.visRcp = NP.reshape(self.hduList[1].data['vis_rcp'],(self.freqListLength,self.dataLength,self.visListLength));
            self.visLcp /= float(self.hduList[0].header['VIS_MAX'])
            self.visRcp /= float(self.hduList[0].header['VIS_MAX'])
            self.ampLcp = NP.reshape(self.hduList[1].data['amp_lcp'],(self.freqListLength,self.dataLength,self.antennaNumbers.size));
            self.ampRcp = NP.reshape(self.hduList[1].data['amp_rcp'],(self.freqListLength,self.dataLength,self.antennaNumbers.size));
#            self.visLcp.real = NP.sin(NP.pi/2*self.visLcp.real);
#            self.visRcp.real = NP.sin(NP.pi/2*self.visRcp.real);
#            self.visLcp.imag = NP.sin(NP.pi/2*self.visLcp.imag);
#            self.visRcp.imag = NP.sin(NP.pi/2*self.visRcp.imag);

            self.antZeroRow = self.hduList[3].data['ant_zero_row'][:96]
            self.RAO = BadaryRAO(self.dateObs.split('T')[0])

        except FileNotFoundError:
            print('File %s  not found'%name);
    
    def append(self,name):
        try:
            hduList = fits.open(name);
#            freqTime = NP.reshape(self.hduList[1].data['TIME'],(self.freqListLength,self.dataLength));
            freqTime = hduList[1].data['time']
            dataLength = hduList[1].data['time'].size // self.freqListLength;
            visLcp = NP.reshape(hduList[1].data['vis_lcp'],(self.freqListLength,dataLength,self.visListLength));
            visRcp = NP.reshape(hduList[1].data['vis_rcp'],(self.freqListLength,dataLength,self.visListLength));
            ampLcp = NP.reshape(hduList[1].data['amp_lcp'],(self.freqListLength,dataLength,self.antennaNumbers.size));
            ampRcp = NP.reshape(hduList[1].data['amp_rcp'],(self.freqListLength,dataLength,self.antennaNumbers.size));
#            visLcp.real = NP.sin(NP.pi/2*visLcp.real);
#            visRcp.real = NP.sin(NP.pi/2*visRcp.real);
#            visLcp.imag = NP.sin(NP.pi/2*visLcp.imag);
#            visRcp.imag = NP.sin(NP.pi/2*visRcp.imag);

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
    
#    def updateAntennaAmplitude(self, average = 0):

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

    def updateAntennaPhase(self, baselinesNumber = 1, nonlinear = True):
        if nonlinear:
            self.calculatePhaseLcp_nonlinear(baselinesNumber = baselinesNumber)
#            self.calculatePhaseRcp_nonlinear(baselinesNumber = baselinesNumber)
        else:
            self.calculatePhase_linear(baselinesNumber = baselinesNumber)
        
    def calculatePhase_linear(self, baselinesNumber = 1):
        antNumberN = 31
        antNumberEW = 65
        redIndexesN = []
        for baseline in range(1, baselinesNumber+1):
            redIndexesN.append(NP.where((self.antennaA==98-1+baseline) & (self.antennaB==33))[0][0])
            for i in range(antNumberN - baseline):
                redIndexesN.append(NP.where((self.antennaA==98+i) & (self.antennaB==98+i+baseline))[0][0])
    
        redIndexesEW = []
        for baseline in range(1, baselinesNumber+1):
            for i in range(antNumberEW - baseline):
                redIndexesEW.append(NP.where((self.antennaA==1+i) & (self.antennaB==1+i+baseline))[0][0])
                
        phaMatrix = self.phaMatrixGenPairsEWN(baselinesNumber, antNumberEW, antNumberN)
        redundantVisLcp = self.visLcp[self.frequencyChannel, self.calibIndex, NP.append(redIndexesEW, redIndexesN)]
        sunVisPhases = NP.zeros(2)
        if self.freqList[self.frequencyChannel] > 4e6:
            sunVisPhases = NP.array((NP.pi, NP.pi))
        phasesLcp = NP.concatenate((NP.angle(redundantVisLcp), sunVisPhases))
        antPhaLcp, c, d, e = NP.linalg.lstsq(phaMatrix, phasesLcp, rcond=None)
        self.ewAntPhaLcp = antPhaLcp[baselinesNumber*2:baselinesNumber*2+antNumberEW]
        self.nAntPhaLcp = antPhaLcp[baselinesNumber*2+antNumberEW:]
        
        redundantVisRcp = self.visRcp[self.frequencyChannel, self.calibIndex, NP.append(redIndexesEW, redIndexesN)]
        phasesRcp = NP.concatenate((NP.angle(redundantVisRcp), NP.array((0,0))))
        antPhaRcp, c, d, e = NP.linalg.lstsq(phaMatrix, phasesRcp, rcond=None)
        self.ewAntPhaRcp = antPhaRcp[baselinesNumber*2:baselinesNumber*2+antNumberEW]
        self.nAntPhaRcp = antPhaRcp[baselinesNumber*2+antNumberEW:]
        
    def calculatePhaseLcp_nonlinear(self, baselinesNumber = 2):
        antNumberN = 31
        antNumberEW = 65
        redIndexesN = []
        for baseline in range(1, baselinesNumber+1):
            redIndexesN.append(NP.where((self.antennaA==98-1+baseline) & (self.antennaB==33))[0][0])
            for i in range(antNumberN - baseline):
                redIndexesN.append(NP.where((self.antennaA==98+i) & (self.antennaB==98+i+baseline))[0][0])
    
        redIndexesEW = []
        for baseline in range(1, baselinesNumber+1):
            for i in range(antNumberEW - baseline):
                redIndexesEW.append(NP.where((self.antennaA==1+i) & (self.antennaB==1+i+baseline))[0][0])
             
        redundantVisN = self.visLcp[self.frequencyChannel, self.calibIndex, redIndexesN]
#        redundantVisN[16] *= 0.001 
#        redundantVisN[21] *= 0.001 
#        redundantVisN[28] *= 0.001 
        x_ini = NP.concatenate((NP.ones(baselinesNumber+antNumberN-1), NP.zeros(baselinesNumber+antNumberN-1)))
        ls_res = least_squares(self.northGainsFunc_constrained, x_ini, args = (redundantVisN, antNumberN, baselinesNumber), max_nfev = 50)
        self.n_gains_lcp = self.real_to_complex(ls_res['x'])[baselinesNumber-1:]
        self.nAntPhaLcp = NP.angle(self.n_gains_lcp)
        
        redundantVisEW = self.visLcp[self.frequencyChannel, self.calibIndex, redIndexesEW]
#        redundantVisEW[22] *= 0.001
#        redundantVisEW[24] *= 0.001
#        redundantVisEW[44] *= 0.001
        x_ini = NP.concatenate((NP.ones(baselinesNumber+antNumberEW-2), NP.zeros(baselinesNumber+antNumberEW-2)))
        ls_res = least_squares(self.eastWestGainsFunc_constrained, x_ini, args = (redundantVisEW, antNumberEW, baselinesNumber), max_nfev = 100)
        gains = self.real_to_complex(ls_res['x'])[baselinesNumber-1:]
        self.ew_gains_lcp = NP.insert(gains, antNumberEW//2, (1+0j))
        self.ewAntPhaLcp = NP.angle(self.ew_gains_lcp)
        
    def calculatePhaseRcp_nonlinear(self, baselinesNumber = 2):
        antNumberN = 31
        antNumberEW = 65
        redIndexesN = []
        for baseline in range(1, baselinesNumber+1):
            redIndexesN.append(NP.where((self.antennaA==98-1+baseline) & (self.antennaB==33))[0][0])
            for i in range(antNumberN - baseline):
                redIndexesN.append(NP.where((self.antennaA==98+i) & (self.antennaB==98+i+baseline))[0][0])
    
        redIndexesEW = []
        for baseline in range(1, baselinesNumber+1):
            for i in range(antNumberEW - baseline):
                redIndexesEW.append(NP.where((self.antennaA==1+i) & (self.antennaB==1+i+baseline))[0][0])
             
        redundantVisN = self.visRcp[self.frequencyChannel, self.calibIndex, redIndexesN]
        x_ini = NP.concatenate((NP.ones(baselinesNumber+antNumberN-1), NP.zeros(baselinesNumber+antNumberN-1)))
        ls_res = least_squares(self.northGainsFunc_constrained, x_ini, args = (redundantVisN, antNumberN, baselinesNumber), max_nfev = 150)
        self.n_gains_rcp = self.real_to_complex(ls_res['x'])[baselinesNumber-1:]
        self.nAntPhaRcp = NP.angle(self.n_gains_rcp)
        
        redundantVisEW = self.visRcp[self.frequencyChannel, self.calibIndex, redIndexesEW]
        x_ini = NP.concatenate((NP.ones(baselinesNumber+antNumberEW-2), NP.zeros(baselinesNumber+antNumberEW-2)))
        ls_res = least_squares(self.eastWestGainsFunc_constrained, x_ini, args = (redundantVisEW, antNumberEW, baselinesNumber), max_nfev = 300)
        gains = self.real_to_complex(ls_res['x'])[baselinesNumber-1:]
        self.ew_gains_rcp = NP.insert(gains, antNumberEW//2, (1+0j))
        self.ewAntPhaRcp = NP.angle(self.ew_gains_rcp)
                
    def eastWestGainsFunc_constrained(self, x, obsVis, antNumber, baselineNumber):
        res = NP.zeros_like(obsVis, dtype = complex)
        x_complex = self.real_to_complex(x)
        solVis = NP.append((1+0j), x_complex[:baselineNumber-1])
        gains = NP.insert(x_complex[baselineNumber-1:], antNumber//2, (1+0j))
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
    
    def real_to_complex(self, z):
        return z[:len(z)//2] + 1j * z[len(z)//2:]
    
    def complex_to_real(self, z):
        return NP.concatenate((NP.real(z), NP.imag(z)))
    
    def setCalibIndex(self, calibIndex):
        self.calibIndex = calibIndex;
        self.updateAntennaPhase(baselinesNumber = self.baselines, nonlinear = self.useNonlinearApproach)

    def setFrequencyChannel(self, channel):
        self.frequencyChannel = channel
        self.updateAntennaPhase(baselinesNumber = self.baselines, nonlinear = self.useNonlinearApproach)
        
    def vis2uv(self, scan, phaseCorrect = True):
        self.uvLcp[:,:] = complex(0,0)
        self.uvRcp[:,:] = complex(0,0)
        O = self.sizeOfUv//2
        for i in range(31):
            for j in range(65):
                self.uvLcp[O + (i+1)*2, O + (j-32)*2] = self.visLcp[self.frequencyChannel, scan, i*97+j]
                self.uvRcp[O + (i+1)*2, O + (j-32)*2] = self.visRcp[self.frequencyChannel, scan, i*97+j]
                if (phaseCorrect):
                    self.uvLcp[O + (i+1)*2, O + (j-32)*2] *= NP.exp(1j * (-self.ewAntPhaLcp[j] + self.nAntPhaLcp[i]))
                    self.uvRcp[O + (i+1)*2, O + (j-32)*2] *= NP.exp(1j * (-self.ewAntPhaRcp[j] + self.nAntPhaRcp[i]))
                self.uvLcp[O - (i+1)*2, O - (j-32)*2] = NP.conj(self.uvLcp[O + (i+1)*2, O + (j-32)*2])
                self.uvRcp[O - (i+1)*2, O - (j-32)*2] = NP.conj(self.uvRcp[O + (i+1)*2, O + (j-32)*2])
        for i in range(32):
            self.uvLcp[O, O + (i-32)*2] = self.visLcp[self.frequencyChannel, scan, self.antZeroRow[i]]
            self.uvRcp[O, O + (i-32)*2] = self.visRcp[self.frequencyChannel, scan, self.antZeroRow[i]]
            if (phaseCorrect):
                self.uvLcp[O, O + (i-32)*2] *= NP.exp(1j * (-self.ewAntPhaLcp[i] + self.ewAntPhaLcp[32]))
                self.uvRcp[O, O + (i-32)*2] *= NP.exp(1j * (-self.ewAntPhaRcp[i] + self.ewAntPhaRcp[32]))
            self.uvLcp[O, O + (32-i)*2] = NP.conj(self.uvLcp[O, O + (i-32)*2])
            self.uvRcp[O, O + (32-i)*2] = NP.conj(self.uvRcp[O, O + (i-32)*2])
 
    def uv2lmImage(self):
        self.lcp = NP.fft.fft2(NP.roll(NP.roll(self.uvLcp,self.sizeOfUv//2+1,0),self.sizeOfUv//2+1,1));
        self.lcp = NP.roll(NP.roll(self.lcp,self.sizeOfUv//2-1,0),self.sizeOfUv//2-1,1);
        self.lcp = NP.flip(self.lcp, 1)
        self.rcp = NP.fft.fft2(NP.roll(NP.roll(self.uvRcp,self.sizeOfUv//2+1,0),self.sizeOfUv//2+1,1));
        self.rcp = NP.roll(NP.roll(self.rcp,self.sizeOfUv//2-1,0),self.sizeOfUv//2-1,1);
        self.rcp = NP.flip(self.rcp, 1)
        
    def ij2hd(self, ij):
        return(((ij[0] - self.centerH)*self.deltaH,(ij[1] - self.centerD)*self.deltaD));

    def hd2pq(self, hd):
        return((hd[0],hd[1]));

    def pq2kl(self, pq):
        return((pq[0]/self.deltaP + self.centerP,pq[1]/self.deltaQ + self.centerQ));

    def ij2kl(self, ij):
        return(self.pq2kl(self.hd2pq(self.ij2hd(ij))));
        
    def lmImage2hdImage(self):
        return(ndimage.geometric_transform(self.lcp.real,self.ij2kl));
    
    def getDateObs(self):
        return self.dateObs if self.isOpen else ''

    def getDataLength(self):
        return self.dataLength if self.isOpen else 0

    def getTimesObs(self):
        return self.hduList[1].data['TIME'] if self.isOpen else 0

    def getAntennaA(self):
        return self.antennaA if self.isOpen else 0
        
    def getAntennaB(self):
        return self.antennaB if self.isOpen else 0

    def getVisLcp(self):
        return self.visLcp if self.isOpen else 0
        
    def getVisRcp(self):
        return self.visRcp if self.isOpen else 0
    
    def getFrequencyList(self):
        return self.freqList
    
    def phaseClosure(self, antA, antB, antC):
        pC = NP.zeros(self.dataLength)
        for i in range(self.dataLength):
            pC[i] = NP.angle(self.visLcp[5,i,antA]) + NP.angle(self.visLcp[5,i,antB]) - NP.angle(self.visLcp[5,i,antC])
        return pC

    def phaseAnt(self, antA):
        pC = NP.zeros(self.dataLength)
        for i in range(self.dataLength):
            pC[i] = NP.angle(self.visLcp[5,i,antA]) 
        return pC

    def magnitAnt(self, antA):
        pC = NP.zeros(self.dataLength)
        for i in range(self.dataLength):
            pC[i] = NP.abs(self.visLcp[5,i,antA]) 
        return pC

    def phaseClosureSouthVector(self, visEw, scan, freq):
        pC = NP.zeros(16)
        for i in range(15):
            pC[i + 1] = NP.angle(self.visLcp[freq,scan,512 + 15 + visEw]) + NP.angle(self.visLcp[freq,scan,32*i + visEw + 1]) - NP.angle(self.visLcp[freq,scan,32*i + visEw])
        return pC

    def phaseClosureEastWestVector(self, visS, scan, freq):
        pC = NP.zeros(32)
        pC[0] = 1.
        for i in range(31):
            pC[i + 1] = - NP.angle(self.visLcp[freq,scan,16*visS + i]) + NP.angle(self.visLcp[freq,scan,512 + visS]) - NP.angle(self.visLcp[freq,scan,16*visS + i + 1]) 
        return pC
    
    def solveSouthPhaseClosure(self, pcVector):
        phases, c, d, e = NP.linalg.lstsq(self.sPhaseClosureMatrix,pcVector);
        return phases

    def solveEastWestPhaseClosure(self, pcVector):
        phases, c, d, e = NP.linalg.lstsq(self.ewPhaseClosureMatrix,pcVector);
        return phases
    