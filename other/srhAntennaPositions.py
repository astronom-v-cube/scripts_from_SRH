#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 01:14:18 2022

@author: sergeyvlesovoi
"""

import numpy as NP
from BadaryRAO import BadaryRAO
from srh0612Utils import findFits
from scipy.stats import linregress

from srhFitsFile36 import SrhFitsFile as SrhFits0306
from srhFitsFile612 import SrhFitsFile as SrhFits0612
from srhFitsFile1224 import SrhFitsFile as SrhFits1224

class SrhAntennaPositions():
    def __init__(self, theDate):
        self.date = theDate
        self.deltaT = 300
        self.RAO = BadaryRAO(theDate, 4.9)
        self.earthOmega = NP.deg2rad(15/3600.)
        self.deltaHourAngle = 0.
        
    def setDate(self, theDate):
        self.date = theDate
        self.RAO.setDate(theDate)
        
    def getTimeMinusPi4(self):
        return self.RAO.culmination - (NP.pi/4 / self.earthOmega)

    def getTimePlusPi4(self):
        return self.RAO.culmination + (NP.pi/4 / self.earthOmega)
    
    def getTimeZero(self):
        return self.RAO.culmination
    
    def secsAsHhMmSs(self, t):
        hh = int(t / 3600.)
        t -= hh*3600.
        mm = int(t / 60.)
        t -= mm*60.
        ss = int(t)
        return '%02d:%02d:%02d' % (hh,mm,ss)
    
    def findSuitableFitses(self, srhArray, hAngle0):
#        srhFitsPath = '../SRH' + srhArray + '/' + self.date.replace('-','')
        srhFitsPath = '../SRH_DATA/SRH/SRH' + srhArray + '/' + self.date.replace('-','')
        fitsNames = findFits(srhFitsPath,'*.fit')
        fitsNames.sort()

        firstIndex = -1
        lastIndex = -1
        hPattern0 = ''
        hPattern1 = ''
        if (hAngle0 == 'minus_pi_4'):
            hPattern0 = self.secsAsHhMmSs(self.getTimeMinusPi4() - self.deltaT)[:-3].replace(':','')
            hPattern1 = self.secsAsHhMmSs(self.getTimeMinusPi4() + self.deltaT)[:-3].replace(':','')
        elif (hAngle0 == 'plus_pi_4'):
            hPattern0 = self.secsAsHhMmSs(self.getTimePlusPi4() - self.deltaT)[:-3].replace(':','')
            hPattern1 = self.secsAsHhMmSs(self.getTimePlusPi4() + self.deltaT)[:-3].replace(':','')
        elif (hAngle0 == 'zero'):
            hPattern0 = self.secsAsHhMmSs(self.getTimeZero() - self.deltaT)[:-3].replace(':','')
            hPattern1 = self.secsAsHhMmSs(self.getTimeZero() + self.deltaT)[:-3].replace(':','')

        pattern0Secs = int(hPattern0[0:2])*3600 + int(hPattern0[2:4])*60
        pattern1Secs = int(hPattern1[0:2])*3600 + int(hPattern1[2:4])*60

        for name in fitsNames:
#            fitsTime = name.split('T')[1][0:6]
            fitsTime = name.split('1224')[2].split('T')[1]
            fitsSecs = int(fitsTime[0:2])*3600 + int(fitsTime[2:4])*60 + int(fitsTime[4:6])
            if (NP.abs(fitsSecs - pattern0Secs) <= 60):
                firstIndex = fitsNames.index(name)
                break

        for name in fitsNames:
#            fitsTime = name.split('T')[1][0:6]
            fitsTime = name.split('1224')[2].split('T')[1]
            fitsSecs = int(fitsTime[0:2])*3600 + int(fitsTime[2:4])*60 + int(fitsTime[4:6])
            if (NP.abs(fitsSecs - pattern1Secs) <= 60):
                lastIndex = fitsNames.index(name)
                break

        print(firstIndex, lastIndex)
        if (firstIndex != -1 and lastIndex != -1):
            return fitsNames[firstIndex:lastIndex]
        else:
            return []

    def extractBaselines(self, srhArray, hAngle0):
        fitses = self.findSuitableFitses(srhArray, hAngle0)
        if (len(fitses) > 0):
            srhFits = SrhFits1224(fitses[0],1025)
            for fitsName in fitses[1:]:
                print(fitsName)
                srhFits.append(fitsName)
            scanNumber = srhFits.freqTime.shape[1]

            self.hourAngle = NP.zeros(scanNumber)
            for ss in range(scanNumber):
                self.hourAngle[ss] = srhFits.getHourAngle(ss)
            self.deltaHourAngle = self.hourAngle[scanNumber//2 + 1] - self.hourAngle[scanNumber//2]

            lcpVisBaseline = NP.zeros((48,scanNumber))
            rcpVisBaseline = NP.zeros((48,scanNumber))
            
            freqRange = 10
#            freqRange = srhFits.freqListLength
            for vv in range(48):
                for ff in range(freqRange):
                    lcpVisBaseline[vv] += NP.unwrap(NP.angle(srhFits.visLcp[ff,:,11775 + vv]))/(2*NP.pi*srhFits.freqList[ff]*1e3)*3e8
                    rcpVisBaseline[vv] += NP.unwrap(NP.angle(srhFits.visRcp[ff,:,11775 + vv]))/(2*NP.pi*srhFits.freqList[ff]*1e3)*3e8
                lcpVisBaseline[vv] /= freqRange
                rcpVisBaseline[vv] /= freqRange
            srhFits.close()
            return [self.hourAngle, 0.5*(lcpVisBaseline + rcpVisBaseline)]
        else:
            return []

    def calcSlopes(self, baselines):
        baselineSlope = []
        for bb in range(baselines[1].shape[0]):
            bSlope, interc, A,B,C = linregress(baselines[0],baselines[1][bb])
            baselineSlope.append(bSlope)
        return NP.array(baselineSlope)
        