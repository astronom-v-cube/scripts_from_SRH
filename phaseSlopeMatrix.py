#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 01:11:10 2020

@author: svlesovoi
"""
import numpy as NP

class AntennaDelayMatrix():
    def __init__(self, antennaNumber, frequencies):
        frequencyNumber = frequencies.shape[0]
        rows = (antennaNumber - 1)*frequencyNumber #equations
        cols = antennaNumber + frequencyNumber
        self.visAntMatrix = NP.zeros((rows, cols))
        for r in range(rows):
            visPhaseCol = r % frequencyNumber
            freqIndex = visPhaseCol 
            antDelayCol = r // frequencyNumber + frequencyNumber
            self.visAntMatrix[r, visPhaseCol] = 1
            self.visAntMatrix[r, antDelayCol:antDelayCol + 2] = frequencies[freqIndex], -frequencies[freqIndex]
            self.pinvVisAntMatrix = NP.linalg.pinv(self.visAntMatrix)


