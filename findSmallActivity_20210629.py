#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 15:23:27 2021

@author: sergeyvlesovoi
"""
from scipy import signal
import srh36Utils
import ZirinTb

def scaleTo1(array):
    arrMin = array.min()
    arrMax = array.max()
    return (array - arrMin)/(arrMax - arrMin)

try:
    phaseEdit
except NameError:
    print('phaseEdit is not defined')
    
Zi = ZirinTb.ZirinTb()

#startFluxes = [75.5, 83.9, 88.9, 100, 100, 100]
startFluxes = Zi.getSfuAtFrequency(phaseEdit.srhFits.freqList*1e-6)/2.3

imageSize = phaseEdit.uvSize // 2
N = 64
M = imageSize // N

rcpMeanAmps = NP.mean(NP.abs(phaseEdit.srhFits.ampRcp[0,:,:]),axis=0)
validAmpsIndMask = []
invalidAmpsIndMask = []
for i in range(phaseEdit.srhFits.antennaNames.shape[0]):
    validAmpsIndMask.append(rcpMeanAmps[i] < 1e10 and rcpMeanAmps[i] > 1e3)
    invalidAmpsIndMask.append(not(rcpMeanAmps[i] < 1e10 and rcpMeanAmps[i] > 1e3))
    
validAmpsInd = NP.array(NP.where(validAmpsIndMask))
invalidAmpsInd = NP.array(NP.where(invalidAmpsIndMask))
invalidAntNumbers = NP.array(phaseEdit.srhFits.antennaNumbers[invalidAmpsInd],dtype='int')
n2nDict = srh36Utils.number2name()
#for num in range(invalidAntNumbers.size):
#    print(n2nDict[invalidAntNumbers[0][num]])

rCorrFlux = []
lCorrFlux = []
rAmpFlux = []
lAmpFlux = []
rCellFlux = []
lCellFlux = []

for j in range(phaseEdit.srhFits.freqListLength):
    phaseEdit.onFrequencyChannelChanged(j)
    print('frequency %d'%phaseEdit.currentFrequencyChannel)
    rImages = []
    lImages = []
    rCorrPlot = NP.zeros((phaseEdit.srhFits.dataLength))
    lCorrPlot = NP.zeros((phaseEdit.srhFits.dataLength))
    rFluxPlot = NP.zeros((phaseEdit.srhFits.dataLength))
    lFluxPlot = NP.zeros((phaseEdit.srhFits.dataLength))
    rCellPlot = NP.zeros((phaseEdit.srhFits.dataLength,N,N))
    lCellPlot = NP.zeros((phaseEdit.srhFits.dataLength,N,N))
   
    
    # visNumber = phaseEdit.srhFits.visLcp.shape[2]
    # for i in range(visNumber):  
    #     if phaseEdit.srhFits.antennaA[i] in invalidAntNumbers or phaseEdit.srhFits.antennaB[i] in invalidAntNumbers:
    #         phaseEdit.srhFits.visLcp[phaseEdit.currentFrequencyChannel,:,i] *= 0 + 0j
    #         phaseEdit.srhFits.visRcp[phaseEdit.currentFrequencyChannel,:,i] *= 0 + 0j
    
    for i in range(phaseEdit.srhFits.dataLength):
        phaseEdit.currentScan = i; 
        phaseEdit.buildImage(); 
        rImages.append(phaseEdit.rcpData)
        lImages.append(phaseEdit.lcpData)
        lCorrPlot[i] = NP.mean(NP.abs(phaseEdit.srhFits.visLcp[phaseEdit.currentFrequencyChannel, i, :]))
        rCorrPlot[i] = NP.mean(NP.abs(phaseEdit.srhFits.visRcp[phaseEdit.currentFrequencyChannel, i, :]))
        lFluxPlot[i] = NP.mean(NP.abs(phaseEdit.srhFits.ampLcp[phaseEdit.currentFrequencyChannel, i, validAmpsInd]))
        rFluxPlot[i] = NP.mean(NP.abs(phaseEdit.srhFits.ampRcp[phaseEdit.currentFrequencyChannel, i, validAmpsInd]))
        
    rImages = NP.array(rImages)
    lImages = NP.array(lImages)
    
    rCorrPlot -= rCorrPlot.mean()
    lCorrPlot -= lCorrPlot.mean()
    rCorrPlot /= rCorrPlot.max()
    lCorrPlot /= lCorrPlot.max()
    
    for i in range(phaseEdit.srhFits.dataLength):
        for k in range(N):
            for l in range(N):
                rCellPlot[i,k,l] = rImages[i,k*M:(k+1)*M, l*M:(l+1)*M].sum(axis=(0,1))
                lCellPlot[i,k,l] = lImages[i,k*M:(k+1)*M, l*M:(l+1)*M].sum(axis=(0,1))
    
    for k in range(N):
        for l in range(N):
            rCellPlot[:,k,l] = signal.detrend(rCellPlot[:,k,l])
            rCellPlot[:,k,l] -= rCellPlot[:,k,l].mean()
            lCellPlot[:,k,l] = signal.detrend(lCellPlot[:,k,l])
            lCellPlot[:,k,l] -= lCellPlot[:,k,l].mean()
    
    rCorrCoefs = NP.zeros((N,N))
    lCorrCoefs = NP.zeros((N,N))
    for k in range(N):
        for l in range(N):
            rCorrCoefs[k,l] = NP.mean(rCorrPlot*rCellPlot[:,k,l])
            lCorrCoefs[k,l] = NP.mean(lCorrPlot*lCellPlot[:,k,l])
    
    lMax = NP.unravel_index(NP.argmax(lCorrCoefs),lCorrCoefs.shape)
    rMax = NP.unravel_index(NP.argmax(rCorrCoefs),lCorrCoefs.shape)
    
    rCorrFlux.append(rCorrPlot)
    lCorrFlux.append(lCorrPlot)
    rAmpFlux.append(rFluxPlot)
    lAmpFlux.append(lFluxPlot)
    rCellFlux.append(rCellPlot[:,lMax[0],lMax[1]])
    lCellFlux.append(lCellPlot[:,lMax[0],lMax[1]])

rCorrFlux = NP.array(rCorrFlux)
lCorrFlux = NP.array(lCorrFlux)
rAmpFlux = NP.array(rAmpFlux)
lAmpFlux = NP.array(lAmpFlux)
rCellFlux = NP.array(rCellFlux)
lCellFlux = NP.array(lCellFlux)
