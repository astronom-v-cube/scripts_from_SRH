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
startFluxes = Zi.getSfuAtFrequncy(phaseEdit.srhFits.freqList*1e-6)/2.3

imageSize = phaseEdit.uvSize // 2
N = 64
M = imageSize // N

fluxR = []
fluxL = []

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
for num in range(invalidAntNumbers.size):
    print(n2nDict[invalidAntNumbers[0][num]])

for j in range(2):#phaseEdit.srhFits.freqListLength):
    phaseEdit.onFrequencyChannelChanged(j)
    rImages = []
    lImages = []
    rCorrPlot = NP.zeros((phaseEdit.srhFits.dataLength))
    lCorrPlot = NP.zeros((phaseEdit.srhFits.dataLength))
    rFluxPlot = NP.zeros((phaseEdit.srhFits.dataLength))
    lFluxPlot = NP.zeros((phaseEdit.srhFits.dataLength))
    rFullFluxPlot = NP.zeros((phaseEdit.srhFits.dataLength))
    lFullFluxPlot = NP.zeros((phaseEdit.srhFits.dataLength))
    
    visNumber = phaseEdit.srhFits.visLcp.shape[2]
    for i in range(visNumber):  
        if phaseEdit.srhFits.antennaA[i] in invalidAntNumbers or phaseEdit.srhFits.antennaB[i] in invalidAntNumbers:
            phaseEdit.srhFits.visLcp[phaseEdit.currentFrequencyChannel,:,i] *= 0 + 0j
            phaseEdit.srhFits.visRcp[phaseEdit.currentFrequencyChannel,:,i] *= 0 + 0j
    
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
    
    rImagesSmall = NP.zeros((phaseEdit.srhFits.dataLength,N,N))
    lImagesSmall = NP.zeros((phaseEdit.srhFits.dataLength,N,N))
    
    for i in range(phaseEdit.srhFits.dataLength):
        for k in range(N):
            for l in range(N):
                rImagesSmall[i,k,l] = rImages[i,k*M:(k+1)*M, l*M:(l+1)*M].sum(axis=(0,1))
                lImagesSmall[i,k,l] = lImages[i,k*M:(k+1)*M, l*M:(l+1)*M].sum(axis=(0,1))
                rFullFluxPlot[i] = rImages[i].sum()
                lFullFluxPlot[i] = lImages[i].sum()
    
    for k in range(N):
        for l in range(N):
            rImagesSmall[:,k,l] = signal.detrend(rImagesSmall[:,k,l])
            rImagesSmall[:,k,l] -= rImagesSmall[:,k,l].mean()
            lImagesSmall[:,k,l] = signal.detrend(lImagesSmall[:,k,l])
            lImagesSmall[:,k,l] -= lImagesSmall[:,k,l].mean()
    
    rCorrCoefs = NP.zeros((N,N))
    lCorrCoefs = NP.zeros((N,N))
    for k in range(N):
        for l in range(N):
            rCorrCoefs[k,l] = NP.mean(rCorrPlot*rImagesSmall[:,k,l])
            lCorrCoefs[k,l] = NP.mean(lCorrPlot*lImagesSmall[:,k,l])
    
    lMax = NP.unravel_index(NP.argmax(lCorrCoefs),lCorrCoefs.shape)
    rMax = NP.unravel_index(NP.argmax(rCorrCoefs),lCorrCoefs.shape)
    
    fluxR.append(rImagesSmall[:,lMax[0],lMax[1]])
    fluxL.append(lImagesSmall[:,lMax[0],lMax[1]])

# PL.figure()
# PL.plot(scaleTo1(rFluxPlot))
# PL.plot(scaleTo1(rCorrPlot))
# PL.plot(scaleTo1(rImagesSmall[:,lMax[0],lMax[1]]))

# PL.figure()
# PL.imshow(lImages[0],origin='lower')
# PL.plot([lMax[1]*M, lMax[1]*M, (lMax[1]+1)*M, (lMax[1]+1)*M, lMax[1]*M],
#         [lMax[0]*M,(lMax[0]+1)*M,(lMax[0]+1)*M, lMax[0]*M, lMax[0]*M],color='white')

# meanFlux = rFluxPlot[0:75].mean()
# rCalFlux = rFluxPlot - 0.5*meanFlux

#freq x
# PL.figure()
# PL.xlabel('UTC')
# PL.ylabel('s.f.u.')
# PL.title('%d MHz'%(phaseEdit.srhFits.freqList[phaseEdit.currentFrequencyChannel]/1000))
# PL.plot(rCalFlux/rCalFlux[0:75].mean()*startFluxes[phaseEdit.currentFrequencyChannel])
# PL.plot((lImagesSmall[:,lMax[0],lMax[1]]+rImagesSmall[:,lMax[0],lMax[1]])*.005+startFluxes[phaseEdit.currentFrequencyChannel])
# PL.grid()
