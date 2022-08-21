#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 15:23:27 2021

@author: sergeyvlesovoi
"""
from scipy import signal
import scipy.constants
import srh36Utils
import ZirinTb
import numpy as NP

def scaleTo1(array):
    arrMin = array.min()
    arrMax = array.max()
    return (array - arrMin)/(arrMax - arrMin)

def hms_format(t, pos):
  hh = int(t / 3600.)
  t -= hh*3600.
  mm = int(t / 60.)
  t -= mm*60.
  ss = int(t)
  return '%02d:%02d:%02d' % (hh,mm,ss);

try:
    phaseEdit
except NameError:
    print('phaseEdit is not defined')
    
Zi = ZirinTb.ZirinTb()

#startFluxes = [75.5, 83.9, 88.9, 100, 100, 100]
startFluxes = Zi.getSfuAtFrequncy(phaseEdit.srhFits.freqList*1e-6)/2.3

imageSize = phaseEdit.uvSize // 2
N = 32
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
rMaxList = []
lMaxList = []
rCell = NP.zeros((phaseEdit.srhFits.freqListLength, phaseEdit.srhFits.dataLength,M,M))
lCell = NP.zeros((phaseEdit.srhFits.freqListLength, phaseEdit.srhFits.dataLength,M,M))

srcYX  = NP.array([[308, 188],  [314, 194], [319, 216], [333, 216],  [338, 119], [300, 144]])
qSunYX = NP.array([[420, 280], [420, 280], [212, 283], [212, 314], [169, 261], [169, 261]])
skyYX  = NP.array([[471, 161],  [471, 161],  [471, 161], [345, 138],  [244, 240],   [244, 240]])

# for i in range(srcYX.shape[0]):
#     srcYX[i,0] = 512 - srcYX[i,0]
#     qSunYX[i,0] = 512 - qSunYX[i,0]
#     skyYX[i,0] = 512 - skyYX[i,0]
    
t0 = [30,30,30,30,30,30]
F0 = [1.2, 1.2, 1.2, 1.2, 1.2, 1.2]
dYdX = NP.array([[15, 12], [15, 12], [15, 12], [15, 12], [15, 12], [15, 12]])
Su = [1,1,1,1,2,1]
Sk = [1,1,1,1,0,1]

lCellSun = NP.zeros((phaseEdit.srhFits.freqListLength))
rCellSun = NP.zeros((phaseEdit.srhFits.freqListLength))
lCellSky = NP.zeros((phaseEdit.srhFits.freqListLength))
rCellSky = NP.zeros((phaseEdit.srhFits.freqListLength))
lCellSrc = NP.zeros((phaseEdit.srhFits.freqListLength,phaseEdit.srhFits.dataLength))
rCellSrc = NP.zeros((phaseEdit.srhFits.freqListLength,phaseEdit.srhFits.dataLength))
iCellSrc = NP.zeros((phaseEdit.srhFits.freqListLength,phaseEdit.srhFits.dataLength))
iCellMax = NP.zeros((phaseEdit.srhFits.freqListLength,phaseEdit.srhFits.dataLength))
vCellMax = NP.zeros((phaseEdit.srhFits.freqListLength,phaseEdit.srhFits.dataLength))

for j in NP.arange(0,2,1):
#for j in NP.arange(0,phaseEdit.srhFits.freqListLength,1):
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
    rCellScaled = NP.zeros((phaseEdit.srhFits.dataLength,N,N))
    lCellScaled = NP.zeros((phaseEdit.srhFits.dataLength,N,N))
   
    
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
            rCellScaled[:,k,l] = signal.detrend(rCellPlot[:,k,l])
            rCellScaled[:,k,l] -= rCellPlot[:,k,l].mean()
            lCellScaled[:,k,l] = signal.detrend(lCellPlot[:,k,l])
            lCellScaled[:,k,l] -= lCellPlot[:,k,l].mean()
    
    rCorrCoefs = NP.zeros((N,N))
    lCorrCoefs = NP.zeros((N,N))
    for k in range(N):
        for l in range(N):
            rCorrCoefs[k,l] = NP.mean(rCorrPlot*rCellScaled[:,k,l])
            lCorrCoefs[k,l] = NP.mean(lCorrPlot*lCellScaled[:,k,l])
    
    lMax = NP.unravel_index(NP.argmax(lCorrCoefs),lCorrCoefs.shape)
    rMax = NP.unravel_index(NP.argmax(rCorrCoefs),lCorrCoefs.shape)
    lMaxList.append(lMax)
    rMaxList.append(rMax)
    for i in range(phaseEdit.srhFits.dataLength):
        if j == 4:
            lCell[j,i,:,:] = lImages[i,lMax[0]*M + 4:(lMax[0]+1)*M + 4, lMax[1]*M + 4:(lMax[1]+1)*M + 4]
            rCell[j,i,:,:] = rImages[i,rMax[0]*M + 4:(rMax[0]+1)*M + 4, rMax[1]*M + 4:(rMax[1]+1)*M + 4]
        else:            
            lCell[j,i,:,:] = lImages[i,lMax[0]*M:(lMax[0]+1)*M, lMax[1]*M:(lMax[1]+1)*M]
            rCell[j,i,:,:] = rImages[i,rMax[0]*M:(rMax[0]+1)*M, rMax[1]*M:(rMax[1]+1)*M]
    
    rCorrFlux.append(rCorrPlot)
    lCorrFlux.append(lCorrPlot)
    rAmpFlux.append(rFluxPlot)
    lAmpFlux.append(lFluxPlot)
    rCellFlux.append(rCellPlot[:,lMax[0],lMax[1]])
    lCellFlux.append(lCellPlot[:,lMax[0],lMax[1]])

    lCellSun[j] = lImages[:,qSunYX[j][0]:qSunYX[j][0] + dYdX[j][0], qSunYX[j][1]:qSunYX[j][1] + dYdX[j][1]].mean()
    rCellSun[j] = rImages[:,qSunYX[j][0]:qSunYX[j][0] + dYdX[j][0], qSunYX[j][1]:qSunYX[j][1] + dYdX[j][1]].mean()
    lCellSky[j] = lImages[:,skyYX[j][0] :skyYX[j][0]  + dYdX[j][0], skyYX[j][1]: skyYX[j][1]  + dYdX[j][1]].mean()
    rCellSky[j] = rImages[:,skyYX[j][0] :skyYX[j][0]  + dYdX[j][0], skyYX[j][1]: skyYX[j][1]  + dYdX[j][1]].mean()

    qSunTb = 0.5*Zi.getTbAtFrequncy(phaseEdit.srhFits.freqList[j]*1e-6)*1e3
    lImages[:] -= lCellSky[j]
    lImages[:] /= (lCellSun[j] - lCellSky[j])
    lImages[:] *= qSunTb
    rImages[:] -= rCellSky[j]
    rImages[:] /= (rCellSun[j] - rCellSky[j])
    rImages[:] *= qSunTb
    
    lCellSrc[j] = (lImages[:,srcYX[j][0]:srcYX[j][0] + dYdX[j][0], srcYX[j][1]:srcYX[j][1] + dYdX[j][1]]).sum(axis=(1,2))
    rCellSrc[j] = (rImages[:,srcYX[j][0]:srcYX[j][0] + dYdX[j][0], srcYX[j][1]:srcYX[j][1] + dYdX[j][1]]).sum(axis=(1,2))

    iCellSrc[j] = lCellSrc[j] + rCellSrc[j]
    iCellSrc[j] = 2*iCellSrc[j]*scipy.constants.k/(scipy.constants.c/(phaseEdit.srhFits.freqList[j]*1e3))**2 * NP.deg2rad(4.911/3600)**2 / 1e-22

    for ii in range(phaseEdit.srhFits.dataLength):
        lFoundCell = lImages[ii,srcYX[j][0]:srcYX[j][0] + dYdX[j][0], srcYX[j][1]:srcYX[j][1] + dYdX[j][1]]
        curMax = NP.array(NP.unravel_index(lFoundCell.argmax(),lFoundCell.shape))
        curMax[0] += srcYX[j][0]
        curMax[1] += srcYX[j][1]
        lCurVal = lImages[ii,curMax[0]-3:curMax[0]+4,curMax[1]-3:curMax[1]+4].sum(axis=(0,1))
        rFoundCell = rImages[ii,srcYX[j][0]:srcYX[j][0] + dYdX[j][0], srcYX[j][1]:srcYX[j][1] + dYdX[j][1]]
        curMax = NP.array(NP.unravel_index(rFoundCell.argmax(),rFoundCell.shape))
        curMax[0] += srcYX[j][0]
        curMax[1] += srcYX[j][1]
        rCurVal = rImages[ii,curMax[0]-3:curMax[0]+4,curMax[1]-3:curMax[1]+4].sum(axis=(0,1))
        iCellMax[j,ii] = 2*(lCurVal + rCurVal)*scipy.constants.k/(scipy.constants.c/(phaseEdit.srhFits.freqList[j]*1e3))**2 * NP.deg2rad(4.911/3600)**2 / 1e-22
        vCellMax[j,ii] = 2*(lCurVal - rCurVal)*scipy.constants.k/(scipy.constants.c/(phaseEdit.srhFits.freqList[j]*1e3))**2 * NP.deg2rad(4.911/3600)**2 / 1e-22

rCorrFlux = NP.array(rCorrFlux)
lCorrFlux = NP.array(lCorrFlux)
rAmpFlux = NP.array(rAmpFlux)
lAmpFlux = NP.array(lAmpFlux)
rCellFlux = NP.array(rCellFlux)
lCellFlux = NP.array(lCellFlux)

fig = PL.figure(figsize=(12,10))
fig.suptitle('SRH ' + srhCpFitFile[0].header['DATE-OBS'] + ' flux')
pl = fig.subplots(nrows=1,ncols=1)
pl.xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
pl.xaxis.set_major_locator(PL.MultipleLocator(60))
for f in range(phaseEdit.srhFits.freqListLength):
    pl.plot(phaseEdit.srhFits.freqTime[f],iCellMax[f],label='%d MHz'%int(phaseEdit.srhFits.freqList[f]*1e-3))
pl.set_xlabel('Time UT')
pl.set_ylabel('s.f.u.')
pl.legend()
pl.grid()
pl.set_xlim(3*3600 + 18*60,3*3600 + 21*60)
pl.set_ylim(0,1.5)

fig = PL.figure(figsize=(12,10))
fig.suptitle('SRH ' + srhCpFitFile[0].header['DATE-OBS'] + ' flux - flux(04:32:52)')
pl = fig.subplots(nrows=1,ncols=1)
pl.xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
pl.xaxis.set_major_locator(PL.MultipleLocator(30))
pl.yaxis.set_major_locator(PL.MultipleLocator(.025))
for f in range(phaseEdit.srhFits.freqListLength):
    pl.plot(phaseEdit.srhFits.freqTime[f],iCellMax[f] - iCellMax[f,30],label='%d MHz'%int(phaseEdit.srhFits.freqList[f]*1e-3))
pl.set_xlabel('Time UT')
pl.set_ylabel('s.f.u.')
pl.legend()
pl.grid()
pl.set_xlim(3*3600 + 18*60,3*3600 + 21*60)
pl.set_ylim(-0.01,.39)

#lMaxVal = NP.zeros((phaseEdit.srhFits.freqListLength,phaseEdit.srhFits.dataLength))
#rMaxVal = NP.zeros((phaseEdit.srhFits.freqListLength,phaseEdit.srhFits.dataLength))
#for j in NP.arange(0,phaseEdit.srhFits.freqListLength,1):
#    for i in range(phaseEdit.srhFits.dataLength):
#        curMax = NP.unravel_index(lCell[j,i].argmax(),(M,M))
#        for kk in NP.arange(-1,2,1):
#            for ll in NP.arange(-1,2,1):
#                lMaxVal[j,i] += lCell[j,i,curMax[0]+kk,curMax[1]+ll]
#        curMax = NP.unravel_index(rCell[j,i].argmax(),(M,M))
#        for kk in NP.arange(-1,2,1):
#            for ll in NP.arange(-1,2,1):
#                rMaxVal[j,i] += rCell[j,i,curMax[0]+kk,curMax[1]+ll]
#

# PL.figure()
# PL.imshow(lImages[0],origin='lower')
# PL.plot([lMax[1]*M, lMax[1]*M, (lMax[1]+1)*M, (lMax[1]+1)*M, lMax[1]*M],
#         [lMax[0]*M,(lMax[0]+1)*M,(lMax[0]+1)*M, lMax[0]*M, lMax[0]*M],color='white')

#flux = 2*gIImage[row][256-200:256+200,256-200:256+200].sum()*scipy.constants.k/(scipy.constants.c/gFrequency[row]*1e-6)**2 * NP.deg2rad(gArcsecPerPixel/3600)**2 / 1e-22

#----------------------------------------------------
#yS = 277
#xS = 208
#
#yQ = 244
#xQ = 185
#
#y0 = 183
#x0 = 53
#
#t0 = 30
#F0 = 1.2
#lCellSun5600 = lImages[:,yQ:yQ+15,xQ:xQ+12].mean()
#rCellSun5600 = rImages[:,yQ:yQ+15,xQ:xQ+12].mean()
#lCellSky5600 = lImages[:,y0:y0+15,x0:x0+12].mean()
#rCellSky5600 = rImages[:,y0:y0+15,x0:x0+12].mean()
#
#lCellCurve5600 = (lImages[:,yS:yS+15,xS:xS+12] - lCellSky5600 - lCellSun5600).sum(axis=(1,2))
#rCellCurve5600 = (rImages[:,yS:yS+15,xS:xS+12] - rCellSky5600 - rCellSun5600).sum(axis=(1,2))
#
#iCellCurve5600 = lCellCurve5600 + rCellCurve5600
#iCellCurve5600  = iCellCurve5600 / iCellCurve5600[t0] * F0 #s.f.u.
#----------------------------------------------------
#yS = 396
#xS = 136

#yQ = 323
#xQ = 210

#y0 = 0
#x0 = 0
#lCellSun4700 = lImages[:,yQ:yQ+15,xQ:xQ+12].mean()
#rCellSun4700 = rImages[:,yQ:yQ+15,xQ:xQ+12].mean()
#lCellSky4700 = lImages[:,y0:y0+15,x0:x0+12].mean()
#rCellSky4700 = rImages[:,y0:y0+15,x0:x0+12].mean()
#
#lCellCurve4700 = (lImages[:,yS:yS+15,xS:xS+12] - 2*lCellSun4700).sum(axis=(1,2))
#rCellCurve4700 = (rImages[:,yS:yS+15,xS:xS+12] - 2*rCellSun4700).sum(axis=(1,2))
#
#iCellCurve4700 = lCellCurve4700 + rCellCurve4700
#iCellCurve4700  = iCellCurve4700 / iCellCurve4700[0] * 1.2 #s.f.u.
#
#----------------------------------------------------
#yS = 222
#xS = 94

#yQ = 360
#xQ = 190

#y0 = 400
#x0 = 97

#lCellSun3900 = lImages[:,yQ:yQ+15,xQ:xQ+12].mean()
#rCellSun3900 = rImages[:,yQ:yQ+15,xQ:xQ+12].mean()
#lCellSky3900 = lImages[:,y0:y0+15,x0:x0+12].mean()
#rCellSky3900 = rImages[:,y0:y0+15,x0:x0+12].mean()
#
#lCellCurve3900 = (lImages[:,yS:yS+15,xS:xS+12] - lCellSky3900 - lCellSun3900).sum(axis=(1,2))
#rCellCurve3900 = (rImages[:,yS:yS+15,xS:xS+12] - rCellSky3900 - rCellSun3900).sum(axis=(1,2))
#
#iCellCurve3900 = lCellCurve3900 + rCellCurve3900
#iCellCurve3900  = iCellCurve3900 / iCellCurve3900[0] * 1.2 #s.f.u.

#----------------------------------------------------

#yS = 256
#xS = 212

#yQ = 180
#xQ = 340

#y0 = 395
#x0 = 120

#lCellSun3400 = lImages[:,yQ:yQ+15,xQ:xQ+12].mean()
#rCellSun3400 = rImages[:,yQ:yQ+15,xQ:xQ+12].mean()
#lCellSky3400 = lImages[:,y0:y0+15,x0:x0+12].mean()
#rCellSky3400 = rImages[:,y0:y0+15,x0:x0+12].mean()
#
#lCellCurve3400 = (lImages[:,yS:yS+15,xS:xS+12] - lCellSky3400 - lCellSun3400).sum(axis=(1,2))
#rCellCurve3400 = (rImages[:,yS:yS+15,xS:xS+12] - rCellSky3400 - rCellSun3400).sum(axis=(1,2))
#
#iCellCurve3400 = lCellCurve3400 + rCellCurve3400
#iCellCurve3400  = iCellCurve3400 / iCellCurve3400[0] * 1.2 #s.f.u.
#
#----------------------------------------------------
#
#yS = 243
#xS = 137

#yQ = 125
#xQ = 102

#y0 = 387
#x0 = 34

#lCellSun3100 = lImages[:,yQ:yQ+15,xQ:xQ+12].mean()
#rCellSun3100 = rImages[:,yQ:yQ+15,xQ:xQ+12].mean()
#lCellSky3100 = lImages[:,y0:y0+15,x0:x0+12].mean()
#rCellSky3100 = rImages[:,y0:y0+15,x0:x0+12].mean()
#
#lCellCurve3100 = (lImages[:,yS:yS+15,xS:xS+12] - lCellSky3100 - lCellSun3100).sum(axis=(1,2))
#rCellCurve3100 = (rImages[:,yS:yS+15,xS:xS+12] - rCellSky3100 - rCellSun3100).sum(axis=(1,2))
#
#iCellCurve3100 = lCellCurve3100 + rCellCurve3100
#iCellCurve3100  = iCellCurve3100 / iCellCurve3100[0] * 1.2 #s.f.u.
#
#----------------------------------------------------

#yS = 220
#xS = 39

#yQ = 150
#xQ = 324

#y0 = 374
#x0 = 52

#dy = 15
#dx = 12
#lCellSun2800 = lImages[:,yQ:yQ+dy,xQ:xQ+dx].mean()
#rCellSun2800 = rImages[:,yQ:yQ+dy,xQ:xQ+dx].mean()
#lCellSky2800 = lImages[:,y0:y0+dy,x0:x0+dx].mean()
#rCellSky2800 = rImages[:,y0:y0+dy,x0:x0+dx].mean()
#
#lCellCurve2800 = (lImages[:,yS:yS+dy,xS:xS+dx] - lCellSky2800 - lCellSun2800).sum(axis=(1,2))
#rCellCurve2800 = (rImages[:,yS:yS+dy,xS:xS+dx] - rCellSky2800 - rCellSun2800).sum(axis=(1,2))
#
#iCellCurve2800 = lCellCurve2800 + rCellCurve2800
#iCellCurve2800 = iCellCurve2800 / iCellCurve2800[0] * 1.2 #s.f.u.
