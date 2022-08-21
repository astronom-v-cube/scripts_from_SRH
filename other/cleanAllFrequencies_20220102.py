#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 10:45:32 2021

@author: sergeyvlesovoi
"""
import pylab as PL
from astropy.io import fits
import os
from skimage import measure
import casatasks

def calcEllipseArea(cell, level):
    contours = measure.find_contours(cell, level)
    contourLength = []
    for n, contour in enumerate(contours):
        contourLength.append(len(contour))
    contourLength = NP.array(contourLength)
    
    contourMaxInd = NP.argmax(contourLength)
    srcEll = measure.EllipseModel()
    srcEll.estimate(contours[contourMaxInd])
    exc, eyc, ea, eb, theta = srcEll.params
    return NP.sqrt(NP.pi*ea*eb)

def distanceMeanFrom(array, rc, distance):
    return 0.25* (\
        array[rc[0] - distance,rc[1]] + \
        array[rc[0] + distance,rc[1]] + \
        array[rc[0],           rc[1] - distance] + \
        array[rc[0],           rc[1] + distance])

curDate = '20220102'
freq0 = int(phaseEdit.srhFits.freqList[0] * 1e-3)
dFreq = int((phaseEdit.srhFits.freqList[1] - phaseEdit.srhFits.freqList[0]) * 1e-3)
freqs = phaseEdit.srhFits.freqList*1e3

msPath = '/home/sergey_lesovoi/SRH_06_12_data/20220102/'
for ff in range(16):
    print(ff)
    phaseEdit.onFrequencyChannelChanged(ff)
    phaseEdit.onCenterDisk()
    phaseEdit.saveAsMs2(msPath + '%d.ms'%(int(phaseEdit.srhFits.freqList[ff]*1e-3)))


imagePath = '/home/sergey_lesovoi/SRH_06_12_data/20220102/images/cln_'
for curFreq in range(16):
    visName = msPath + '%d.ms'%(freq0 + curFreq*400)
    imageName = imagePath+'%d'%(freq0 + curFreq*400)
    print(visName)
    print(imageName)
    casatasks.tclean(vis=visName, imagename=imageName ,cell=2.45,imsize=1024,threshold='30000Jy',stokes='RRLL',niter=10000)


clnFitsPath = '/home/sergey_lesovoi/SRH_06_12_data/20220102/images/srh_%s_cln_'%(curDate)
for curFreq in range(16):
    casatasks.exportfits(imagename=imagePath + '%d.image'%(freq0 + curFreq*400), fitsimage=clnFitsPath+'%d.fit'%(freq0 + curFreq*400),history=False)    
    
clnFits = []
lcpImages = []
rcpImages = []
for ff in range(16):
    fitName = ('/home/sergey_lesovoi/SRH_06_12_data/20220102/images/srh_%s_cln_%d.fit'%(curDate,int(phaseEdit.srhFits.freqList[ff]*1e-3)))
    clnFits.append(fits.open(fitName))
    lcpImages.append(clnFits[-1][0].data[0,0])
    rcpImages.append(clnFits[-1][0].data[1,0])

lcpImages_20211230_0 = NP.array(lcpImages)
rcpImages_20211230_0 = NP.array(rcpImages)

Y0 = 340
X0 = 590
dXY = 128
srcLcp = lcpImages_20211230_0[:,Y0:Y0+dXY,X0:X0+dXY]
srcRcp = rcpImages_20211230_0[:,Y0:Y0+dXY,X0:X0+dXY]

rightSrcYX = [55,72]
leftSrcYX = [55,40]
srcDXY = 50
srcDX = 25
srcDY = 50

fig = PL.figure()
fig.suptitle('LCP')
pl = fig.subplots(nrows=4, ncols=4)
for rr in range(4):
    for cc in range(4):
        curImage = srcLcp[rr*4 + cc].copy()
        curImage[rightSrcYX[0]:rightSrcYX[0] + srcDY,rightSrcYX[1]:rightSrcYX[1] + srcDX] += 1e5
        curImage[leftSrcYX[0]:leftSrcYX[0] + srcDY,leftSrcYX[1]:leftSrcYX[1] + srcDX] += 2e5
        pl[rr,cc].imshow(curImage,origin='lower',vmin=2e4,vmax=5e5)
        pl[rr,cc].axis('off')
fig.tight_layout()

fig = PL.figure()
fig.suptitle('RCP')
pl = fig.subplots(nrows=4, ncols=4)
for rr in range(4):
    for cc in range(4):
        curImage = srcRcp[rr*4 + cc].copy()
        curImage[rightSrcYX[0]:rightSrcYX[0] + srcDY,rightSrcYX[1]:rightSrcYX[1] + srcDX] += 1e5
        curImage[leftSrcYX[0]:leftSrcYX[0] + srcDY,leftSrcYX[1]:leftSrcYX[1] + srcDX] += 2e5
        pl[rr,cc].imshow(curImage,origin='lower',vmin=2e4,vmax=5e5)
        pl[rr,cc].axis('off')
fig.tight_layout()

rightSrcLcp = srcLcp[:,rightSrcYX[0]:rightSrcYX[0] + srcDY,rightSrcYX[1]:rightSrcYX[1] + srcDX].max(axis=(1,2))
rightSrcRcp = srcRcp[:,rightSrcYX[0]:rightSrcYX[0] + srcDY,rightSrcYX[1]:rightSrcYX[1] + srcDX].max(axis=(1,2))
leftSrcLcp =  srcLcp[:,leftSrcYX[0]:leftSrcYX[0]   + srcDY,leftSrcYX[1]:leftSrcYX[1] + srcDX].max(axis=(1,2))
leftSrcRcp =  srcRcp[:,leftSrcYX[0]:leftSrcYX[0]   + srcDY,leftSrcYX[1]:leftSrcYX[1] + srcDX].max(axis=(1,2))

rightSrc1Lcp = NP.zeros(16)
rightSrc2Lcp = NP.zeros(16)
rightSrc3Lcp = NP.zeros(16)
rightSrc4Lcp = NP.zeros(16)

rightSrc1Rcp = NP.zeros(16)
rightSrc2Rcp = NP.zeros(16)
rightSrc3Rcp = NP.zeros(16)
rightSrc4Rcp = NP.zeros(16)

leftSrc1Lcp = NP.zeros(16)
leftSrc2Lcp = NP.zeros(16)
leftSrc3Lcp = NP.zeros(16)
leftSrc4Lcp = NP.zeros(16)

leftSrc1Rcp = NP.zeros(16)
leftSrc2Rcp = NP.zeros(16)
leftSrc3Rcp = NP.zeros(16)
leftSrc4Rcp = NP.zeros(16)

rightSrcLcpSize = NP.zeros(16)
rightSrcRcpSize = NP.zeros(16)
leftSrcLcpSize = NP.zeros(16)
leftSrcRcpSize = NP.zeros(16)
sizeLevel = 0.7

for ff in range(16):
    curSrc = srcRcp[ff,rightSrcYX[0]:rightSrcYX[0] + srcDY,rightSrcYX[1]:rightSrcYX[1] + srcDX]
    rightSrcRcpSize[ff] = calcEllipseArea(curSrc,sizeLevel*curSrc.max())
    
    maxInd = NP.unravel_index(NP.argmax(curSrc),curSrc.shape)
    # rightSrc1Rcp[ff] = distanceMeanFrom(curSrc,maxInd,1)
    # rightSrc2Rcp[ff] = distanceMeanFrom(curSrc,maxInd,2)
    #rightSrc3Rcp[ff] = distanceMeanFrom(curSrc,maxInd,3)
    #rightSrc4Rcp[ff] = distanceMeanFrom(curSrc,maxInd,4)

    curSrc = srcLcp[ff,rightSrcYX[0]:rightSrcYX[0] + srcDY,rightSrcYX[1]:rightSrcYX[1] + srcDX]
    rightSrcLcpSize[ff] = calcEllipseArea(curSrc,sizeLevel*curSrc.max())

    maxInd = NP.unravel_index(NP.argmax(curSrc),curSrc.shape)
    # rightSrc1Lcp[ff] = distanceMeanFrom(curSrc,maxInd,1)
    # rightSrc2Lcp[ff] = distanceMeanFrom(curSrc,maxInd,2)
    #rightSrc3Lcp[ff] = distanceMeanFrom(curSrc,maxInd,3)
    #rightSrc4Lcp[ff] = distanceMeanFrom(curSrc,maxInd,4)

    curSrc = srcRcp[ff,leftSrcYX[0]:leftSrcYX[0] + srcDY,leftSrcYX[1]:leftSrcYX[1] + srcDX]
    leftSrcRcpSize[ff] = calcEllipseArea(curSrc,sizeLevel*curSrc.max())

    maxInd = NP.unravel_index(NP.argmax(curSrc),curSrc.shape)
    # leftSrc1Rcp[ff] = distanceMeanFrom(curSrc,maxInd,1)
    # leftSrc2Rcp[ff] = distanceMeanFrom(curSrc,maxInd,2)
    # leftSrc3Rcp[ff] = distanceMeanFrom(curSrc,maxInd,3)
    # leftSrc4Rcp[ff] = distanceMeanFrom(curSrc,maxInd,4)

    curSrc = srcLcp[ff,leftSrcYX[0]:leftSrcYX[0] + srcDY,leftSrcYX[1]:leftSrcYX[1] + srcDX]
    leftSrcLcpSize[ff] = calcEllipseArea(curSrc,sizeLevel*curSrc.max())

    maxInd = NP.unravel_index(NP.argmax(curSrc),curSrc.shape)
    # leftSrc1Lcp[ff] = distanceMeanFrom(curSrc,maxInd,1)
    # leftSrc2Lcp[ff] = distanceMeanFrom(curSrc,maxInd,2)
    # leftSrc3Lcp[ff] = distanceMeanFrom(curSrc,maxInd,3)
    # leftSrc4Lcp[ff] = distanceMeanFrom(curSrc,maxInd,4)

fig = PL.figure()
pl = fig.subplots()
fig.suptitle('right source')
pl.plot(freqs,rightSrcLcp,'-',color='blue')
pl.plot(freqs,rightSrc1Lcp,'--',color='blue',lw=.5)
pl.plot(freqs,rightSrc2Lcp,'--',color='blue',lw=.5)
pl.plot(freqs,rightSrc3Lcp,'--',color='blue',lw=.5)
pl.plot(freqs,rightSrc4Lcp,'--',color='blue',lw=.5)
pl.plot(freqs,rightSrcRcp,'-',color='red')
pl.plot(freqs,rightSrc1Rcp,'--',color='red',lw=.5)
pl.plot(freqs,rightSrc2Rcp,'--',color='red',lw=.5)
pl.plot(freqs,rightSrc3Rcp,'--',color='red',lw=.5)
pl.plot(freqs,rightSrc4Rcp,'--',color='red',lw=.5)
pl.set_ylabel('K')
pl.set_ylim(0,6e5)
pl.set_xlabel('GHz')
pl.grid()

tpl = pl.twinx()
tpl.set_ylabel('arcsec')
tpl.set_ylim(5,20)
tpl.plot(freqs,rightSrcLcpSize,'o-',color='skyblue')
tpl.plot(freqs,rightSrcRcpSize,'o-',color='tomato')

fig = PL.figure()
pl = fig.subplots()
fig.suptitle('left source')
pl.plot(freqs,leftSrcLcp,'-',color='blue')
pl.plot(freqs,leftSrc1Lcp,'--',color='blue',lw=.5)
pl.plot(freqs,leftSrc2Lcp,'--',color='blue',lw=.5)
pl.plot(freqs,leftSrc3Lcp,'--',color='blue',lw=.5)
pl.plot(freqs,leftSrc4Lcp,'--',color='blue',lw=.5)
pl.plot(freqs,leftSrcRcp,'-',color='red')
pl.plot(freqs,leftSrc1Rcp,'--',color='red',lw=.5)
pl.plot(freqs,leftSrc2Rcp,'--',color='red',lw=.5)
pl.plot(freqs,leftSrc3Rcp,'--',color='red',lw=.5)
pl.plot(freqs,leftSrc4Rcp,'--',color='red',lw=.5)
pl.set_ylabel('K')
pl.set_ylim(0,6e5)
pl.set_xlabel('GHz')
pl.grid()

tpl = pl.twinx()
tpl.set_ylabel('arcsec')
tpl.set_ylim(5,20)
tpl.plot(freqs,leftSrcLcpSize,'o-',color='skyblue')
tpl.plot(freqs,leftSrcRcpSize,'o-',color='tomato')

