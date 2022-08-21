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

#def distanceMeanFrom(array, rc, distance):
#    return 0.25* (\
#        array[rc[0] - distance,rc[1]] + \
#        array[rc[0] + distance,rc[1]] + \
#        array[rc[0],           rc[1] - distance] + \
#        array[rc[0],           rc[1] + distance])
#
def distanceMeanFrom(array, rc, distance):
    if distance == 1:
        return (array[rc[0] - 1,rc[1] + 1])
    elif distance == 2:
        return (array[rc[0] - 1,rc[1] - 1])
    elif distance == 3:
        return (array[rc[0] + 1,rc[1] + 1])
    elif distance == 4:
        return (array[rc[0] + 1,rc[1] - 1])

curDate = '20211112'
freq0 = int(phaseEdit.srhFits.freqList[0] * 1e-3)
dFreq = int((phaseEdit.srhFits.freqList[1] - phaseEdit.srhFits.freqList[0]) * 1e-3)
freqs = phaseEdit.srhFits.freqList*1e3

msPath = '/home/sergeyvlesovoi/SRH_0612/20211112_0/'
for ff in range(16):
    print(ff)
    phaseEdit.onFrequencyChannelChanged(ff)
    phaseEdit.onCenterDisk()
    phaseEdit.saveAsMs2(msPath + 'srh%s_%d.ms'%(curDate,int(phaseEdit.srhFits.freqList[ff]*1e-3)))

with open(msPath + 'cleanAllFreqs.py', 'w') as f:
    f.write('freq0 = %d\n'%freq0)
    f.write('curDate = %s\n'%curDate)
    f.write('msPath = \'/home/sergeyvlesovoi/SRH_0612/20211112_0/srh%s_\'\n'%(curDate))
    f.write('imagePath = \'/home/sergeyvlesovoi/SRH_0612/20211112_0/images/cln_\'\n')
    f.write('for curFreq in range(16):\n')
    f.write('\ttclean(vis=msPath+\'%d.ms\'%(freq0 + curFreq*400),imagename=imagePath+\'%d\'%(freq0 + curFreq*400),cell=2.45,imsize=1024,threshold=\'30000Jy\',stokes=\'RRLL\',niter=10000)')
os.system('casa -c ' + msPath + 'cleanAllFreqs.py')

with open(msPath + 'saveFitses.py', 'w') as f:
    f.write('freq0 = %d\n'%freq0)
    f.write('imagePath = \'/home/sergeyvlesovoi/SRH_0612/20211112_0/images/cln_\'\n')
    f.write('clnFitsPath = \'/home/sergeyvlesovoi/SRH_0612/20211112_0/images/srh_%s_cln_\'\n'%(curDate))
    f.write('for curFreq in range(16):\n')
    f.write('\texportfits(imagename=imagePath+\'%d.image\'%(freq0 + curFreq*400),fitsimage=clnFitsPath+\'%d.fit\'%(freq0 + curFreq*400),history=False)')
os.system('casa -c ' + msPath + 'saveFitses.py')

clnFits = []
lcpImages = []
rcpImages = []
for ff in range(16):
    fitName = ('/home/sergeyvlesovoi/SRH_0612/20211112_0/images/srh_%s_cln_%d.fit'%(curDate,int(phaseEdit.srhFits.freqList[ff]*1e-3)))
    clnFits.append(fits.open(fitName))
    lcpImages.append(clnFits[-1][0].data[0,0])
    rcpImages.append(clnFits[-1][0].data[1,0])

lcpImages_20211112_0 = NP.array(lcpImages)
rcpImages_20211112_0 = NP.array(rcpImages)

Y0S = 286
X0S = 575
dXY = 128
srcLcpS = lcpImages_20211112_0[:,Y0S:Y0S+dXY,X0S:X0S+dXY]
srcRcpS = rcpImages_20211112_0[:,Y0S:Y0S+dXY,X0S:X0S+dXY]
Y0N = 688
X0N = 760
srcLcpN = lcpImages_20211112_0[:,Y0N:Y0N+dXY,X0N:X0N+dXY]
srcRcpN = rcpImages_20211112_0[:,Y0N:Y0N+dXY,X0N:X0N+dXY]

srcYXS = [55,72]
srcYXN = [48,40]
srcDXY = 50
srcDX = 20
srcDY = 30

fig = PL.figure()
fig.suptitle('LCP S')
pl = fig.subplots(nrows=4, ncols=4)
for rr in range(4):
    for cc in range(4):
        curImage = srcLcpS[rr*4 + cc].copy()
        curImage[srcYXS[0]:srcYXS[0] + srcDY,srcYXS[1]:srcYXS[1] + srcDX] += 3e4
        pl[rr,cc].imshow(curImage,origin='lower',vmin=2e4,vmax=1e5)
        pl[rr,cc].axis('off')
fig.tight_layout()

fig = PL.figure()
fig.suptitle('RCP S')
pl = fig.subplots(nrows=4, ncols=4)
for rr in range(4):
    for cc in range(4):
        curImage = srcRcpS[rr*4 + cc].copy()
        curImage[srcYXS[0]:srcYXS[0] + srcDY,srcYXS[1]:srcYXS[1] + srcDX] += 3e4
        pl[rr,cc].imshow(curImage,origin='lower',vmin=2e4,vmax=1e5)
        pl[rr,cc].axis('off')
fig.tight_layout()

fig = PL.figure()
fig.suptitle('LCP N')
pl = fig.subplots(nrows=4, ncols=4)
for rr in range(4):
    for cc in range(4):
        curImage = srcLcpN[rr*4 + cc].copy()
        curImage[srcYXN[0]:srcYXN[0] + srcDY,srcYXN[1]:srcYXN[1] + srcDX] += 3e4
        pl[rr,cc].imshow(curImage,origin='lower',vmin=2e4,vmax=1e5)
        pl[rr,cc].axis('off')
fig.tight_layout()

fig = PL.figure()
fig.suptitle('RCP N')
pl = fig.subplots(nrows=4, ncols=4)
for rr in range(4):
    for cc in range(4):
        curImage = srcRcpN[rr*4 + cc].copy()
        curImage[srcYXN[0]:srcYXN[0] + srcDY,srcYXN[1]:srcYXN[1] + srcDX] += 3e4
        pl[rr,cc].imshow(curImage,origin='lower',vmin=2e4,vmax=1e5)
        pl[rr,cc].axis('off')
fig.tight_layout()

src0LcpS = srcLcpS[:,srcYXS[0]:srcYXS[0] + srcDY,srcYXS[1]:srcYXS[1] + srcDX].max(axis=(1,2))
src0RcpS = srcRcpS[:,srcYXS[0]:srcYXS[0] + srcDY,srcYXS[1]:srcYXS[1] + srcDX].max(axis=(1,2))
src0LcpN = srcLcpN[:,srcYXN[0]:srcYXN[0]   + srcDY,srcYXN[1]:srcYXN[1] + srcDX].max(axis=(1,2))
src0RcpN = srcRcpN[:,srcYXN[0]:srcYXN[0]   + srcDY,srcYXN[1]:srcYXN[1] + srcDX].max(axis=(1,2))

src1LcpS = NP.zeros(16)
src2LcpS = NP.zeros(16)
src3LcpS = NP.zeros(16)
src4LcpS = NP.zeros(16)

src1RcpS = NP.zeros(16)
src2RcpS = NP.zeros(16)
src3RcpS = NP.zeros(16)
src4RcpS = NP.zeros(16)

src1LcpN = NP.zeros(16)
src2LcpN = NP.zeros(16)
src3LcpN = NP.zeros(16)
src4LcpN = NP.zeros(16)

src1RcpN = NP.zeros(16)
src2RcpN = NP.zeros(16)
src3RcpN = NP.zeros(16)
src4RcpN = NP.zeros(16)

srcLcpSSize = NP.zeros(16)
srcRcpSSize = NP.zeros(16)
srcLcpNSize = NP.zeros(16)
srcRcpNSize = NP.zeros(16)
sizeLevel = 0.7

for ff in range(16):
    curSrc = srcRcpS[ff,srcYXS[0]:srcYXS[0] + srcDY,srcYXS[1]:srcYXS[1] + srcDX]
    srcRcpSSize[ff] = calcEllipseArea(curSrc,sizeLevel*curSrc.max())
    
    maxInd = NP.unravel_index(NP.argmax(curSrc),curSrc.shape)
    src1RcpS[ff] = distanceMeanFrom(curSrc,maxInd,1)
    src2RcpS[ff] = distanceMeanFrom(curSrc,maxInd,2)
    src3RcpS[ff] = distanceMeanFrom(curSrc,maxInd,3)
    src4RcpS[ff] = distanceMeanFrom(curSrc,maxInd,4)

    curSrc = srcLcpS[ff,srcYXS[0]:srcYXS[0] + srcDY,srcYXS[1]:srcYXS[1] + srcDX]
    srcLcpSSize[ff] = calcEllipseArea(curSrc,sizeLevel*curSrc.max())

    maxInd = NP.unravel_index(NP.argmax(curSrc),curSrc.shape)
    src1LcpS[ff] = distanceMeanFrom(curSrc,maxInd,1)
    src2LcpS[ff] = distanceMeanFrom(curSrc,maxInd,2)
    src3LcpS[ff] = distanceMeanFrom(curSrc,maxInd,3)
    src4LcpS[ff] = distanceMeanFrom(curSrc,maxInd,4)

    curSrc = srcRcpN[ff,srcYXN[0]:srcYXN[0] + srcDY,srcYXN[1]:srcYXN[1] + srcDX]
    srcRcpNSize[ff] = calcEllipseArea(curSrc,sizeLevel*curSrc.max())

    maxInd = NP.unravel_index(NP.argmax(curSrc),curSrc.shape)
    src1RcpN[ff] = distanceMeanFrom(curSrc,maxInd,1)
    src2RcpN[ff] = distanceMeanFrom(curSrc,maxInd,2)
    src3RcpN[ff] = distanceMeanFrom(curSrc,maxInd,3)
    src4RcpN[ff] = distanceMeanFrom(curSrc,maxInd,4)

    curSrc = srcLcpN[ff,srcYXN[0]:srcYXN[0] + srcDY,srcYXN[1]:srcYXN[1] + srcDX]
    srcLcpNSize[ff] = calcEllipseArea(curSrc,sizeLevel*curSrc.max())

    maxInd = NP.unravel_index(NP.argmax(curSrc),curSrc.shape)
    src1LcpN[ff] = distanceMeanFrom(curSrc,maxInd,1)
    src2LcpN[ff] = distanceMeanFrom(curSrc,maxInd,2)
    src3LcpN[ff] = distanceMeanFrom(curSrc,maxInd,3)
    src4LcpN[ff] = distanceMeanFrom(curSrc,maxInd,4)

fig = PL.figure()
pl = fig.subplots()
fig.suptitle('S source')
pl.plot(freqs,src0LcpS,'-',color='blue')
pl.plot(freqs,src1LcpS,'--',color='blue',lw=.5)
pl.plot(freqs,src2LcpS,'--',color='blue',lw=.5)
pl.plot(freqs,src3LcpS,'--',color='blue',lw=.5)
pl.plot(freqs,src4LcpS,'--',color='blue',lw=.5)
pl.plot(freqs,src0RcpS,'-',color='red')
pl.plot(freqs,src1RcpS,'--',color='red',lw=.5)
pl.plot(freqs,src2RcpS,'--',color='red',lw=.5)
pl.plot(freqs,src3RcpS,'--',color='red',lw=.5)
pl.plot(freqs,src4RcpS,'--',color='red',lw=.5)
pl.set_ylabel('K')
pl.set_ylim(0,5e5)
pl.set_xlabel('GHz')
pl.grid()

tpl = pl.twinx()
tpl.set_ylabel('arcsec')
tpl.plot(freqs,srcLcpSSize,'o-',color='skyblue')
tpl.plot(freqs,srcRcpSSize,'o-',color='tomato')

fig = PL.figure()
pl = fig.subplots()
fig.suptitle('N source')
pl.plot(freqs,src0LcpN,'-',color='blue')
pl.plot(freqs,src1LcpN,'--',color='blue',lw=.5)
pl.plot(freqs,src2LcpN,'--',color='blue',lw=.5)
pl.plot(freqs,src3LcpN,'--',color='blue',lw=.5)
pl.plot(freqs,src4LcpN,'--',color='blue',lw=.5)
pl.plot(freqs,src0RcpN,'-',color='red')
pl.plot(freqs,src1RcpN,'--',color='red',lw=.5)
pl.plot(freqs,src2RcpN,'--',color='red',lw=.5)
pl.plot(freqs,src3RcpN,'--',color='red',lw=.5)
pl.plot(freqs,src4RcpN,'--',color='red',lw=.5)
pl.set_ylabel('K')
pl.set_ylim(0,5e5)
pl.set_xlabel('GHz')
pl.grid()

tpl = pl.twinx()
tpl.set_ylabel('arcsec')
tpl.plot(freqs,srcLcpNSize,'o-',color='skyblue')
tpl.plot(freqs,srcRcpNSize,'o-',color='tomato')

#-------------------------------------------------------------------------------------------------------
arcsecPerPixel = NP.abs(clnFits[15][0].header['CDELT1'] * 3600)
N = 1440
H = int(240/arcsecPerPixel)
XY0 = int(clnFits[15][0].header['CRPIX1'])
R0 = int(930/arcsecPerPixel)
overLimb = NP.zeros((16,N,H))
hfOverLimb = NP.zeros((16,N,H))

alpha = 2*NP.pi*NP.linspace(0,N,N)/N
for ff in range(16):
    for radius in range(H):
        indX = XY0 + NP.array(NP.ceil((radius + R0)*NP.cos(alpha)),dtype='int')
        indY = XY0 + NP.array(NP.ceil((radius + R0)*NP.sin(alpha)),dtype='int')
        overLimb[ff,:,radius] = 0.5*(lcpImages[ff][indY,indX] + rcpImages[ff][indY,indX])

for ff in range(16):
    hfOverLimb[ff] = overLimb[15]
    
PL.figure()
PL.imshow(NP.concatenate(overLimb,axis=1).T,origin='lower',vmin=0,vmax=3e4)
PL.contour(NP.concatenate(hfOverLimb,axis=1).T,origin='lower',levels=[3e4])

saveFitsOverLimbHdu = fits.PrimaryHDU(header=phaseEdit.srhFits.hduList[0].header)
saveFitsOverLimbPath = 'srh_0612_overlimb_' + phaseEdit.srhFits.dateObs + '.fit'

overLimbFreqColumn = fits.Column(name='frequency', format='D', array = phaseEdit.srhFits.freqList)
saveFitsOverLimbFreqExtHdu = fits.BinTableHDU.from_columns([overLimbFreqColumn])

overLimbDataColumn = fits.Column(name='data', format=('%dD'%(N*H)), array = overLimb)
saveFitsOverLimbDataExtHdu = fits.BinTableHDU.from_columns([overLimbDataColumn])

hduList = fits.HDUList([saveFitsOverLimbHdu, saveFitsOverLimbFreqExtHdu, saveFitsOverLimbDataExtHdu])
hduList.writeto(saveFitsOverLimbPath)
#-------------------------------------------------------------------------------------------------------
iSpectrum = NP.zeros((1024,1024,3))
iSpectrum[:,:,0] = (lcpImages+rcpImages)[0:5].mean(axis=0)
iSpectrum[:,:,1] = (lcpImages+rcpImages)[5:10].mean(axis=0)
iSpectrum[:,:,2] = (lcpImages+rcpImages)[10:16].mean(axis=0)
iSpectrum /= iSpectrum.max()

vSpectrum = NP.zeros((1024,1024,3))
vSpectrum[:,:,0] = (lcpImages-rcpImages)[0:5].mean(axis=0)
vSpectrum[:,:,1] = (lcpImages-rcpImages)[5:10].mean(axis=0)
vSpectrum[:,:,2] = (lcpImages-rcpImages)[10:16].mean(axis=0)
vSpectrum /= vSpectrum.max()
vSpectrum /= 2
vSpectrum += .5

fig = PL.figure(figsize=(10,10))
fig.suptitle(phaseEdit.srhFits.dateObs)
pl = fig.subplots(2,2)
PL.tight_layout()
pl[0,0].imshow(iSpectrum[:,:,0],origin='lower',cmap='Reds',vmin=0.0,vmax=0.3)
pl[0,0].axis('off')
pl[0,1].imshow(iSpectrum[:,:,1],origin='lower',cmap='Greens',vmin=0.0,vmax=0.3)
pl[0,1].axis('off')
pl[1,0].imshow(iSpectrum[:,:,2],origin='lower',cmap='Blues',vmin=0.0,vmax=0.3)
pl[1,0].axis('off')
pl[1,1].imshow(iSpectrum + .2,origin='lower',vmin=0.2,vmax=0.01)
pl[1,1].axis('off')


fig = PL.figure(figsize=(10,10))
fig.suptitle(phaseEdit.srhFits.dateObs)
pl = fig.subplots(2,2)
PL.tight_layout()
pl[0,0].imshow(vSpectrum[:,:,0],origin='lower',cmap='Reds',vmin=0.0,vmax=0.6)
pl[0,0].axis('off')
pl[0,1].imshow(vSpectrum[:,:,1],origin='lower',cmap='Greens',vmin=0.0,vmax=0.6)
pl[0,1].axis('off')
pl[1,0].imshow(vSpectrum[:,:,2],origin='lower',cmap='Blues',vmin=0.0,vmax=0.6)
pl[1,0].axis('off')
pl[1,1].imshow(vSpectrum,origin='lower',vmin=0.0,vmax=1.0)
pl[1,1].axis('off')

