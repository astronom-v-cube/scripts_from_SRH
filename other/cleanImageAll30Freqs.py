#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 01:37:52 2021

@author: svlesovoi
"""

import numpy as NP
import pylab as PL
import srh36Utils
from astropy.io import fits
from scipy.stats import linregress
from skimage.transform import warp, AffineTransform
import scipy
from ZirinTb import ZirinTb

def arcmin_format(xy, pos):
  return '%2d' % ((xy - 512/2) * 4.911 / 60);

def drawRLPointSpectrum(rSpectrum, lSpectrum, x, y, freqs, label, ylimit):
    fig = PL.figure(facecolor='black',edgecolor='black')
    fig.subplots_adjust(left=0.1,right=0.95,top=0.95,bottom=0.05)
    ax = PL.axes()
    ax.set_facecolor('black')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['bottom'].set_color('white')
    PL.title(label, color='white')
    PL.xlabel('GHz')
    PL.ylabel('brightness temeperature x 1000')
    PL.ylim(ylimit[0],ylimit[1])
    PL.grid(linestyle='--')
    PL.plot(1e-9*freqs,1e-3*0.5*rSpectrum[:,y,x],'.',color='red',label='RCP')
    PL.plot(1e-9*freqs,1e-3*0.5*lSpectrum[:,y,x],'.',color='dodgerblue',label='LCP')
    PL.plot(1e-9*freqs,1e-3*0.5*(rSpectrum[:,y,x] - lSpectrum[:,y,x]),'-',color='yellow',label='V')
    PL.legend()
    PL.savefig('/home/svlesovoi/Pictures/VAK2021/spectrum_%s.png'%label, facecolor='black', edgecolor='black')

    
iFitses = srh36Utils.findFits('SRH36_temp_20210815_30freq/clean/','*I*.fit')
vFitses = srh36Utils.findFits('SRH36_temp_20210815_30freq/clean/','*V*.fit')
iFitses.sort()
vFitses.sort()
iFitses = iFitses#[0:13]
vFitses = vFitses#[0:13]

iImages = []
vImages = []
iT0 = 1e5
vT0 = 1e4

for fitsName in iFitses:
    fitsHandle = fits.open(fitsName)
    iImages.append(fitsHandle[0].data)
    fitsHandle.close()
for fitsName in vFitses:
    fitsHandle = fits.open(fitsName)
    vImages.append(fitsHandle[0].data)
    fitsHandle.close()
    
iImages = NP.array(iImages)
vImages = NP.array(vImages)
rImages = iImages + vImages
lImages = iImages - vImages
rImages[15] *= -1
rImages[16] *= -1

r15Mean = rImages[15].mean()
rImages[15] = (rImages[15] - r15Mean)*20 + r15Mean
r16Mean = rImages[16].mean()
rImages[16] = (rImages[16] - r16Mean)*45 + r16Mean

rImages[15] = NP.roll(rImages[15],(-327+390,-257+208),axis=(0,1))
rImages[16] = NP.roll(rImages[16],(-401+462,-309+320),axis=(0,1))

iImages = (rImages + lImages).copy()
vImages = (rImages - lImages).copy()

iHist = []
for i in range(len(iFitses)):
    iI = iImages[i].copy()
    iHist.append(NP.histogram(iImages[i],500))

skySunInd = []
PL.figure()
for i in range(len(iFitses)): 
    PL.plot(iHist[i][1][1:],iHist[i][0] + i*10000)
    peaks, _ = scipy.signal.find_peaks(iHist[i][0],prominence = 400)
    ind1 = peaks[0]
    skySunInd.append(peaks)
    PL.plot(iHist[i][1][peaks+1],iHist[i][0][peaks] + i*10000,'.',color='green')
    
sunStair = []
skyLevel = []
for i in range(len(iFitses)): 
    if(skySunInd[i].shape[0] > 1):
        sunStair.append([i,iHist[i][1][skySunInd[i][1]] - iHist[i][1][skySunInd[i][0]]])
        if(i > 15):
            skyLevel.append([i,iHist[i][1][skySunInd[i][0]] - sunStair[-1][1]])
        else:
            skyLevel.append([i,iHist[i][1][skySunInd[i][0]]])

sunStair = NP.array(sunStair)
skyLevel = NP.array(skyLevel)
PL.figure()
PL.plot(sunStair[:,0],sunStair[:,1],'.')
PL.plot(skyLevel[:,0],skyLevel[:,1],'.')
levelsArg = NP.linspace(0,len(iFitses)-1,len(iFitses))

slope,inter,A,B,C = linregress(sunStair[:,0],sunStair[:,1])
sunCalLevels = slope*levelsArg + inter
PL.plot(levelsArg,sunCalLevels,'x')

slope,inter,A,B,C = linregress(skyLevel[:,0],skyLevel[:,1])
skyCalLevels = slope*levelsArg + inter
PL.plot(levelsArg,skyCalLevels,'x')

ZirinQSunTb = ZirinTb()
freqArg = NP.linspace(0,len(iFitses),len(iFitses))*0.1e9 + 2.8e9

for i in range(len(iFitses)): 
    iImages[i] -= skyCalLevels[i]
    iImages[i] /= sunCalLevels[i]
    iImages[i] *= ZirinQSunTb.getTbAtFrequncy(freqArg[i]*1e-9)*1e3
    vImages[i] /= sunCalLevels[i]
    vImages[i] *= ZirinQSunTb.getTbAtFrequncy(freqArg[i]*1e-9)*1e3
      
iMeanImage = iImages.mean(axis=0)
vMeanImage = vImages.mean(axis=0)

vT0 = 3e3
maxV = 1.*NP.max(NP.abs(vMeanImage))
nLevels = NP.linspace(-maxV,-vT0,5)
pLevels = NP.linspace(vT0,maxV,5)

PL.figure()
PL.imshow(iMeanImage,origin='lower',cmap='hot',vmin=0,vmax=iT0)
PL.contour(vMeanImage,levels=nLevels,colors='blue',linewidths=0.5)
PL.contour(vMeanImage,levels=pLevels,colors='red',linewidths=0.5)

PL.figure()
x0 = 50
y0 = 50
dL = 350
iImages[14] = NP.roll(iImages[14],[-100,0],axis=(0,1))
vImages[14] = NP.roll(vImages[14],[-100,0],axis=(0,1))
iImages[15] = NP.roll(iImages[15],[-100,0],axis=(0,1))
vImages[15] = NP.roll(vImages[15],[-100,0],axis=(0,1))
iImages[16] = NP.roll(iImages[16],[100,0],axis=(0,1))
vImages[16] = NP.roll(vImages[16],[100,0],axis=(0,1))
iImages[17] = NP.roll(iImages[17],[-100,0],axis=(0,1))
vImages[17] = NP.roll(vImages[17],[-100,0],axis=(0,1))
iImages[18] = NP.roll(iImages[18],[-100,0],axis=(0,1))
vImages[18] = NP.roll(vImages[18],[-100,0],axis=(0,1))

maxInd0 = NP.unravel_index(NP.argmax(iImages[0,y0:y0+dL,x0:x0+dL]),(dL,dL))
maxInd = []
iShifted = NP.zeros_like(iImages)
vShifted = NP.zeros_like(vImages)
for i in range(len(iFitses)):
    maxInd.append(NP.unravel_index(NP.argmax(iImages[i][y0:y0+dL,x0:x0+dL]),(dL,dL)))
    if (i < 20):
        if (maxInd[-1][0] > 200):
            iShifted[i] = NP.roll(iImages[i],NP.array(maxInd0) - NP.array(maxInd[-1]),axis=(0,1))
            vShifted[i] = NP.roll(vImages[i],NP.array(maxInd0) - NP.array(maxInd[-1]),axis=(0,1))
        else:
            print(i)
            iShifted[i] = iImages[i]
            vShifted[i] = vImages[i]
    else:
        iShifted[i] = iImages[i]
        vShifted[i] = vImages[i]

iShifted[14] = NP.roll(iShifted[14],[-13,-62],axis=(0,1))
vShifted[14] = NP.roll(vShifted[14],[-13,-62],axis=(0,1))

for i in range(len(iFitses)):
    maxI = iShifted[i].max()
    PL.contourf(iShifted[i],levels=[0.5*maxI,0.51*maxI],origin='lower',colors=[(i/len(iFitses)/2,i/len(iFitses)/2,i/len(iFitses)/2)])

shiftedMaxInd = []
for i in range(len(iFitses)):
    shiftedMaxInd.append(NP.unravel_index(NP.argmax(iShifted[i][y0:y0+dL,x0:x0+dL]),(dL,dL)))

spectrumIslope = NP.zeros_like(iMeanImage)
spectrumIinter = NP.zeros_like(iMeanImage)
spectrumVslope = NP.zeros_like(iMeanImage)
spectrumVinter = NP.zeros_like(iMeanImage)
#for i in range(spectrumIslope.shape[0]):
#    for j in range(spectrumIslope.shape[1]):
#        slope, inter, A, B, C = linregress(freqArg,iImages[:,i,j])
#        spectrumIslope[i,j] = slope
#        spectrumIinter[i,j] = inter
#        slope, inter, A, B, C = linregress(freqArg,vImages[:,i,j])
#        spectrumVslope[i,j] = slope
#        spectrumVinter[i,j] = inter

rMean = 0.95*(iShifted.mean(axis=0) - vShifted.mean(axis=0))
lMean = iShifted.mean(axis=0) + vShifted.mean(axis=0)
iMean = lMean + rMean
vMean = lMean - rMean

iT0 = 100000
vT0 = 3e3
maxV = 1.*NP.max(NP.abs(vMeanImage))
nLevels = NP.linspace(-maxV,-vT0,5)
pLevels = NP.linspace(vT0,maxV,5)

PL.figure()
PL.imshow(iMean,origin='lower',cmap='hot',vmin=3e4,vmax=iT0)
PL.contour(vMean,levels=nLevels,colors='blue',linewidths=0.5)
PL.contour(vMean,levels=pLevels,colors='red',linewidths=0.5)

nRows = 3
nCols = 10
fig, pl = PL.subplots(nrows=nRows,ncols=nCols,figsize=(12,4))
fig.subplots_adjust(left=0.01,right=0.99,top=0.99,bottom=0.01,hspace=0.01,wspace=0.01)
for row in range(nRows):
    for col in range(nCols):
        pl[row,col].axis('off')
        pl[row,col].imshow(iShifted[row*nCols + col][275:325,240:290],origin='lower',vmin=0,vmax=2e5,cmap='hot')
        pl[row,col].contour(iShifted[row*nCols + col][275:325,240:290],origin='lower',levels=[.9e5,1.2e5,1.6e5,2.3e5],colors='black',linewidths=0.5)
        pl[row,col].text(3,43,'%1.1f GHz'%(freqArg[row*6 + col]*1e-9),color='white')

fig, pl = PL.subplots(nrows=nRows,ncols=nCols,figsize=(12,4))
fig.subplots_adjust(left=0.01,right=0.99,top=0.99,bottom=0.01,hspace=0.01,wspace=0.01)
for row in range(nRows):
    for col in range(nCols):
        pl[row,col].axis('off')
        pl[row,col].imshow(vShifted[row*nCols + col][275:325,240:290],origin='lower',vmin=-2e4,vmax=2e4,cmap='gray')
        pl[row,col].contour(iShifted[row*nCols + col][275:325,240:290],origin='lower',levels=[.9e5,1.2e5],colors='white',linewidths=0.5)

rShifted = iShifted + vShifted
lShifted = iShifted - vShifted

xSize = rShifted.shape[1]
ySize = rShifted.shape[2]

rSpectrumResult = NP.zeros((xSize,ySize,3))
lSpectrumResult = NP.zeros((xSize,ySize,3))
iSpectrum = NP.zeros((xSize,ySize,3))
vSpectrum = NP.zeros((xSize,ySize,3))
iSpectrumSqrt = NP.zeros((xSize,ySize,3))
vSpectrumSqrt = NP.zeros((xSize,ySize,3))
iSpectrumNorm = NP.zeros((xSize,ySize,3))
vSpectrumNorm = NP.zeros((xSize,ySize,3))
for i in range(ySize):
    for j in range(xSize):
        spectrum = rShifted[:,i,j].copy()
        spectrum /= spectrum.max()
        rSpectrumResult[i,j,0] = NP.sqrt(NP.abs(spectrum[0:10].mean()))
        rSpectrumResult[i,j,1] = NP.sqrt(NP.abs(spectrum[10:20].mean()))
        rSpectrumResult[i,j,2] = NP.sqrt(NP.abs(spectrum[20:30].mean()))

        spectrum = lShifted[:,i,j].copy()
        spectrum /= spectrum.max()
        lSpectrumResult[i,j,0] = NP.sqrt(NP.abs(spectrum[0:10].mean()))
        lSpectrumResult[i,j,1] = NP.sqrt(NP.abs(spectrum[10:20].mean()))
        lSpectrumResult[i,j,2] = NP.sqrt(NP.abs(spectrum[20:30].mean()))

        spectrum = NP.abs(iShifted[:,i,j].copy())
        iSpectrum[i,j,0] = spectrum[0:10].mean()
        iSpectrum[i,j,1] = spectrum[10:20].mean()
        iSpectrum[i,j,2] = spectrum[20:30].mean()
        iSpectrumSqrt[i,j,0] = NP.sqrt(spectrum[0:10].mean())
        iSpectrumSqrt[i,j,1] = NP.sqrt(spectrum[10:20].mean())
        iSpectrumSqrt[i,j,2] = NP.sqrt(spectrum[20:30].mean())
#        spectrum /= spectrum.max()
#        iSpectrumNorm[i,j,0] = spectrum[0:10].mean()
#        iSpectrumNorm[i,j,1] = spectrum[10:20].mean()
#        iSpectrumNorm[i,j,2] = spectrum[20:30].mean()

        spectrum[0:12] /= spectrum[0:12].max()
        iSpectrumNorm[i,j,0] = spectrum[0:4].mean()
        iSpectrumNorm[i,j,1] = spectrum[4:8].mean()
        iSpectrumNorm[i,j,2] = spectrum[8:12].mean()

        spectrum = NP.abs(vShifted[:,i,j].copy())
        vSpectrum[i,j,0] = spectrum[0:10].mean()
        vSpectrum[i,j,1] = spectrum[10:20].mean()
        vSpectrum[i,j,2] = spectrum[20:30].mean()
        vSpectrumSqrt[i,j,0] = NP.sqrt(spectrum[0:10].mean())
        vSpectrumSqrt[i,j,1] = NP.sqrt(spectrum[10:20].mean())
        vSpectrumSqrt[i,j,2] = NP.sqrt(spectrum[20:30].mean())
#        spectrum /= spectrum.max()
#        vSpectrumNorm[i,j,0] = spectrum[0:10].mean()
#        vSpectrumNorm[i,j,1] = spectrum[10:20].mean()
#        vSpectrumNorm[i,j,2] = spectrum[20:30].mean()

        spectrum[0:12] /= spectrum[0:12].max()
        vSpectrumNorm[i,j,0] = spectrum[0:4].mean()
        vSpectrumNorm[i,j,1] = spectrum[4:8].mean()
        vSpectrumNorm[i,j,2] = spectrum[8:12].mean()
    
rSpectrumResult /= rSpectrumResult.max()
lSpectrumResult /= lSpectrumResult.max()
iSpectrum /= iSpectrum.max()
vSpectrum /= vSpectrum.max()
iSpectrumSqrt /= iSpectrumSqrt.max()
vSpectrumSqrt /= vSpectrumSqrt.max()

rSpectrumResult = rSpectrumResult.clip(min = 0,max = 1)
lSpectrumResult = lSpectrumResult.clip(min = 0,max = 1)
iSpectrum = iSpectrum.clip(min = 0,max = 1)
vSpectrum = vSpectrum.clip(min = 0,max = 1)
iSpectrumSqrt = iSpectrum.clip(min = 0,max = 1)
vSpectrumSqrt = vSpectrum.clip(min = 0,max = 1)

fig, pl = PL.subplots(figsize=(8,8))
fig.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.1)
pl.xaxis.set_major_locator(PL.MultipleLocator(128));
pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
pl.yaxis.set_major_locator(PL.MultipleLocator(128));
pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
pl.set_title('SRH 20210815 07:30 UT multifrequency (2.8-5.8 GHz) I image')
pl.set_xlabel('arcmin')
pl.set_ylabel('arcmin')
pl.grid(linestyle='--')
pl.imshow(iShifted.mean(axis=0),origin='lower',cmap='hot',vmin=10000)

fig, pl = PL.subplots(figsize=(8,8))
fig.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.1)
pl.xaxis.set_major_locator(PL.MultipleLocator(128));
pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
pl.yaxis.set_major_locator(PL.MultipleLocator(128));
pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
pl.set_title('SRH 20210815 07:30 UT multifrequency (2.8-5.8 GHz) V image')
pl.set_xlabel('arcmin')
pl.set_ylabel('arcmin')
pl.grid(linestyle='--')
pl.imshow(vShifted.mean(axis=0),origin='lower',cmap='gray')

fig, pl = PL.subplots(figsize=(8,8))
fig.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.1)
pl.xaxis.set_major_locator(PL.MultipleLocator(128));
pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
pl.yaxis.set_major_locator(PL.MultipleLocator(128));
pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
pl.set_title('SRH 20210815 07:30 UT RGB spectrum I')
pl.set_xlabel('arcmin')
pl.set_ylabel('arcmin')
pl.grid(linestyle='--')
pl.imshow(iSpectrum,origin='lower')

fig, pl = PL.subplots(figsize=(8,8))
fig.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.1)
pl.xaxis.set_major_locator(PL.MultipleLocator(128));
pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
pl.yaxis.set_major_locator(PL.MultipleLocator(128));
pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
pl.set_title('SRH 20210815 07:30 UT RGB spectrum V')
pl.set_xlabel('arcmin')
pl.set_ylabel('arcmin')
pl.grid(linestyle='--')
pl.imshow(vSpectrum,origin='lower')

fig, pl = PL.subplots(figsize=(8,8))
fig.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.1)
pl.xaxis.set_major_locator(PL.MultipleLocator(128));
pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
pl.yaxis.set_major_locator(PL.MultipleLocator(128));
pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
pl.set_title('SRH 20210815 07:30 UT nonlinear RGB spectrum I')
pl.set_xlabel('arcmin')
pl.set_ylabel('arcmin')
pl.grid(linestyle='--')
pl.imshow(iSpectrumSqrt,origin='lower')

fig, pl = PL.subplots(figsize=(8,8))
fig.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.1)
pl.xaxis.set_major_locator(PL.MultipleLocator(128));
pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
pl.yaxis.set_major_locator(PL.MultipleLocator(128));
pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
pl.set_title('SRH 20210815 07:30 UT nonlinear RGB spectrum V')
pl.set_xlabel('arcmin')
pl.set_ylabel('arcmin')
pl.grid(linestyle='--')
pl.imshow(vSpectrumSqrt,origin='lower')

fig, pl = PL.subplots(figsize=(8,8))
fig.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.1)
pl.xaxis.set_major_locator(PL.MultipleLocator(128));
pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
pl.yaxis.set_major_locator(PL.MultipleLocator(128));
pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
pl.set_title('SRH 20210815 07:30 UT normalized RGB spectrum I')
pl.set_xlabel('arcmin')
pl.set_ylabel('arcmin')
pl.grid(linestyle='--')
pl.imshow(iSpectrumNorm,origin='lower')

fig, pl = PL.subplots(figsize=(8,8))
fig.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.1)
pl.xaxis.set_major_locator(PL.MultipleLocator(128));
pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
pl.yaxis.set_major_locator(PL.MultipleLocator(128));
pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
pl.set_title('SRH 20210815 07:30 UT normalized RGB spectrum V')
pl.set_xlabel('arcmin')
pl.set_ylabel('arcmin')
pl.grid(linestyle='--')
pl.imshow(vSpectrumNorm,origin='lower')

fig = PL.figure(facecolor='black',edgecolor='black')
ax = PL.axes()
ax.set_facecolor('black')
ax.xaxis.label.set_color('white')
ax.yaxis.label.set_color('white')
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')
ax.spines['left'].set_color('white')
PL.title('source A')
PL.xlabel('Hz')
PL.ylabel('brightness temeperature')
PL.ylim(0,2e5)
PL.grid(linestyle='--')
src = [294,266]#src A
dL = 1
PL.plot(freqArg,0.5*rShifted[:,src[0]-dL:src[0]+dL,src[1]].mean(axis=1),'.',color='red')
PL.plot(freqArg,0.5*lShifted[:,src[0]-dL:src[0]+dL,src[1]].mean(axis=1),'.',color='blue')
