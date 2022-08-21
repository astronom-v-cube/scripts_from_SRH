#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 07:23:57 2022

@author: sergey_lesovoi
"""

import numpy as NP
import pylab as PL
from skimage import measure#, resize
from astropy.io import fits
from scipy import interpolate
import scipy.optimize as opt
from PIL import Image

def create_pl(title,MJ=64):
    fig = PL.figure(figsize=(5,5))
    fig.tight_layout()
    fig.suptitle(title)
    pl = fig.subplots(nrows=1,ncols=1)
    pl.xaxis.set_major_locator(PL.MultipleLocator(MJ))
    pl.xaxis.set_minor_locator(PL.MultipleLocator(MJ//8))
    pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format))
    pl.yaxis.set_major_locator(PL.MultipleLocator(MJ))
    pl.xaxis.set_minor_locator(PL.MultipleLocator(MJ//8))
    pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format))
    pl.set_xlabel('arcmin')
    pl.set_ylabel('arcmin')
    pl.grid(linestyle='--')
    return pl

def arcmin_format(xy, pos):
  return '%2d' % ((xy - 1023/2) * 2.45 / 60);

lcpImagesList = []
rcpImagesList = []

scanNumber = phaseEdit.srhFits.freqTime.shape[1]
average = 3
for ss in NP.arange(0,scanNumber,3):
    print('Scan %d'%(ss))
    
    lcpImages = NP.zeros((16,512,512))
    rcpImages = NP.zeros((16,512,512))

    for ff in range(16):
        print('Frequency %d'%(ff))
        phaseEdit.srhFits.setFrequencyChannel(ff) 
        phaseEdit.srhFits.getHourAngle(ss)
        phaseEdit.srhFits.vis2uv(ss,phaseCorrect=True,average=average)
        phaseEdit.srhFits.uv2lmImage()
        phaseEdit.srhFits.lm2Heliocentric()
        lcpImages[ff] = phaseEdit.srhFits.lcp
        rcpImages[ff] = phaseEdit.srhFits.rcp
    
    lcpImagesList.append(lcpImages)
    rcpImagesList.append(rcpImages)
    lcpImages = 0
    rcpImages = 0

lcpImagesList = NP.array(lcpImagesList)
rcpImagesList = NP.array(rcpImagesList)

r0 = 170
c0 = 235
dS = 35

lcpCells = NP.zeros((lcpImagesList.shape[0],16,dS,dS))
rcpCells = NP.zeros((lcpImagesList.shape[0],16,dS,dS))
for ss in range(lcpImagesList.shape[0]):
    lcpCells[ss] = lcpImagesList[ss,:,r0:r0+dS,c0:c0+dS]
    rcpCells[ss] = rcpImagesList[ss,:,r0:r0+dS,c0:c0+dS]

# PL.figure()
# for ss in range(lcpImagesList.shape[0]):
#     PL.imshow(lcpCells[ss,0,:,:].T,origin='lower',cmap='jet',vmin=1e3,vmax=1e5,interpolation='bessel')
#     PL.savefig('lcp%s.png'%ss)

meanLcpCells = lcpCells.mean(axis=0)
meanRcpCells = rcpCells.mean(axis=0)

# meanLcpCells[:,17,23] = 0
# meanLcpCells[:,17,19] = 0
# meanLcpCells[:,19,14] = 0
# meanLcpCells[:,15,5] = 0

# meanRcpCells[:,17,23] = 0
# meanRcpCells[:,17,19] = 0
# meanRcpCells[:,19,14] = 0
# meanRcpCells[:,15,5] = 0
    
fig = PL.figure(figsize=(4,8))
pl = fig.subplots(nrows=1,ncols=2)
pl[0].imshow(NP.concatenate(meanLcpCells,axis=0),aspect=.5,vmin=1.5e4,vmax=4e4,origin='lower',cmap='jet',interpolation='bessel')
pl[1].imshow(NP.concatenate(meanRcpCells,axis=0),aspect=.5,vmin=1.5e4,vmax=4e4,origin='lower',cmap='jet',interpolation='bessel')

fig = PL.figure(figsize=(8,8))
pl = fig.subplots(nrows=2,ncols=2)
pl[0,0].plot(meanLcpCells[:,17,23],color='blue',linewidth=1)
pl[0,0].plot(meanRcpCells[:,17,23],color='red',linewidth=1)
pl[0,0].set_ylim(0,1.5e5)

pl[0,1].plot(meanLcpCells[:,17,19],color='blue',linewidth=1)
pl[0,1].plot(meanRcpCells[:,17,19],color='red',linewidth=1)
pl[0,1].set_ylim(0,1.5e5)

pl[1,0].plot(meanLcpCells[:,19,14],color='blue',linewidth=1)
pl[1,0].plot(meanRcpCells[:,19,14],color='red',linewidth=1)
pl[1,0].set_ylim(0,1.5e5)

pl[1,1].plot(meanLcpCells[:,15,5],color='blue',linewidth=1)
pl[1,1].plot(meanRcpCells[:,15,5],color='red',linewidth=1)
pl[1,1].set_ylim(0,1.5e5)

sdo171 = fits.open('../sunpy/data/aia20220219_033000_0171.fits')
sdo304 = fits.open('../sunpy/data/aia20220219_033000_0304.fits')
lcpImagesMean = lcpImagesList.mean(axis=0)
rcpImagesMean = rcpImagesList.mean(axis=0)

iMeans = NP.zeros((16,1024,1024))
vMeans = NP.zeros((16,1024,1024))

for ff in range(16):
    iLcp = NP.array(Image.fromarray(lcpImagesMean[ff]).resize((1024,1024)))
    iRcp = NP.array(Image.fromarray(rcpImagesMean[ff]).resize((1024,1024)))
    iMeans[ff] = 0.5*(iRcp + iLcp)
    vMeans[ff] = 0.5*(iRcp - iLcp)
    
#------------------------------------------------------------------------------
pl = create_pl('20220219 03:30 SDO 171, SRH 7.8 I',MJ=128)
pl.imshow(sdo171[0].data,origin='lower',cmap='jet',vmin=0,vmax=1000)
pl.contour(iMeans[5],cmap='Greens_r',levels=[5e3,1.8e4,2.0e4,2.4e4,3.0e4,3.2e4,3.4e4,5e4])

pl = create_pl('20220219 03:30 SDO 171, SRH 11.8 I,V')
pl.imshow(sdo171[0].data,origin='lower',cmap='jet',vmin=0,vmax=2000)
pl.contourf(iMeans[15],cmap='Greens_r',levels=[1.8e4,2.0e4,2.4e4,3.0e4,4e4])
pl.contour(vMeans[15],cmap='bwr',levels=[-2.e3,-1.e3,1.0e3,2.0e3],linestyles='--')
pl.set_xlim(450,550)
pl.set_ylim(320,420)

pl = create_pl('20220219 03:30 SDO 171, SRH 11.8 I,V')
pl.imshow(sdo171[0].data,origin='lower',cmap='jet',vmin=0,vmax=4000)
pl.contour(iMeans[15],cmap='gray',levels=[1.8e4,2.0e4,2.4e4,3.0e4,4e4])
pl.contour(vMeans[15],cmap='bwr',levels=[-2.e3,-1.e3,1.0e3,2.0e3],linestyles='--')
pl.set_xlim(450,550)
pl.set_ylim(320,420)

#------------------------------------------------------------------------------
pl = create_pl('20220219 03:30 SDO 304, SRH 7.8 I',MJ=128)
pl.imshow(sdo304[0].data,origin='lower',cmap='jet',vmin=1,vmax=50)
pl.contour(iMeans[5],cmap='Greens_r',levels=[5e3,1.8e4,2.4e4,3.0e4,3.4e4,5e4])

pl = create_pl('20220219 03:30 SDO 304, SRH 11.8 I,V')
pl.imshow(sdo304[0].data,origin='lower',cmap='jet',vmin=0,vmax=50)
pl.contourf(iMeans[15],cmap='Greens_r',levels=[1.8e4,2.0e4,2.4e4,3.0e4,4e4])
pl.contour(vMeans[15],cmap='bwr',levels=[-2.e3,-1.e3,1.0e3,2.0e3],linestyles='--')
pl.set_xlim(450,550)
pl.set_ylim(320,420)

pl = create_pl('20220219 03:30 SDO 304, SRH 11.8 I,V')
pl.imshow(sdo304[0].data,origin='lower',cmap='jet',vmin=0,vmax=200)
pl.contour(iMeans[15],cmap='gray',levels=[1.8e4,2.0e4,2.4e4,3.0e4,4e4])
pl.contour(vMeans[15],cmap='bwr',levels=[-2.e3,-1.e3,1.0e3,2.0e3],linestyles='--')
pl.set_xlim(450,550)
pl.set_ylim(320,420)
#------------------------------------------------------------------------------
fig = PL.figure(figsize=(4,8))
fig.suptitle('SRH 20220219 03:30 LCP, RCP 5.8-11.8 GHz')
pl = fig.subplots(nrows=1,ncols=2)
pl[0].imshow(NP.concatenate(lcpCells.mean(axis=0),axis=0),aspect=.5,vmin=1.5e4,vmax=4e4,origin='lower',cmap='jet')
pl[1].imshow(NP.concatenate(rcpCells.mean(axis=0),axis=0),aspect=.5,vmin=1.5e4,vmax=4e4,origin='lower',cmap='jet')

#------------------------------------------------------------------------------
# r0 = 275
# c0 = 165
# dS = 35
r0 = 300
c0 = 55
dS = 45

lcpCells = NP.zeros((lcpImagesList.shape[0],16,dS,dS))
rcpCells = NP.zeros((lcpImagesList.shape[0],16,dS,dS))
for ss in range(lcpImagesList.shape[0]):
    lcpCells[ss] = lcpImagesList[ss,:,r0:r0+dS,c0:c0+dS]
    rcpCells[ss] = rcpImagesList[ss,:,r0:r0+dS,c0:c0+dS]

fig = PL.figure(figsize=(4,8))
fig.suptitle('SRH 20220219 camp 03:30 LCP, RCP 5.8-11.8 GHz')
pl = fig.subplots(nrows=1,ncols=2)
pl[0].imshow(NP.concatenate(lcpCells.mean(axis=0),axis=0),aspect=.5,vmin=1.5e4,vmax=2e4,origin='lower',cmap='jet')
pl[1].imshow(NP.concatenate(rcpCells.mean(axis=0),axis=0),aspect=.5,vmin=1.5e4,vmax=2e4,origin='lower',cmap='jet')

lcpCellsMean_t = lcpCells.mean(axis=0)
rcpCellsMean_t = rcpCells.mean(axis=0)
lcpCellsMean_f = lcpCells.mean(axis=1)
rcpCellsMean_f = rcpCells.mean(axis=1)

lcpCellMaxT_f = []
rcpCellMaxT_f = []
lcpCellMeanT_f = []
rcpCellMeanT_f = []
for ff in range(16):
    maxInd = NP.unravel_index(NP.argmax(lcpCellsMean_t[ff]),lcpCellsMean_t[ff].shape)
    lcpCellMaxT_f.append(lcpCellsMean_t[ff,maxInd[0]-1:maxInd[0]+2,maxInd[1]-1:maxInd[1]+2].mean())

    maxInd = NP.unravel_index(NP.argmax(rcpCellsMean_t[ff]),rcpCellsMean_t[ff].shape)
    rcpCellMaxT_f.append(rcpCellsMean_t[ff,maxInd[0]-1:maxInd[0]+2,maxInd[1]-1:maxInd[1]+2].mean())

    lcpCellMeanT_f.append(lcpCellsMean_t[ff,0:5,0:25].mean())
    rcpCellMeanT_f.append(rcpCellsMean_t[ff,0:5,0:25].mean())

lcpCellMaxT_f = NP.array(lcpCellMaxT_f)
rcpCellMaxT_f = NP.array(rcpCellMaxT_f)
lcpCellMeanT_f = NP.array(lcpCellMeanT_f)
rcpCellMeanT_f = NP.array(rcpCellMeanT_f)

lcpCellMaxT_t = []
rcpCellMaxT_t = []
lcpCellMeanT_t = []
rcpCellMeanT_t = []
for tt in range(lcpCellsMean_f.shape[0]):
    maxInd = NP.unravel_index(NP.argmax(lcpCellsMean_f[tt]),lcpCellsMean_f[tt].shape)
    lcpCellMaxT_t.append(lcpCellsMean_f[tt,maxInd[0]-1:maxInd[0]+2,maxInd[1]-1:maxInd[1]+2].mean())

    maxInd = NP.unravel_index(NP.argmax(rcpCellsMean_f[tt]),rcpCellsMean_f[tt].shape)
    rcpCellMaxT_t.append(rcpCellsMean_f[tt,maxInd[0]-1:maxInd[0]+2,maxInd[1]-1:maxInd[1]+2].mean())

    lcpCellMeanT_t.append(lcpCellsMean_f[tt,0:3,0:3].mean())
    rcpCellMeanT_t.append(rcpCellsMean_f[tt,0:3,0:3].mean())

lcpCellMaxT_t = NP.array(lcpCellMaxT_t)
rcpCellMaxT_t = NP.array(rcpCellMaxT_t)
lcpCellMeanT_t = NP.array(lcpCellMeanT_t)
rcpCellMeanT_t = NP.array(rcpCellMeanT_t)

campSize = []
waveLen = 3e8/(phaseEdit.srhFits.freqList*1e3)

for ff in range(16):
    curI = lcpCellsMean_t[ff] + rcpCellsMean_t[ff]
    curMax = curI.max()
    half = NP.where(curI[12:22,12:22] > 0.9*curMax)
#    half = NP.where(curI > 0.98*curMax)
    campSize.append((NP.sqrt(half[0].shape[0]) - 1.)*4.911)
    curI = 0
                   