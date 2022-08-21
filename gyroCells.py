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

sdo171 = fits.open('../sunpy/data/aia20220215_033000_0171.fits')
sdo304 = fits.open('../sunpy/data/aia20220215_033000_0304.fits')
sdoHmi = fits.open('../sunpy/data/hmi_m_45s_2022_02_15_03_31_30_tai_magnetogram.fits')

hmiData = NP.flipud(NP.fliplr(NP.array(Image.fromarray(sdoHmi[1].data).resize((1024,1024)))))

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
pl = create_pl('20220215 03:30 SDO 171, SRH 7.8 I',MJ=128)
pl.imshow(sdo171[0].data,origin='lower',cmap='jet',vmin=0,vmax=1000)
pl.contour(iMeans[5],cmap='Greens_r',levels=[5e3,1.8e4,2.0e4,2.4e4,3.0e4,3.2e4,3.4e4,5e4])

pl = create_pl('20220215 03:30 SDO 171, SRH 11.8 I,V')
pl.imshow(sdo171[0].data,origin='lower',cmap='jet',vmin=0,vmax=2000)
pl.contourf(iMeans[15],cmap='Greens_r',levels=[2.5e4,3e4,4e4,5.e4,6e4,8e4,12e4])
pl.contour(vMeans[15],cmap='bwr',levels=[-2.e4,-5.e3,5.0e3,2.0e4],linestyles='--')
pl.set_xlim(620,720)
pl.set_ylim(350,450)

pl = create_pl('20220215 03:30 SDO 171, SRH 11.8 I,V')
pl.imshow(sdo171[0].data,origin='lower',cmap='jet',vmin=0,vmax=4000)
pl.contour(iMeans[15],cmap='gray',levels=[2.5e4,3e4,4e4,5.e4,6e4,8e4,12e4])
pl.contour(vMeans[15],cmap='bwr',levels=[-2.e4,-8.e3,8.0e3,2.0e4],linestyles='--')
pl.set_xlim(620,720)
pl.set_ylim(350,450)

#------------------------------------------------------------------------------
pl = create_pl('20220215 03:30 SDO 304, SRH 7.8 I',MJ=128)
pl.imshow(sdo304[0].data,origin='lower',cmap='jet',vmin=1,vmax=50)
pl.contour(iMeans[5],cmap='Greens_r',levels=[3e4,4e4,5.e4,8e4,10e4,15e4])

pl = create_pl('20220215 03:30 SDO 304, SRH 11.8 I,V')
pl.imshow(sdo304[0].data,origin='lower',cmap='jet',vmin=0,vmax=50)
pl.contourf(iMeans[15],cmap='Greens_r',levels=[2e4,3e4,4e4,5.e4,6e4,8e4,12e4])
pl.contour(vMeans[15],cmap='bwr',levels=[-2.e4,-5.e3,5.0e3,2.0e4],linestyles='--')
pl.set_xlim(620,720)
pl.set_ylim(350,450)

pl = create_pl('20220215 03:30 SDO 304, SRH 11.8 I,V')
pl.imshow(sdo304[0].data,origin='lower',cmap='jet',vmin=0,vmax=200)
pl.contour(iMeans[15],cmap='gray',levels=[2e4,3e4,4e4,5.e4,6e4,8e4,12e4])
pl.contour(vMeans[15],cmap='bwr',levels=[-2.e4,-5.e3,5.0e3,2.0e4],linestyles='--')
pl.set_xlim(620,720)
pl.set_ylim(350,450)
#------------------------------------------------------------------------------

FOV_AIA = float(sdo304[0].header['CDELT1']*1023)
FOV_HMI = float(sdoHmi[1].header['CDELT1']*4095)
FOV_SRH = 4.911*511

dSdo = int(1024*FOV_HMI/FOV_AIA/2 + .5)
dSrh = int(1024*FOV_HMI/FOV_SRH/2 + .5)

new304 = NP.array(Image.fromarray(sdo304[0].data[512-dSdo:512+dSdo,512-dSdo:512+dSdo]).resize((1024,1024)))
new171 = NP.array(Image.fromarray(sdo171[0].data[512-dSdo:512+dSdo,512-dSdo:512+dSdo]).resize((1024,1024)))
newSrhI = NP.array(Image.fromarray(iMeans[13,512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))
newSrhV = NP.array(Image.fromarray(vMeans[13,512-dSrh:512+dSrh,512-dSrh:512+dSrh]).resize((1024,1024)))

pl = create_pl('20220215 03:30 SDO 171, SRH 11.0 I',MJ=128)
#pl.imshow(hmiData,origin='lower',cmap='gray',vmin=-300,vmax=300)
#pl.imshow(new304,origin='lower',cmap='jet',vmin=1,vmax=70)
pl.imshow(new171,origin='lower',cmap='jet',vmin=1,vmax=1000)
pl.contour(newSrhI,cmap='Greens_r',levels=[3e4,5e4,8e4,12e4,15e4,20e4])
pl.contour(newSrhV,cmap='bwr',levels=[-4.e4,-2.e4,2.0e4,4.0e4],linestyles='--')
pl.set_xlim(650,750)
pl.set_ylim(340,440)
