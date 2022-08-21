#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 11:49:23 2021

@author: sergey_lesovoi
"""
#20211111 scan0 04:13:39 scans 300

import numpy as NP
import pylab as PL
import ZirinTb
from skimage import measure, resize
from astropy.io import fits

def freqGHzFormat(f, pos):
    return '%3.2f'%(f * 0.4 + 5.8)

def hms_format(t, pos):
#  t *= 10
#  t += cpTime[0,0]
  t = int(t)
  if (t < cpTime.shape[0]):
      if (t < 10):
          t = cpTime[int(t),0]
          hh = int(t / 3600.)
          t -= hh*3600.
          mm = int(t / 60.)
          t -= mm*60.
          ss = int(t)
          return '%02d:%02d:%02d' % (hh,mm,ss)
      else:
#          t = cpTime[int(t),0] - cpTime[0,0]
          t *= 10
          mm = int(t / 60.)
          t -= mm*60.
          ss = int(t)
          return '%02d:%02d' % (mm,ss)
  else:
      t *= 10
      t += cpTime[0,0]
      hh = int(t / 3600.)
      t -= hh*3600.
      mm = int(t / 60.)
      t -= mm*60.
      ss = int(t)
      return '%02d:%02d:%02d' % (hh,mm,ss)
#  return '%02d:%02d' % (hh,mm);

def calcEllipseArea(cell, level):
    contours = measure.find_contours(cell, level)
    contourLength = []
    for n, contour in enumerate(contours):
        contourLength.append(len(contour))
    contourLength = NP.array(contourLength)
    
    contourMaxInd = NP.argmax(contourLength)
    srcEll.estimate(contours[contourMaxInd])
    exc, eyc, ea, eb, theta = srcEll.params
    return NP.sqrt(NP.pi*ea*eb)


Zi = ZirinTb.ZirinTb()

#for ff in range(16):
#    print('Centering frequency %d'%(ff))
#    phaseEdit.onFrequencyChannelChanged(ff)
#    phaseEdit.onCenterDisk()

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
        lcpImages[ff] = NP.roll(phaseEdit.srhFits.lcp,(int(ss/3.5),int(-ss/30)),axis=(0,1))
        rcpImages[ff] = NP.roll(phaseEdit.srhFits.rcp,(int(ss/3.5),int(-ss/30)),axis=(0,1))
    
    lcpImagesList.append(lcpImages)
    rcpImagesList.append(rcpImages)
    lcpImages = 0
    rcpImages = 0

cpTime = phaseEdit.srhFits.freqTime[:,NP.arange(0,scanNumber,3)].T

lcpImagesList = NP.array(lcpImagesList,dtype='float32')
rcpImagesList = NP.array(rcpImagesList,dtype='float32')

S = 10
R = 6
arcTheta = NP.deg2rad(9)
diameter = lcpImagesList.shape[3]
arcDX = 20
arcDY = 1
cellMax = False

arcs = []

for radiusIndex in range(R):
    lcpMaxT_arc1_t = NP.zeros((S,lcpImagesList.shape[0], 16))
    rcpMaxT_arc1_t = NP.zeros((S,lcpImagesList.shape[0], 16))
    
    arc1radius = 210 + radiusIndex*15
    arc1theta = NP.deg2rad(115 - radiusIndex*1.2)
    arc1 = NP.zeros((S,2),dtype='int')
    for i in range(S):
        arc1[i,0] = int(diameter/2 + arc1radius * NP.sin((i - S/2)*arcTheta/S + arc1theta + .5))
        arc1[i,1] = int(diameter/2 + arc1radius * NP.cos((i - S/2)*arcTheta/S + arc1theta + .5))
    
    arcs.append(arc1)
    
    lcpArcStrip0_t = NP.zeros((S*16,lcpImagesList.shape[0]))
    rcpArcStrip0_t = NP.zeros((S*16,lcpImagesList.shape[0]))
    
    for ss in range(lcpImagesList.shape[0]):
        for ii in range(16):
            for p in range(S):
                cell = lcpImagesList[ss,ii,arc1[p,0]:arc1[p,0] + arcDY,arc1[p,1]:arc1[p,1] + arcDX]
                if (cellMax):
                    lcpArcStrip0_t[p*16+ii,ss] = cell.max()
                else:
                    lcpArcStrip0_t[p*16+ii,ss] = cell.mean()
                cell = rcpImagesList[ss,ii,arc1[p,0]:arc1[p,0] + arcDY,arc1[p,1]:arc1[p,1] + arcDX]
                if (cellMax):
                    rcpArcStrip0_t[p*16+ii,ss] = cell.max()
                else:
                    rcpArcStrip0_t[p*16+ii,ss] = cell.mean()
    #------------------------------------------------------------------------------------------------------------
    #PL.figure()
    #PL.title('arc1 %d'%(arc1radius))
    #PL.imshow(lcpImagesList[40,1],vmax=1e5,vmin=0)
    #PL.plot(arc1[:,1], arc1[:,0],'.',color='white')
    #------------------------------------------------------------------------------------------------------------
    Tsat = 5e4
    Tmin = 0e3
    cmap = 'rainbow'
    fig = PL.figure(figsize=(3,8))
    fig.suptitle('arc1 %d'%(arc1radius))
    fig.subplots_adjust(hspace=0.001,wspace=0.001,top=0.95,bottom=0.01,left=0.01,right=0.99)
    pl = fig.subplots(nrows=1,ncols=1)
    pl.xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
    pl.xaxis.set_major_locator(PL.MultipleLocator(18))
    pl.yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
    pl.yaxis.set_major_locator(PL.MultipleLocator(6))
    pl.set_xlabel('UT')
    pl.set_ylabel('GHz')
    pl.axis('off')
    pl.imshow(rcpArcStrip0_t + lcpArcStrip0_t,aspect=2., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap=cmap)
#    pl.imshow(rcpArcStrip0_t,aspect=2., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap=cmap)
    
#------------------------------------------------------------------------------------------------------------
PL.figure()
PL.imshow(rcpImagesList.mean(axis=(0,1)),vmax=5e3,vmin=0,origin='lower')
#PL.imshow(lcpImagesList.mean(axis=(0))[15],vmax=5e3,vmin=0,origin='lower')
for rr in range(R):
    PL.plot(arcs[rr][:,1],arcs[rr][:,0],'.',markersize=0.7)

#------------------------------------------------------------------------------------------------------------
cellMax = True
lcpArcPoint0_t = NP.zeros((16,lcpImagesList.shape[0]))
rcpArcPoint0_t = NP.zeros((16,lcpImagesList.shape[0]))
arcDX = 30
arcDY = 30

for ss in range(lcpImagesList.shape[0]):
    arc1radius = int(230 + ss*0.7)
    arc1theta = NP.deg2rad(116 - ss*0.01)
    arc1 = NP.zeros(2,dtype='int')
    arc1[0] = int(diameter/2 + arc1radius * NP.sin(arc1theta + .5))
    arc1[1] = int(diameter/2 + arc1radius * NP.cos(arc1theta + .5))
    for ii in range(16):
        cell = lcpImagesList[ss,ii,arc1[0]:arc1[0] + arcDY,arc1[1]:arc1[1] + arcDX]
        if (cellMax):
            lcpArcPoint0_t[ii,ss] = cell.max()
        else:
            lcpArcPoint0_t[ii,ss] = cell.mean()
        cell = rcpImagesList[ss,ii,arc1[0]:arc1[0] + arcDY,arc1[1]:arc1[1] + arcDX]
        if (cellMax):
            rcpArcPoint0_t[ii,ss] = cell.max()
        else:
            rcpArcPoint0_t[ii,ss] = cell.mean()

iImagesList = 0.5*(lcpImagesList + rcpImagesList)

PL.figure()
NN = 120
MM = 100
ss0 = 40
R0 = 390
C0 = 0
#R0 = 110
#C0 = 250
for ss0 in range(120):
    ss1 = ss0 + 4
    if (ss0 == 0):
        iImag0 = iImagesList[ss0:ss1,0].mean(axis=0)
        maxInd = NP.unravel_index(NP.argmax(iImag0[135:173,274:291]),(173-135,291-274))
        maxInd00 = maxInd
    pointSpectrum = NP.zeros((NN*16,MM))
    iImag = NP.zeros((16,512,512))
    for ff in range(16):
        iImag[ff] = iImagesList[ss0:ss1,ff].mean(axis=0)
        if (ff == 0):
            maxInd = NP.unravel_index(NP.argmax(iImag[ff,135:173,274:291]),(173-135,291-274))
            iImag[ff] = NP.roll(iImag[ff],NP.array(maxInd00) - NP.array(maxInd),axis=(0,1))
            maxInd = NP.unravel_index(NP.argmax(iImag[ff,135:173,274:291]),(173-135,291-274))
            maxInd0 = maxInd
        else:
            maxInd = NP.unravel_index(NP.argmax(iImag[ff,135:173,274:291]),(173-135,291-274))
            iImag[ff] = NP.roll(iImag[ff],NP.array(maxInd0) - NP.array(maxInd),axis=(0,1))
    
    #PL.figure()
    #PL.title('scan %d'%ss0)
    #PL.contour(iImag[0],levels=[5e4],cmap='rainbow')
    #for ff in range(16):
    #    PL.contourf(iImag[ff],levels=[5e3,5e4,7e4,1e5],colors=[(ff/16,0,0)])
        
    for rr in range(NN):
        for cc in range(MM):
            for ff in range(16):
                pointSpectrum[rr*16 + ff,cc] = iImag[ff,R0+rr,C0+cc]
    
    #PL.figure()
    #PL.title('scan %d'%ss0)
    #PL.imshow(pointSpectrum,vmax=5e4,origin='lower',aspect=.1)
    
    rgbPointSpectrum = NP.zeros((NN,MM,3))
    for rr in range(NN):
        for cc in range(MM):
            rgbPointSpectrum[rr,cc,0] = pointSpectrum[rr*16:rr*16+5,cc].mean()
            rgbPointSpectrum[rr,cc,1] = pointSpectrum[rr*16+5:rr*16+10,cc].mean()
            rgbPointSpectrum[rr,cc,2] = pointSpectrum[rr*16+10:rr*16+15,cc].mean()
    rgbPointSpectrum /= 1e4
    rgbPointSpectrum = rgbPointSpectrum.clip(0,1)
    
#    PL.figure()
    PL.clf()
    PL.title('scan %03d'%ss0)
    PL.imshow(rgbPointSpectrum,origin='lower')
    PL.savefig('scan_%03d.png'%ss0)
    print(ss0)
