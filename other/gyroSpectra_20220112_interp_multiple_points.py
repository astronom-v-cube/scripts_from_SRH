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

lcpMaxT_SE_t = []
rcpMaxT_SE_t = []
srcSE = [184,134]

lcpMaxT_NE_t = []
rcpMaxT_NE_t = []
srcNE = [155,380]

lcpMaxT_NE1_t = []
rcpMaxT_NE1_t = []
srcNE1 = [140,340]

lcpMaxT_QS0_t = []
rcpMaxT_QS0_t = []
srcQS0 = [160,160]

lcpMaxT_QS1_t = []
rcpMaxT_QS1_t = []
srcQS1 = [370,200]

lcpMaxT_PR_t = []
rcpMaxT_PR_t = []
srcPR = [370,60]

lcpMaxT_PR1_t = NP.zeros((10,lcpImagesList.shape[0], 16))
rcpMaxT_PR1_t = NP.zeros((10,lcpImagesList.shape[0], 16))
lcpMaxT_PR2_t = NP.zeros((10,lcpImagesList.shape[0], 16))
rcpMaxT_PR2_t = NP.zeros((10,lcpImagesList.shape[0], 16))
srcPR1 = NP.zeros((10,2),dtype='int')
for i in range(10):
    srcPR1[i,0] = 380 + 12*i
    srcPR1[i,1] = 80  - 6*i

S = 10
SDX = 5
SDY = 5
lcpStrip0_t = NP.zeros((S,lcpImagesList.shape[0], 16))
rcpStrip0_t = NP.zeros((S,lcpImagesList.shape[0], 16))
strip0_d = NP.zeros((S,2),dtype='int')
for i in range(S):
#    strip0_d[i,0] = 380 - 5 + SDY*i
#    strip0_d[i,1] = 60 + 5  + SDX*i
    strip0_d[i,0] = 350 + SDY*i
    strip0_d[i,1] = 40  + SDX*i

lcpStrip1_t = NP.zeros((S,lcpImagesList.shape[0], 16))
rcpStrip1_t = NP.zeros((S,lcpImagesList.shape[0], 16))
strip1_d = NP.zeros((S,2),dtype='int')
for i in range(S):
    strip1_d[i,0] = 380 + 5 + SDY*i
    strip1_d[i,1] = 60 - 5  + SDX*i

dX = 20
dY = 30
dXPR = 60
dYPR = 60
PRdY = 10
PRdX = 5

for ss in range(lcpImagesList.shape[0]):
    lcpMaxT_SE = []
    rcpMaxT_SE = []
    lcpMaxT_NE = []
    rcpMaxT_NE = []
    lcpMaxT_NE1 = []
    rcpMaxT_NE1 = []
    lcpMaxT_QS0 = []
    rcpMaxT_QS0 = []
    lcpMaxT_QS1 = []
    rcpMaxT_QS1 = []
    lcpMaxT_PR = []
    rcpMaxT_PR = []
    lcpMaxT_PR1 = []
    rcpMaxT_PR1 = []
    for ii in range(16):
        lcpMaxT_QS0.append(lcpImagesList[ss,ii,srcQS0[0]:srcQS0[0]+dY,srcQS0[1]:srcQS0[1]+dX].max())
        rcpMaxT_QS0.append(rcpImagesList[ss,ii,srcQS0[0]:srcQS0[0]+dY,srcQS0[1]:srcQS0[1]+dX].max())
        lcpMaxT_QS1.append(lcpImagesList[ss,ii,srcQS1[0]:srcQS1[0]+dY,srcQS1[1]:srcQS1[1]+dX].max())
        rcpMaxT_QS1.append(rcpImagesList[ss,ii,srcQS1[0]:srcQS1[0]+dY,srcQS1[1]:srcQS1[1]+dX].max())

        cell = lcpImagesList[ss,ii,srcSE[0]:srcSE[0]+dY,srcSE[1]:srcSE[1]+dX] - lcpImagesList[ss,ii,srcQS0[0]:srcQS0[0]+dY,srcQS0[1]:srcQS0[1]+dX]
        lcpMaxT_SE.append(resize(cell,(2*dY,2*dX)).max())

        cell = rcpImagesList[ss,ii,srcSE[0]:srcSE[0]+dY,srcSE[1]:srcSE[1]+dX] - rcpImagesList[ss,ii,srcQS0[0]:srcQS0[0]+dY,srcQS0[1]:srcQS0[1]+dX]
        rcpMaxT_SE.append(resize(cell,(2*dY,2*dX)).max())

        cell = lcpImagesList[ss,ii,srcNE[0]:srcNE[0]+dY,srcNE[1]:srcNE[1]+dX] - lcpImagesList[ss,ii,srcQS1[0]:srcQS1[0]+dY,srcQS1[1]:srcQS1[1]+dX]
        lcpMaxT_NE.append(resize(cell,(2*dY,2*dX)).max())

        cell = rcpImagesList[ss,ii,srcNE[0]:srcNE[0]+dY,srcNE[1]:srcNE[1]+dX] - rcpImagesList[ss,ii,srcQS1[0]:srcQS1[0]+dY,srcQS1[1]:srcQS1[1]+dX]
        rcpMaxT_NE.append(resize(cell,(2*dY,2*dX)).max())

        cell = lcpImagesList[ss,ii,srcNE1[0]:srcNE1[0]+dY,srcNE1[1]:srcNE1[1]+dX] - lcpImagesList[ss,ii,srcQS1[0]:srcQS1[0]+dY,srcQS1[1]:srcQS1[1]+dX]
        lcpMaxT_NE1.append(resize(cell,(2*dY,2*dX)).max())

        cell = rcpImagesList[ss,ii,srcNE1[0]:srcNE1[0]+dY,srcNE1[1]:srcNE1[1]+dX] - rcpImagesList[ss,ii,srcQS1[0]:srcQS1[0]+dY,srcQS1[1]:srcQS1[1]+dX]
        rcpMaxT_NE1.append(resize(cell,(2*dY,2*dX)).max())

        cell = lcpImagesList[ss,ii,srcPR[0]:srcPR[0]+dYPR,srcPR[1]:srcPR[1]+dXPR] - rcpImagesList[ss,ii,srcQS1[0]:srcQS1[0]+dYPR,srcQS1[1]:srcQS1[1]+dXPR]
        lcpMaxT_PR.append(resize(cell,(2*dY,2*dX)).max())

        cell = rcpImagesList[ss,ii,srcPR[0]:srcPR[0]+dYPR,srcPR[1]:srcPR[1]+dXPR] - rcpImagesList[ss,ii,srcQS1[0]:srcQS1[0]+dYPR,srcQS1[1]:srcQS1[1]+dXPR]
        rcpMaxT_PR.append(resize(cell,(2*dY,2*dX)).max())

        for p in range(S):
            cell = lcpImagesList[ss,ii,srcPR1[p,0]:srcPR1[p,0]+PRdY,srcPR1[p,1]:srcPR1[p,1]+PRdX]
            lcpMaxT_PR1_t[p,ss,ii] = cell.max()
#            lcpMaxT_PR2_t[p,ss,ii] = cell[PRdY//2,PRdX//2]
            lcpMaxT_PR2_t[p,ss,ii] = cell.mean()
            cell = rcpImagesList[ss,ii,srcPR1[p,0]:srcPR1[p,0]+PRdY,srcPR1[p,1]:srcPR1[p,1]+PRdX]
            rcpMaxT_PR1_t[p,ss,ii] = cell.max()
#            rcpMaxT_PR2_t[p,ss,ii] = cell[PRdY//2,PRdX//2]
            rcpMaxT_PR2_t[p,ss,ii] = cell.mean()

            cell = lcpImagesList[ss,ii,strip0_d[p,0]:strip0_d[p,0]+SDY,strip0_d[p,1]:strip0_d[p,1]+SDX]
            lcpStrip0_t[p,ss,ii] = cell.max()
#            lcpStrip0_t[p,ss,ii] = cell.mean()
            cell = rcpImagesList[ss,ii,strip0_d[p,0]:strip0_d[p,0]+SDY,strip0_d[p,1]:strip0_d[p,1]+SDX]
            rcpStrip0_t[p,ss,ii] = cell.max()
#            rcpStrip0_t[p,ss,ii] = cell.mean()
            cell = lcpImagesList[ss,ii,strip1_d[p,0]:strip1_d[p,0]+SDY,strip1_d[p,1]:strip1_d[p,1]+SDX]
            lcpStrip1_t[p,ss,ii] = cell.max()
#            lcpStrip1_t[p,ss,ii] = cell.mean()
            cell = rcpImagesList[ss,ii,strip1_d[p,0]:strip1_d[p,0]+SDY,strip1_d[p,1]:strip1_d[p,1]+SDX]
            rcpStrip1_t[p,ss,ii] = cell.max()
#            rcpStrip1_t[p,ss,ii] = cell.mean()


    lcpMaxT_SE_t.append(lcpMaxT_SE)
    rcpMaxT_SE_t.append(rcpMaxT_SE)
    lcpMaxT_NE_t.append(lcpMaxT_NE)
    rcpMaxT_NE_t.append(rcpMaxT_NE)
    lcpMaxT_NE1_t.append(lcpMaxT_NE1)
    rcpMaxT_NE1_t.append(rcpMaxT_NE1)
    lcpMaxT_QS0_t.append(lcpMaxT_QS0)
    rcpMaxT_QS0_t.append(rcpMaxT_QS0)
    lcpMaxT_QS1_t.append(lcpMaxT_QS1)
    rcpMaxT_QS1_t.append(rcpMaxT_QS1)
    lcpMaxT_PR_t.append(lcpMaxT_PR)
    rcpMaxT_PR_t.append(rcpMaxT_PR)

maxScan = 127
lcpMaxT_SE_t = NP.array(lcpMaxT_SE_t[0:maxScan])
rcpMaxT_SE_t = NP.array(rcpMaxT_SE_t[0:maxScan])
lcpMaxT_NE_t = NP.array(lcpMaxT_NE_t[0:maxScan])
rcpMaxT_NE_t = NP.array(rcpMaxT_NE_t[0:maxScan])
lcpMaxT_NE1_t = NP.array(lcpMaxT_NE1_t[0:maxScan])
rcpMaxT_NE1_t = NP.array(rcpMaxT_NE1_t[0:maxScan])
lcpMaxT_QS0_t = NP.array(lcpMaxT_QS0_t[0:maxScan])
rcpMaxT_QS0_t = NP.array(rcpMaxT_QS0_t[0:maxScan])
lcpMaxT_QS1_t = NP.array(lcpMaxT_QS1_t[0:maxScan])
rcpMaxT_QS1_t = NP.array(rcpMaxT_QS1_t[0:maxScan])
lcpMaxT_PR_t = NP.array(lcpMaxT_PR_t[0:maxScan])
rcpMaxT_PR_t = NP.array(rcpMaxT_PR_t[0:maxScan])

#------------------------------------------------------------------------------------------------------------
fig = PL.figure(figsize=(12,6))
pl = fig.subplots(nrows=2,ncols=2)
for rr in range(2):
    for cc in range(2):
        pl[rr,cc].xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
        pl[rr,cc].xaxis.set_major_locator(PL.MultipleLocator(18))
        pl[rr,cc].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
        pl[rr,cc].yaxis.set_major_locator(PL.MultipleLocator(3))
        pl[rr,cc].set_xlabel('UT')
        pl[rr,cc].set_ylabel('GHz')

Tsat = 1.5e5
Tmin = 3e4
pl[0,0].set_title('x-mode (LCP SE)')
pl[1,0].set_title('o-mode (RCP SE)')
pos0 = pl[0,0].imshow(lcpMaxT_SE_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
pos1 = pl[1,0].imshow(rcpMaxT_SE_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')

Tsat = 4.e5
Tmin = 1e5
pl[0,1].set_title('o-mode (LCP NE1)')
pl[1,1].set_title('x-mode (RCP NE1)')
pos0 = pl[0,1].imshow(lcpMaxT_NE1_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
pos1 = pl[1,1].imshow(rcpMaxT_NE1_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
#fig.colorbar(pos0,ax=pl[0,1])
#fig.colorbar(pos1,ax=pl[1,1])

#------------------------------------------------------------------------------------------------------------
fig = PL.figure(figsize=(12,6))
pl = fig.subplots(nrows=2,ncols=2)
for rr in range(2):
    for cc in range(2):
        pl[rr,cc].xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
        pl[rr,cc].xaxis.set_major_locator(PL.MultipleLocator(18))
        pl[rr,cc].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
        pl[rr,cc].yaxis.set_major_locator(PL.MultipleLocator(3))
        pl[rr,cc].set_xlabel('UT')
        pl[rr,cc].set_ylabel('GHz')

Tsat = 3.5e5
Tmin = 1e5
cmap = 'hot'
pl[0,0].set_title('o-mode (LCP) NE1')
pos0 = pl[0,0].imshow(lcpMaxT_NE1_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap=cmap)

pl[1,0].set_title('x-mode (RCP) NE1')
pos1 = pl[1,0].imshow(rcpMaxT_NE1_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap=cmap)

Tsat = 3.e5
Tmin = .7e5
pl[0,1].set_title('x-mode (LCP) NE')
pos0 = pl[0,1].imshow(lcpMaxT_NE_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap=cmap)

pl[1,1].set_title('o-mode (RCP) NE')
pos1 = pl[1,1].imshow(rcpMaxT_NE_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap=cmap)

#------------------------------------------------------------------------------------------------------------

Tsat = 1e5
Tmin = 0e3
fig = PL.figure(figsize=(12,6))
pl = fig.subplots(nrows=2,ncols=2)
for rr in range(2):
    for cc in range(2):
        pl[rr,cc].xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
        pl[rr,cc].xaxis.set_major_locator(PL.MultipleLocator(18))
        pl[rr,cc].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
        pl[rr,cc].yaxis.set_major_locator(PL.MultipleLocator(3))
        pl[rr,cc].set_xlabel('UT')
        pl[rr,cc].set_ylabel('GHz')

pl[0,0].set_title('x-mode')
pl[1,0].set_title('o-mode')
pos0 = pl[0,0].imshow(lcpMaxT_QS0_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
pos1 = pl[1,0].imshow(rcpMaxT_QS0_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
pos0 = pl[0,1].imshow(lcpMaxT_QS1_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
pos1 = pl[1,1].imshow(rcpMaxT_QS1_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
fig.colorbar(pos0,ax=pl[0])
fig.colorbar(pos1,ax=pl[1])

#------------------------------------------------------------------------------------------------------------
Tsat = 1e5
Tmin = 0e3
fig = PL.figure(figsize=(12,6))
pl = fig.subplots(nrows=2,ncols=2)
for rr in range(2):
    for cc in range(2):
        pl[rr,cc].xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
        pl[rr,cc].xaxis.set_major_locator(PL.MultipleLocator(18))
        pl[rr,cc].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
        pl[rr,cc].yaxis.set_major_locator(PL.MultipleLocator(3))
        pl[rr,cc].set_xlabel('UT')
        pl[rr,cc].set_ylabel('GHz')

pl[0,0].set_title('LCP SE (x-mode)')
pl[1,0].set_title('RCP SE (o-mode)')
pos0 = pl[0,0].imshow(lcpMaxT_SE_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
pos1 = pl[1,0].imshow(rcpMaxT_SE_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
pl[0,1].set_title('LCP PR')
pl[1,1].set_title('RCP PR')
pos0 = pl[0,1].imshow(lcpMaxT_PR_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
pos1 = pl[1,1].imshow(rcpMaxT_PR_t.T,aspect=3., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')

#------------------------------------------------------------------------------------------------------------
Tsat = 1e5
Tmin = 0e3
fig = PL.figure(figsize=(12,6))
pl = fig.subplots(nrows=10,ncols=2)
for rr in range(10):
    for cc in range(2):
        pl[rr,cc].xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
        pl[rr,cc].xaxis.set_major_locator(PL.MultipleLocator(18))
        pl[rr,cc].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
        pl[rr,cc].yaxis.set_major_locator(PL.MultipleLocator(3))
        pl[rr,cc].set_xlabel('UT')
        pl[rr,cc].set_ylabel('GHz')
#    pl[rr,0].imshow(lcpMaxT_PR1_t[rr].T,aspect=1., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
#    pl[rr,1].imshow(rcpMaxT_PR1_t[rr].T,aspect=1., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
    pl[rr,0].imshow(rcpMaxT_PR1_t[rr].T + lcpMaxT_PR1_t[rr].T,aspect=1., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')
    pl[rr,1].imshow(rcpMaxT_PR1_t[rr].T - lcpMaxT_PR1_t[rr].T,aspect=1., interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap='rainbow')

#------------------------------------------------------------------------------------------------------------
Tsat = 2e4
Tmin = 0e3
cmap = 'gray_r'
fig = PL.figure(figsize=(10,8))
pl = fig.subplots(nrows=10,ncols=1)
for rr in range(10):
    pl[rr].xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
    pl[rr].xaxis.set_major_locator(PL.MultipleLocator(18))
    pl[rr].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
    pl[rr].yaxis.set_major_locator(PL.MultipleLocator(6))
    pl[rr].set_xlabel('UT')
    pl[rr].set_ylabel('GHz')
    if (rr > 0 and rr < 9):
        pl[rr].axis('off')
    pl[rr].imshow(rcpMaxT_PR2_t[rr,0:90].T + lcpMaxT_PR2_t[rr,0:90].T,aspect=.5, interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap=cmap)


#------------------------------------------------------------------------------------------------------------
Tsat = 5e4
Tmin = 0e3
cmap = 'gray_r'
fig = PL.figure(figsize=(10,8))
fig.subplots_adjust(hspace=0.001,wspace=0.001,top=0.95,bottom=0.01,left=0.01,right=0.99)
pl = fig.subplots(nrows=10,ncols=1)
for rr in range(10):
    pl[rr].xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
    pl[rr].xaxis.set_major_locator(PL.MultipleLocator(18))
    pl[rr].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
    pl[rr].yaxis.set_major_locator(PL.MultipleLocator(6))
    pl[rr].set_xlabel('UT')
    pl[rr].set_ylabel('GHz')
    pl[rr].axis('off')
#    pl[rr].imshow(rcpStrip0_t[rr,0:90].T + lcpStrip0_t[rr,0:90].T,aspect=.5, interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap=cmap)
    pl[rr].imshow(rcpStrip1_t[rr,0:90].T + lcpStrip1_t[rr,0:90].T,aspect=.5, interpolation='bessel', vmin=Tmin, vmax=Tsat, cmap=cmap)

