#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 11:49:23 2021

@author: sergey_lesovoi
"""
#20211111 scan0 04:13:39 scans 300

import numpy as NP
import pylab as PL
from skimage import measure#, resize
from astropy.io import fits
from scipy import interpolate
import scipy.optimize as opt

#def fitTb(t,A,B,C):
#    return A + B/t**C

def fitTb(t,A,B,C):
    return A + (t-B)**C

def fitTb1(t,A,B,C,D):
    return A + B*(t-C)**D

def fitTb2(t,A,B,C):
    return A*(t+B)**C

def fitTb3(t,A,B):
    return A*(t+220)**B
    
def fitTb4(l,A,B):
    return A*(l+.4)**B
    
def freqGHzFormat(f, pos):
    return '%3.2f'%(f * 0.4 + 5.8)

def times1_format(t, pos):
  t = int(t)
  if (t < cpTime.shape[0]):
      if (t < 10):
          t = cpTime[t + scans1[0],0]
          hh = int(t / 3600.)
          t -= hh*3600.
          mm = int(t / 60.)
          t -= mm*60.
          ss = int(t)
          return '%02d:%02d:%02d' % (hh,mm,ss)
      else:
          t *= 10
          mm = int(t / 60.)
          t -= mm*60.
          ss = int(t)
          return '%02d:%02d' % (mm,ss)
  else:
      t *= 10
      t += cpTime[scans1[0],0]
      hh = int(t / 3600.)
      t -= hh*3600.
      mm = int(t / 60.)
      t -= mm*60.
      ss = int(t)
      return '%02d:%02d:%02d' % (hh,mm,ss)

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

def formatCpTime(t):
    t = cpTime[int(t),0]
    hh = int(t / 3600.)
    t -= hh*3600.
    mm = int(t / 60.)
    t -= mm*60.
    ss = int(t)
    return '%02d:%02d:%02d' % (hh,mm,ss)
    
    
def findMaxPath(iArray, frequency, scanFirst,scanLast,fromPoint,findWindow):
    maxPath = []
    currentWindow = NP.zeros((2,2),dtype='int')
    currentWindow[0,0] = fromPoint[0] - findWindow[0]
    currentWindow[0,1] = fromPoint[0] + findWindow[0]
    currentWindow[1,0] = fromPoint[1] - findWindow[1]
    currentWindow[1,1] = fromPoint[1] + findWindow[1]
    windowShape = [findWindow[0]*2, findWindow[1]*2]
    
    for point in range(scanLast - scanFirst):
        maxInd = NP.array(NP.unravel_index(NP.argmax(iArray[point + scanFirst,frequency,currentWindow[0,0]:currentWindow[0,1],currentWindow[1,0]:currentWindow[1,1]]),windowShape))
        maxInd[0] += currentWindow[0,0]
        maxInd[1] += currentWindow[1,0]
        maxPath.append(maxInd)
        currentWindow[0,0] = maxInd[0] - findWindow[0]
        currentWindow[0,1] = maxInd[0] + findWindow[0]
        currentWindow[1,0] = maxInd[1] - findWindow[1]
        currentWindow[1,1] = maxInd[1] + findWindow[1]
    
    return NP.array(maxPath)


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

#------------------------------------------------------------------------------------------------------------
for ss in range(lcpImagesList.shape[0]):
    for ff in range(16):
        lcpMaxInd = NP.unravel_index(NP.argmax(lcpImagesList[ss,ff,135:173,274:291]),(173-135,291-274))
        rcpMaxInd = NP.unravel_index(NP.argmax(rcpImagesList[ss,ff,135:173,274:291]),(173-135,291-274))
        rcpImagesList[ss,ff] = NP.roll(rcpImagesList[ss,ff],NP.array(lcpMaxInd) - NP.array(rcpMaxInd),axis=(0,1))
    
iImagesList = 0.5*(lcpImagesList + rcpImagesList)

PL.figure()
NN = 130
MM = 90
ss0 = 40
R0 = 360
C0 = 0
#R0 = 165
#0 = 110
dynamicPointSpectrum = NP.zeros_like(iImagesList)
for ss0 in range(iImagesList.shape[0]):
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
            dynamicPointSpectrum[ss0,ff,:,:] = iImag[ff]
            maxInd = NP.unravel_index(NP.argmax(iImag[ff,135:173,274:291]),(173-135,291-274))
            maxInd0 = maxInd
        else:
            maxInd = NP.unravel_index(NP.argmax(iImag[ff,135:173,274:291]),(173-135,291-274))
            iImag[ff] = NP.roll(iImag[ff],NP.array(maxInd0) - NP.array(maxInd),axis=(0,1))
            dynamicPointSpectrum[ss0,ff,:,:] = iImag[ff]
    
    for rr in range(NN):
        for cc in range(MM):
            for ff in range(16):
                pointSpectrum[rr*16 + ff,cc] = iImag[ff,R0+rr,C0+cc]
    
pointNumber = 10
#point00 = NP.array([380,90])
point00 = NP.array([385,85])
point01 = NP.array([460,40])
point02 = NP.array([460,50])
point03 = NP.array([470,60])
point04 = NP.array([470,80])

dP01 = (point01 - point00)/(pointNumber - 1)
dP02 = (point02 - point00)/(pointNumber - 1)
dP03 = (point03 - point00)/(pointNumber - 1)
dP04 = (point04 - point00)/(pointNumber - 1)

spectrumPoints = NP.zeros((pointNumber,4,2))
for p in range(pointNumber):
    spectrumPoints[p,0] = point00 + [p*dP01[0], p*dP01[1]]
    spectrumPoints[p,1] = point00 + [p*dP02[0], p*dP02[1]]
    spectrumPoints[p,2] = point00 + [p*dP03[0], p*dP03[1]]
    spectrumPoints[p,3] = point00 + [p*dP03[0], p*dP04[1]]

fig = PL.figure(figsize=(12,10))
pl = fig.subplots(nrows=pointNumber,ncols=4)
for p in range(pointNumber):
    for l in range(4):
        pl[p,l].imshow(dynamicPointSpectrum[:,:,int(spectrumPoints[p,l,0]+.5),int(spectrumPoints[p,l,1]+.5)].T,interpolation='bessel',aspect=2,vmax=1e4,vmin=1e3)
        
PL.figure()
PL.imshow(dynamicPointSpectrum.mean(axis=(0,1)),origin='lower',vmin=1e2,vmax=5e3)
PL.plot(spectrumPoints[:,0,1],spectrumPoints[:,0,0],'.')
PL.plot(spectrumPoints[:,1,1],spectrumPoints[:,1,0],'.')
PL.plot(spectrumPoints[:,2,1],spectrumPoints[:,2,0],'.')
PL.plot(spectrumPoints[:,3,1],spectrumPoints[:,3,0],'.')

piece1 = NP.array([\
                  [32,403,80],\
                  [36,412,73],\
                  [40,421,66],\
                  [44,426,60],\
                  [48,437,53],\
                  [52,443,46],\
                  [56,450,40],\
                  [60,457,35],\
                  [64,464,27],\
                  [68,469,24],\
                  [72,481,20]\
        ])
scans1 = NP.linspace(piece1.T[0,0],piece1.T[0,-1],piece1.T[0,-1]-piece1.T[0,0]+1,dtype='int')
rows1 = NP.array(NP.interp(scans1, piece1.T[0], piece1.T[1] + 4)+.5,dtype='int')
cols1 = NP.array(NP.interp(scans1, piece1.T[0], piece1.T[2])+.5,dtype='int')
dSpectrum1 = NP.zeros((70,16))
dSpectrum1[32-32:73-32,:] = dynamicPointSpectrum[scans1,:,rows1,cols1]

piece2 = NP.array([\
                  [36,401,78],\
                  [40,417,76],\
                  [44,424,71],\
                  [48,432,68],\
                  [52,441,60],\
                  [56,448,53],\
                  [60,460,48],\
                  [64,468,42],\
                  [68,482,38],\
                  [72,492,31],\
                  [76,507,26]\
        ])

scans2 = NP.linspace(piece2.T[0,0],piece2.T[0,-1],piece2.T[0,-1]-piece2.T[0,0]+1,dtype='int')
rows2 = NP.array(NP.interp(scans2, piece2.T[0], piece2.T[1] + 4)+.5,dtype='int')
cols2 = NP.array(NP.interp(scans2, piece2.T[0], piece2.T[2])+.5,dtype='int')
dSpectrum2 = NP.zeros((70,16))
dSpectrum2[32-32:73-32,:] = dynamicPointSpectrum[scans2,:,rows2,cols2]

piece3 = NP.array([\
                  [32, 393,83],\
                  [42, 393,83],\
                  [45, 387,83],\
                  [50, 387,81],\
                  [60, 399,72],\
                  [70, 411,61],\
                  [80, 419,53],\
                  [90, 431,43],\
                  [100,444,34]\
        ])

scans3 = NP.linspace(piece3.T[0,0],piece3.T[0,-1],piece3.T[0,-1]-piece3.T[0,0]+1,dtype='int')
rows3 = NP.array(NP.interp(scans3, piece3.T[0], piece3.T[1] + 4)+.5,dtype='int')
cols3 = NP.array(NP.interp(scans3, piece3.T[0], piece3.T[2])+.5,dtype='int')
dSpectrum3 = NP.zeros((70,16))
dSpectrum3[32-32:101-32,:] = dynamicPointSpectrum[scans3,:,rows3,cols3]

piece4 = NP.array([\
                  [32, 395,87],\
                  [42, 394,87],\
                  [45, 395,81],\
                  [50, 403,77],\
                  [60, 417,69],\
                  [70, 428,60],\
                  [80, 443,51],\
                  [90, 457,44],\
                  [100,471,37]\
        ])

scans4 = NP.linspace(piece4.T[0,0],piece4.T[0,-1],piece4.T[0,-1]-piece3.T[0,0]+1,dtype='int')
rows4 = NP.array(NP.interp(scans4, piece4.T[0], piece4.T[1] + 4)+.5,dtype='int')
cols4 = NP.array(NP.interp(scans4, piece4.T[0], piece4.T[2])+.5,dtype='int')
dSpectrum4 = NP.zeros((70,16))
dSpectrum4[32-32:101-32,:] = dynamicPointSpectrum[scans4,:,rows4,cols4]

seconds1 = scans1*3*3.2
seconds2 = scans2*3*3.2
seconds3 = scans3*3*3.2
seconds4 = scans4*3*3.2
RR0 = 345
CC0 = 93
distance1 = NP.sqrt((rows1-RR0)**2 + (cols1-CC0)**2)*4.911*700*1e-6
distance2 = NP.sqrt((rows2-RR0)**2 + (cols2-CC0)**2)*4.911*700*1e-6
distance3 = NP.sqrt((rows3-RR0)**2 + (cols3-CC0)**2)*4.911*700*1e-6
distance4 = NP.sqrt((rows4-RR0)**2 + (cols4-CC0)**2)*4.911*700*1e-6

iSeconds = NP.linspace(1,1000,200)
iDistance = NP.linspace(.1e6,.7e6,200)*1e-6
#------------------------------------------------------------------------------
cmap='jet'
times1 = cpTime[scans1,0]
times2 = cpTime[scans2,0]
times3 = cpTime[scans3,0]
times4 = cpTime[scans4,0]
Tmax1 = 3e4
Tmax2 = 3e4
Tmax3 = 3e4
Tmax4 = 3e4
fig = PL.figure()
pl = fig.subplots(nrows=4,ncols=1)
pl[0].axes.get_xaxis().set_visible(False)
pl[1].axes.get_xaxis().set_visible(False)
pl[2].axes.get_xaxis().set_visible(False)
pl[3].xaxis.set_major_formatter(PL.FuncFormatter(times1_format))
pl[3].xaxis.set_major_locator(PL.MultipleLocator(18))

pl[0].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[0].yaxis.set_major_locator(PL.MultipleLocator(6))
pl[0].set_ylabel('GHz')
pl[1].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[1].yaxis.set_major_locator(PL.MultipleLocator(6))
pl[1].set_ylabel('GHz')
pl[2].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[2].yaxis.set_major_locator(PL.MultipleLocator(6))
pl[2].set_ylabel('GHz')
pl[3].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[3].yaxis.set_major_locator(PL.MultipleLocator(6))
pl[3].set_xlabel('UT')
pl[3].set_ylabel('GHz')

pos1 = pl[0].imshow(dSpectrum1.T,vmin=0,vmax=Tmax1,interpolation='bessel',cmap=cmap)
pos2 = pl[1].imshow(dSpectrum2.T,vmin=0,vmax=Tmax2,interpolation='bessel',cmap=cmap)
pos3 = pl[2].imshow(dSpectrum3.T,vmin=0,vmax=Tmax3,interpolation='bessel',cmap=cmap)
pos4 = pl[3].imshow(dSpectrum4.T,vmin=0,vmax=Tmax4,interpolation='bessel',cmap=cmap)
fig.colorbar(pos1,ax=pl[0])
fig.colorbar(pos2,ax=pl[1])
fig.colorbar(pos3,ax=pl[2])
fig.colorbar(pos4,ax=pl[3])

#------------------------------------------------------------------------------
fitPar1 = opt.curve_fit(fitTb3, seconds1,dSpectrum1[0:41,0:3].mean(axis=1),p0=[4e8,-4.])
fitPar2 = opt.curve_fit(fitTb3, seconds2,dSpectrum2[0:41,0:3].mean(axis=1),p0=[4e8,-4.])
fitPar3 = opt.curve_fit(fitTb3, seconds3,dSpectrum3[0:69,0:3].mean(axis=1),p0=[4e8,-4.])
fitPar4 = opt.curve_fit(fitTb3, seconds4,dSpectrum4[0:69,0:3].mean(axis=1),p0=[4e8,-4.])

PL.figure(figsize=(8,4))
PL.plot(iSeconds, fitTb3(iSeconds,*fitPar1[0]),label='%.2f'%fitPar1[0][1], color='red')
PL.plot(iSeconds, fitTb3(iSeconds,*fitPar2[0]),label='%.2f'%fitPar2[0][1], color='green')
PL.plot(iSeconds, fitTb3(iSeconds,*fitPar3[0]),label='%.2f'%fitPar3[0][1], color='blue')
PL.plot(iSeconds, fitTb3(iSeconds,*fitPar4[0]),label='%.2f'%fitPar4[0][1], color='orange')

PL.plot(seconds1, dSpectrum1[0:41,0:3].mean(axis=1),'.', markersize=3, color='red')
PL.plot(seconds2, dSpectrum2[0:41,0:3].mean(axis=1),'.', markersize=3, color='green')
PL.plot(seconds3, dSpectrum3[0:69,0:3].mean(axis=1),'.', markersize=3, color='blue')
PL.plot(seconds4, dSpectrum4[0:69,0:3].mean(axis=1),'.', markersize=3, color='orange')

PL.legend()
PL.xlim(0,1e3)
PL.ylim(0,1.5e4)
PL.xlabel('seconds')
PL.ylabel('Tb')
PL.grid()
#------------------------------------------------------------------------------
fitPar1 = opt.curve_fit(fitTb4, distance1,dSpectrum1[0:41,0:3].mean(axis=1),p0=[1e5,-5.])
fitPar2 = opt.curve_fit(fitTb4, distance2,dSpectrum2[0:41,0:3].mean(axis=1),p0=[1e5,-5.])
fitPar3 = opt.curve_fit(fitTb4, distance3,dSpectrum3[0:69,0:3].mean(axis=1),p0=[1e5,-5.])
fitPar4 = opt.curve_fit(fitTb4, distance4,dSpectrum4[0:69,0:3].mean(axis=1),p0=[1e5,-5.])

PL.figure(figsize=(8,4))
PL.plot(iDistance, fitTb4(iDistance,*fitPar1[0]),label='%.2f'%fitPar1[0][1], color='red')
PL.plot(iDistance, fitTb4(iDistance,*fitPar2[0]),label='%.2f'%fitPar2[0][1], color='green')
PL.plot(iDistance, fitTb4(iDistance,*fitPar3[0]),label='%.2f'%fitPar3[0][1], color='blue')
PL.plot(iDistance, fitTb4(iDistance,*fitPar4[0]),label='%.2f'%fitPar4[0][1], color='orange')

PL.plot(distance1, dSpectrum1[0:41,0:3].mean(axis=1),'.', markersize=3, color='red')
PL.plot(distance2, dSpectrum2[0:41,0:3].mean(axis=1),'.', markersize=3, color='green')
PL.plot(distance3, dSpectrum3[0:69,0:3].mean(axis=1),'.', markersize=3, color='blue')
PL.plot(distance4, dSpectrum4[0:69,0:3].mean(axis=1),'.', markersize=3, color='orange')

PL.legend()
#PL.xlim(0,1e3)
#PL.ylim(0,1.5e4)
PL.xlabel('km')
PL.ylabel('Tb')
PL.grid()
#------------------------------------------------------------------------------

PL.figure()
PL.imshow(dynamicPointSpectrum.mean(axis=(0,1)),origin='lower',vmin=0,vmax=1e4)
PL.plot(piece1[:,2], piece1[:,1], '.')
PL.plot(piece2[:,2], piece2[:,1], '.')
PL.plot(piece3[:,2], piece3[:,1], '.')
PL.plot(piece4[:,2], piece4[:,1], '.')

#---------------------------------------------------------------------------
srcRC = [188,135]
srcDS = 10
lcpSrcDynSpectrum = NP.zeros((lcpImagesList.shape[0],lcpImagesList.shape[1]))
rcpSrcDynSpectrum = NP.zeros((lcpImagesList.shape[0],lcpImagesList.shape[1]))
for tt in range(lcpImagesList.shape[0]):
    for ff in range(lcpImagesList.shape[1]):
        lcpSrcDynSpectrum[tt,ff] = lcpImagesList[tt,ff,srcRC[0]:srcRC[0]+srcDS,srcRC[1]:srcRC[1]+srcDS].max()
        rcpSrcDynSpectrum[tt,ff] = rcpImagesList[tt,ff,srcRC[0]:srcRC[0]+srcDS,srcRC[1]:srcRC[1]+srcDS].max()
        
fig = PL.figure(figsize=(8,8))
pl = fig.subplots(nrows=2,ncols=1) 
pl[0].xaxis.set_major_formatter(PL.FuncFormatter(times1_format))
pl[0].xaxis.set_major_locator(PL.MultipleLocator(18))
pl[0].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[0].yaxis.set_major_locator(PL.MultipleLocator(6))
pos = pl[0].imshow(lcpSrcDynSpectrum.T,aspect=5,interpolation='bessel',vmax=1.5e5,cmap='jet')
fig.colorbar(pos,ax=pl[0])

pl[1].xaxis.set_major_formatter(PL.FuncFormatter(times1_format))
pl[1].xaxis.set_major_locator(PL.MultipleLocator(18))
pl[1].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[1].yaxis.set_major_locator(PL.MultipleLocator(6))
pos = pl[1].imshow(rcpSrcDynSpectrum.T,aspect=5,interpolation='bessel',vmax=1.5e5,cmap='jet')
fig.colorbar(pos,ax=pl[1])
#---------------------------------------------------------------------------
srcRC = [315,370]
srcDS = 10
lcpSrcDynSpectrum = NP.zeros((lcpImagesList.shape[0],lcpImagesList.shape[1]))
rcpSrcDynSpectrum = NP.zeros((lcpImagesList.shape[0],lcpImagesList.shape[1]))
for tt in range(lcpImagesList.shape[0]):
    for ff in range(lcpImagesList.shape[1]):
        lcpSrcDynSpectrum[tt,ff] = lcpImagesList[tt,ff,srcRC[0]:srcRC[0]+srcDS,srcRC[1]:srcRC[1]+srcDS].max()
        rcpSrcDynSpectrum[tt,ff] = rcpImagesList[tt,ff,srcRC[0]:srcRC[0]+srcDS,srcRC[1]:srcRC[1]+srcDS].max()
        
fig = PL.figure(figsize=(8,8))
pl = fig.subplots(nrows=2,ncols=1) 
pl[0].xaxis.set_major_formatter(PL.FuncFormatter(times1_format))
pl[0].xaxis.set_major_locator(PL.MultipleLocator(18))
pl[0].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[0].yaxis.set_major_locator(PL.MultipleLocator(6))
pos = pl[0].imshow(lcpSrcDynSpectrum.T,aspect=5,interpolation='bessel',vmax=1.5e5,cmap='jet')
fig.colorbar(pos,ax=pl[0])

pl[1].xaxis.set_major_formatter(PL.FuncFormatter(times1_format))
pl[1].xaxis.set_major_locator(PL.MultipleLocator(18))
pl[1].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[1].yaxis.set_major_locator(PL.MultipleLocator(6))
pos = pl[1].imshow(rcpSrcDynSpectrum.T,aspect=5,interpolation='bessel',vmax=1.5e5,cmap='jet')
fig.colorbar(pos,ax=pl[1])
#---------------------------------------------------------------------------
srcRC = [145,273]
srcDS = 10
lcpSrcDynSpectrum = NP.zeros((lcpImagesList.shape[0],lcpImagesList.shape[1]))
rcpSrcDynSpectrum = NP.zeros((lcpImagesList.shape[0],lcpImagesList.shape[1]))
for tt in range(lcpImagesList.shape[0]):
    for ff in range(lcpImagesList.shape[1]):
        lcpSrcDynSpectrum[tt,ff] = lcpImagesList[tt,ff,srcRC[0]:srcRC[0]+srcDS,srcRC[1]:srcRC[1]+srcDS].max()
        rcpSrcDynSpectrum[tt,ff] = rcpImagesList[tt,ff,srcRC[0]:srcRC[0]+srcDS,srcRC[1]:srcRC[1]+srcDS].max()
        
fig = PL.figure(figsize=(8,8))
pl = fig.subplots(nrows=2,ncols=1) 
pl[0].xaxis.set_major_formatter(PL.FuncFormatter(times1_format))
pl[0].xaxis.set_major_locator(PL.MultipleLocator(18))
pl[0].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[0].yaxis.set_major_locator(PL.MultipleLocator(6))
pos = pl[0].imshow(lcpSrcDynSpectrum.T,aspect=5,interpolation='bessel',vmax=1.5e5,cmap='jet')
fig.colorbar(pos,ax=pl[0])

pl[1].xaxis.set_major_formatter(PL.FuncFormatter(times1_format))
pl[1].xaxis.set_major_locator(PL.MultipleLocator(18))
pl[1].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[1].yaxis.set_major_locator(PL.MultipleLocator(6))
pos = pl[1].imshow(rcpSrcDynSpectrum.T,aspect=5,interpolation='bessel',vmax=1.5e5,cmap='jet')
fig.colorbar(pos,ax=pl[1])
#---------------------------------------------------------------------------
srcRC = [140,345]
srcDS = 15
lcpSrcDynSpectrum = NP.zeros((lcpImagesList.shape[0],lcpImagesList.shape[1]))
rcpSrcDynSpectrum = NP.zeros((lcpImagesList.shape[0],lcpImagesList.shape[1]))
for tt in range(lcpImagesList.shape[0]):
    for ff in range(lcpImagesList.shape[1]):
        lcpSrcDynSpectrum[tt,ff] = lcpImagesList[tt,ff,srcRC[0]:srcRC[0]+srcDS,srcRC[1]:srcRC[1]+srcDS].max()
        rcpSrcDynSpectrum[tt,ff] = rcpImagesList[tt,ff,srcRC[0]:srcRC[0]+srcDS,srcRC[1]:srcRC[1]+srcDS].max()
        
fig = PL.figure(figsize=(8,8))
pl = fig.subplots(nrows=2,ncols=1) 
pl[0].xaxis.set_major_formatter(PL.FuncFormatter(times1_format))
pl[0].xaxis.set_major_locator(PL.MultipleLocator(18))
pl[0].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[0].yaxis.set_major_locator(PL.MultipleLocator(6))
pos = pl[0].imshow(lcpSrcDynSpectrum.T,aspect=5,interpolation='bessel',vmax=4e5,cmap='jet')
fig.colorbar(pos,ax=pl[0])

pl[1].xaxis.set_major_formatter(PL.FuncFormatter(times1_format))
pl[1].xaxis.set_major_locator(PL.MultipleLocator(18))
pl[1].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[1].yaxis.set_major_locator(PL.MultipleLocator(6))
pos = pl[1].imshow(rcpSrcDynSpectrum.T,aspect=5,interpolation='bessel',vmax=4e5,cmap='jet')
fig.colorbar(pos,ax=pl[1])
#---------------------------------------------------------------------------
srcRC = [150,380]
srcDS = 15
lcpSrcDynSpectrum = NP.zeros((lcpImagesList.shape[0],lcpImagesList.shape[1]))
rcpSrcDynSpectrum = NP.zeros((lcpImagesList.shape[0],lcpImagesList.shape[1]))
for tt in range(lcpImagesList.shape[0]):
    for ff in range(lcpImagesList.shape[1]):
        lcpSrcDynSpectrum[tt,ff] = lcpImagesList[tt,ff,srcRC[0]:srcRC[0]+srcDS,srcRC[1]:srcRC[1]+srcDS].max()
        rcpSrcDynSpectrum[tt,ff] = rcpImagesList[tt,ff,srcRC[0]:srcRC[0]+srcDS,srcRC[1]:srcRC[1]+srcDS].max()
        
fig = PL.figure(figsize=(8,8))
pl = fig.subplots(nrows=2,ncols=1) 
pl[0].xaxis.set_major_formatter(PL.FuncFormatter(times1_format))
pl[0].xaxis.set_major_locator(PL.MultipleLocator(18))
pl[0].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[0].yaxis.set_major_locator(PL.MultipleLocator(6))
pos = pl[0].imshow(lcpSrcDynSpectrum.T,aspect=5,interpolation='bessel',vmax=3e5,cmap='jet')
fig.colorbar(pos,ax=pl[0])


pl[1].xaxis.set_major_formatter(PL.FuncFormatter(times1_format))
pl[1].xaxis.set_major_locator(PL.MultipleLocator(18))
pl[1].yaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl[1].yaxis.set_major_locator(PL.MultipleLocator(6))
pos = pl[1].imshow(rcpSrcDynSpectrum.T,aspect=5,interpolation='bessel',vmax=3e5,cmap='jet')
fig.colorbar(pos,ax=pl[1])
#---------------------------------------------------------------------------
PL.figure()
PL.imshow(iImagesList[:,:].mean(axis=(0,1)),origin='lower',vmin=0,vmax=5e4,cmap='jet')
PL.axis('off')
#---------------------------------------------------------------------------
R0 = 360
C0 = 0
NN = 150
MM = 110
for ss in range(lcpImagesList.shape[0]):
    rgbPointSpectrum = NP.zeros((NN,MM,3))
    for rr in range(NN):
        for cc in range(MM):
            rgbPointSpectrum[rr,cc,0] = dynamicPointSpectrum[ss,0:5,rr + R0,cc + C0].mean()
            rgbPointSpectrum[rr,cc,1] = dynamicPointSpectrum[ss,5:10,rr + R0,cc + C0].mean()
            rgbPointSpectrum[rr,cc,2] = dynamicPointSpectrum[ss,10:16,rr + R0,cc + C0].mean()
    rgbPointSpectrum /= 1e4
    rgbPointSpectrum = rgbPointSpectrum.clip(0,1)
    
    PL.figure(figsize=(4,6))
    PL.title(formatCpTime(ss))
    PL.tight_layout()
    PL.ylim(0,NN)
    PL.xlim(0,MM)
    PL.axis('off')
    PL.imshow(rgbPointSpectrum,origin='lower')
    if (ss >= scans1[0] and ss <= scans1[-1]):
        PL.plot(cols1[ss - scans1[0]] - C0,rows1[ss - scans1[0]] - R0 + 0,'.',markersize=15,color='red')
    if (ss >= scans2[0] and ss <= scans2[-1]):
        PL.plot(cols2[ss - scans2[0]] - C0,rows2[ss - scans2[0]] - R0 + 0,'.',markersize=15,color='green')
    if (ss >= scans3[0] and ss <= scans3[-1]):
        PL.plot(cols3[ss - scans3[0]] - C0,rows3[ss - scans3[0]] - R0 + 0,'.',markersize=15,color='blue')
    if (ss >= scans4[0] and ss <= scans4[-1]):
        PL.plot(cols4[ss - scans4[0]] - C0,rows4[ss - scans4[0]] - R0 + 0,'.',markersize=15,color='orange')
    PL.savefig('scan_%03d.png'%ss)
    PL.close()
    print(ss)

#---------------------------------------------------------------------------
fig = PL.figure(figsize=(7.4,10))
fig.suptitle('SRH 20220112 04:30 LCP, 5.8-11.8 GHz')
pl = fig.subplots(nrows=4,ncols=4)
fig.subplots_adjust(left=0, bottom=0, right=1, top=0.95, wspace=0, hspace=0)
for rr in range(4):
    for cc in range(4):
        pl[rr,cc].imshow(lcpImagesList[0:10,rr*4 + cc].mean(axis=0),vmin=1e3,vmax=2e5,cmap='jet',origin='lower')
        pl[rr,cc].axes.get_xaxis().set_visible(False)
        pl[rr,cc].axes.get_yaxis().set_visible(False)

fig = PL.figure(figsize=(7.4,10))
fig.suptitle('SRH 20220112 04:30 RCP, 5.8-11.8 GHz')
fig.subplots_adjust(left=0, bottom=0, right=1, top=0.95, wspace=0, hspace=0)
pl = fig.subplots(nrows=4,ncols=4)
for rr in range(4):
    for cc in range(4):
        pl[rr,cc].imshow(rcpImagesList[0:10,rr*4 + cc].mean(axis=0),vmin=1e3,vmax=2e5,cmap='jet',origin='lower')
        pl[rr,cc].axes.get_xaxis().set_visible(False)
        pl[rr,cc].axes.get_yaxis().set_visible(False)
