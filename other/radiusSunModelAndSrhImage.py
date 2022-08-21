#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 01:20:56 2018

@author: sergey
"""
import os, fnmatch;
from astropy.io import fits
import pylab as PL
import numpy as NP
import scipy.signal
import scipy.special
import skimage
from skimage import filters

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename))
    result.sort()
    return result


N = 256
L = 300
pix=9.8 # arcsec per pixel
arcsecRadius = 1040 #firstFreq
degRadius = NP.deg2rad(arcsecRadius/3600)
radius = int(arcsecRadius/pix +0.5)
FOV = NP.deg2rad(pix * N / 3600)

qSun = NP.zeros((N, N))
qRing = NP.zeros((N, N))

theta_arg = FOV*NP.linspace(-.5,.5,N)
gx, gy = NP.meshgrid(theta_arg,theta_arg)
HBW = NP.deg2rad(1.)# half beam width [degree]
EBW = NP.deg2rad(.4)# half beam width [degree]
primaryBeam = NP.exp(-((gx/HBW)**2 + (gy/HBW)**2))
edgeBright = 1.2*NP.exp(-((gx/EBW)**2 + (gy/EBW)**2))

dL = 13
kern = NP.ones((dL,dL))
dLring = 21
kern_arg = NP.linspace(-1,1,dLring)
kernRingX, kernRingY = NP.meshgrid(kern_arg, kern_arg)
kernSun = NP.exp(-((kernRingX/.7)**2 + (kernRingY/.7)**2))
kernRing = NP.exp(-((kernRingX/.3)**2 + (kernRingY/.3)**2))
sunEll = skimage.measure.EllipseModel()

#srhFitsPath = '/home/svlesovoi/Documents/Python Scripts/srhFits'
srhFitsPath = '/home/svlesovoi/Documents/Python Scripts/srhFits_082002'
#srhFitsPath = '/home/svlesovoi/Documents/Python Scripts/srhFits_082003'
#srhFitsPath = '/home/svlesovoi/Documents/Python Scripts/srhFits_082004'
fitNames = findFits(srhFitsPath, '*LCP.fit')
srhImages = []
srhEdges = []
srhSunRadiusViaSobel = []
qSunModelsRadiusViaSobel = []
srhSunRadiusTheta = []
srhFreqs = []
qSunModels = []
srhImageMaxInd = []

for i in range(len(fitNames)):
    qSun[:, :] = 0.
    qRing[:, :] = 0.
    arcsecRadius = 1020 - i*2
    qRingValue = 8. - i*.5
    arcsecRadiusRing = 1010 - i*2.
    degRadius = NP.deg2rad(arcsecRadius/3600)
    degRadiusRing = NP.deg2rad(arcsecRadiusRing/3600)
    radius = int(arcsecRadius/pix +0.5)
    radiusRing = int(arcsecRadiusRing/pix +0.5)
    for i in range(N):
        for j in range(N):
            x=i - N/2
            y=j - N/2
            if (NP.sqrt(x**2 + y**2) < 1.0*radiusRing and NP.sqrt(x**2 + y**2) > 0.89*radiusRing and NP.abs(y) > NP.abs(1.*x)):
                qRing[i, j] = qRingValue
            if (NP.sqrt(x**2 + y**2) < radius):
                qSun[i, j] = 3.
    cqSun = scipy.signal.fftconvolve(qSun,kernSun) / dLring**2
    cqSun = cqSun[dLring//2:dLring//2+N,dLring//2:dLring//2+N]
    cqRing = scipy.signal.fftconvolve(qRing,kernRing) / dLring**2
    cqRing = cqRing[dLring//2:dLring//2+N,dLring//2:dLring//2+N]
    qSunModels.append((cqSun + cqRing)*primaryBeam)
    
for fitName in fitNames:
    srhFits = fits.open(fitName)
    srhFreqs.append(float(srhFits[0].header['FREQUENC']))
    srh_x_size = srhFits[0].header['NAXIS1']
    srh_y_size = srhFits[0].header['NAXIS2']
    srh_x_delt = srhFits[0].header['CDELT1'] / 60
    srh_y_delt = srhFits[0].header['CDELT2'] / 60
    srhFitsData = srhFits[0].data
    srhFitsData = srhFitsData + NP.abs(srhFitsData)
    srhImages.append(srhFitsData)
    x0 = int(srh_x_size//2)
    y0 = int(srh_y_size//2)
    dx = int(srh_x_size//4)
    dy = int(srh_y_size//4)
    srhFitsDataMean = NP.mean(srhFitsData[x0 - dx:x0 + dx, y0 - dy:y0 + dy])
    srhFitsEdge = filters.sobel(NP.clip(srhFitsData,srhFitsDataMean/1.5 - srhFitsDataMean/50,srhFitsDataMean/1.5 + srhFitsDataMean/50))
    srhEdges.append(srhFitsEdge)
    ellInd = NP.where(srhFitsEdge > .5*srhFitsEdge.max())
    ellXY = NP.zeros((len(ellInd[1]),2))
    ellXY[:,0] = ellInd[1]
    ellXY[:,1] = ellInd[0]
    ellXY[:,0] = (ellXY[:,0] - x0) * srh_x_delt
    ellXY[:,1] = (ellXY[:,1] - y0) * srh_y_delt
    sunEll.estimate(ellXY)
    xc, yc, a, b, theta = sunEll.params
    srhSunRadiusViaSobel.append(a*60)
    srhSunRadiusTheta.append(theta)
    srhImageMaxInd.append(NP.argmax(srhFitsData))


for fitNumber in range(len(fitNames)):
    qSunModelsMean = NP.mean(qSunModels[fitNumber][x0 - dx:x0 + dx, y0 - dy:y0 + dy])
    qSunModelsEdge = filters.sobel(NP.clip(qSunModels[fitNumber],qSunModelsMean/1.5 - qSunModelsMean/50,qSunModelsMean/1.5 + qSunModelsMean/50))
    ellInd = NP.where(qSunModelsEdge > NP.max(qSunModelsEdge)*.5)
    ellXY = NP.zeros((len(ellInd[1]),2))
    ellXY[:,0] = ellInd[0]
    ellXY[:,1] = ellInd[1]
    sunEll.estimate(ellXY)
    xc, yc, b, a, theta = sunEll.params
    qSunModelsRadiusViaSobel.append(b*pix)

L = 6000
srhSunVis = NP.zeros((len(fitNames),L))
srhSunRadius = NP.zeros((len(fitNames)))
qSunModelsVis = NP.zeros((len(fitNames),L))
qSunModelsRadius = NP.zeros((len(fitNames)))

uv = NP.zeros(L)
scale = 1/40
for i in range(L):
    uv[i] = i*scale
    uv_arg = 2.*NP.pi*uv[i]*FOV*NP.linspace(-.5,.5,N)
    for fitNumber in range(len(fitNames)):
        qSunModelsVis[fitNumber, i] = (NP.cos(uv_arg)*qSunModels[fitNumber]).sum()
        srhSunVis[fitNumber, i] = (NP.cos(uv_arg)*srhImages[fitNumber]).sum()

for fitNumber in range(len(fitNames)):
    qSunModelsVis[fitNumber] = NP.abs(qSunModelsVis[fitNumber]/qSunModelsVis[fitNumber].max())
    srhSunVis[fitNumber] = NP.abs(srhSunVis[fitNumber]/srhSunVis[fitNumber].max())

_A_ = 1.22
for fitNumber in range(len(fitNames)):
    qSunModelsRadius[fitNumber] = NP.rad2deg(_A_/(NP.argmin(qSunModelsVis[fitNumber])*scale))*3600/2
    srhSunRadius[fitNumber] = NP.rad2deg(_A_/(NP.argmin(srhSunVis[fitNumber])*scale))*3600/2

theta_arg = NP.rad2deg(theta_arg)*3600
radiusColors = PL.get_cmap('rainbow')
freqColors = radiusColors(NP.linspace(1.,0.,len(fitNames)))

fig = PL.figure(figsize=(12,6));
sp0 = fig.add_subplot(2,2,1,title='qSunModel')
sp1 = fig.add_subplot(2,2,2,title='srhImage')
sp2 = fig.add_subplot(2,2,3)
sp3 = fig.add_subplot(2,2,4)
sp2.set_xlim((0,L*scale))
sp3.set_xlim((0,L*scale))
sp2.set_ylim((-0.2,1.2))
sp3.set_ylim((-0.2,1.2))

for fitNumber in range(len(fitNames)):
    sp0.plot(theta_arg,qSunModels[fitNumber][N//2,:], color=freqColors[fitNumber])

    sp1.plot(theta_arg,3*srhImages[fitNumber][srh_x_size//2,:]/srhImages[fitNumber].max(),color=freqColors[fitNumber])
#    sp1.plot(theta_arg,3*srhImages[fitNumber][:,srh_y_size//2]/srhImages[fitNumber].max(),color=freqColors[fitNumber])
        
    sp1.plot(theta_arg,qSunModels[fitNumber][N//2,:], color='green')
#    sp1.plot(theta_arg,qSunModels[fitNumber][:,N//2], color='green')

    sp2.plot(uv, qSunModelsVis[fitNumber], label='%.2f, %.2f'%(qSunModelsRadius[fitNumber], qSunModelsRadiusViaSobel[fitNumber]), color=freqColors[fitNumber])
#    sp3.plot(uv, qSunModelsVis[fitNumber], label='%.2f, %.2f'%(qSunModelsRadius[fitNumber], qSunModelsRadiusViaSobel[fitNumber]), color='green')
    sp3.plot(uv, srhSunVis[fitNumber], label='%.2f, %.2f'%(srhSunRadius[fitNumber], srhSunRadiusViaSobel[fitNumber]),color=freqColors[fitNumber])

sp2.legend()
sp3.legend()
sp2.grid()
sp3.grid()
