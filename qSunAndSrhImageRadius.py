#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 01:20:56 2018

@author: sergey
"""
from astropy.io import fits
import pylab as PL
import numpy as NP
import scipy.signal
import datetime as DT
import ephem
import base2uvw as b2u
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
#arcsecRadius = 990 #lastFreq
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

for i in range(N):
    for j in range(N):
        x=i - N/2
        y=j - N/2
#        if (NP.sqrt(x**2 + y**2) < 1.0*radius and NP.sqrt(x**2 + y**2) > 0.91*radius and NP.abs(y) > NP.abs(1.6*x)):
        if (NP.sqrt(x**2 + y**2) < 1.0*radius and NP.sqrt(x**2 + y**2) > 0.95*radius and NP.abs(y) > NP.abs(1.6*x)):
            qRing[i, j] = .85
        if (NP.sqrt(x**2 + y**2) < radius):
            qSun[i, j] = 1.
            
#qSun[int(N/2)-5:int(N/2)+5,174:182] += 10
qSunPB = qSun*primaryBeam/edgeBright
            
#dL = 1*( 25//2) + 1
dL = 1*( 19//2) + 1
kern = NP.ones((dL,dL))
cqSun = scipy.signal.fftconvolve(qSun,kern) / dL**2
cqSun = cqSun[dL//2:dL//2+N,dL//2:dL//2+N]
cqRing = scipy.signal.fftconvolve(qRing,kern) / dL**2
cqRing = cqRing[dL//2:dL//2+N,dL//2:dL//2+N]

#cqSunPB = cqSun*primaryBeam/edgeBright
cqSunPB = (cqSun + cqRing)*primaryBeam

qSunEdge = filters.sobel(qSun)
cqSunEdge = filters.sobel(cqSun)
qSunPBEdge = filters.sobel(qSunPB)
cqSunPBEdge = filters.sobel(cqSunPB)

sunEll = skimage.measure.EllipseModel()

srhFitsPath = '/home/svlesovoi/Documents/Python Scripts/srhFits'
fitNames = findFits(srhFitsPath, '*.fits')
srhImages = []
srhEdges = []
srhSunRadiusViaSobel = []
qSunModelsRadiusViaSobel = []
srhFreqs = []

qSunModels = []
for i in range(len(fitNames)):
    qSun[:, :] = 0.
    qRing[:, :] = 0.
    arcsecRadius = 1040 - i*20
    degRadius = NP.deg2rad(arcsecRadius/3600)
    radius = int(arcsecRadius/pix +0.5)
    for i in range(N):
        for j in range(N):
            x=i - N/2
            y=j - N/2
            if (NP.sqrt(x**2 + y**2) < 1.0*radius and NP.sqrt(x**2 + y**2) > 0.95*radius and NP.abs(y) > NP.abs(1.6*x)):
                qRing[i, j] = .85
            if (NP.sqrt(x**2 + y**2) < radius):
                qSun[i, j] = 1.
    cqSun = scipy.signal.fftconvolve(qSun,kern) / dL**2
    cqSun = cqSun[dL//2:dL//2+N,dL//2:dL//2+N]
    cqRing = scipy.signal.fftconvolve(qRing,kern) / dL**2
    cqRing = cqRing[dL//2:dL//2+N,dL//2:dL//2+N]
    qSunModels.append((cqSun + cqRing)*primaryBeam)
    
for fitName in fitNames:
    srhFits = fits.open(fitName)
    srhFreqs.append(float(srhFits[0].header['FREQUENC']))
    srh_x_size = srhFits[0].header['NAXIS1']
    srh_y_size = srhFits[0].header['NAXIS2']
    srh_x_delt = srhFits[0].header['CDELT1'] / 60 * 2
    srh_y_delt = srhFits[0].header['CDELT2'] / 60 * 2
    srhFitsData = srhFits[0].data
    srhImages.append(srhFitsData + NP.abs(srhFitsData))
    x0 = int(srh_x_size//2)
    y0 = int(srh_y_size//2)
    dx = int(srh_x_size//4)
    dy = int(srh_y_size//4)
    srhFitsDataMean = NP.mean(srhFitsData[x0 - dx:x0 + dx, y0 - dy:y0 + dy])
    srhFitsEdge = filters.sobel(NP.clip(srhFitsData,srhFitsDataMean/3 - srhFitsDataMean/100,srhFitsDataMean/3 + srhFitsDataMean/100))
    srhEdges.append(srhFitsEdge)
    ellInd = NP.where(srhFitsEdge > .5*srhFitsEdge.max())
    ellXY = NP.zeros((len(ellInd[1]),2))
    ellXY[:,0] = ellInd[1]
    ellXY[:,1] = ellInd[0]
    ellXY[:,0] = (ellXY[:,0] - x0) * srhFits[0].header['CDELT1'] / 60 * 2
    ellXY[:,1] = (ellXY[:,1] - y0) * srhFits[0].header['CDELT2'] / 60 * 2
    sunEll.estimate(ellXY)
    xc, yc, a, b, theta = sunEll.params
    srhSunRadiusViaSobel.append(b*60)


for fitNumber in range(len(fitNames)):
    qSunModelsEdge = filters.sobel(qSunModels[fitNumber])
    ellInd = NP.where(qSunModelsEdge > NP.max(qSunModelsEdge)*.5)
    ellXY = NP.zeros((len(ellInd[1]),2))
    ellXY[:,0] = ellInd[0]
    ellXY[:,1] = ellInd[1]
    sunEll.estimate(ellXY)
    xc, yc, b, a, theta = sunEll.params
    qSunModelsRadiusViaSobel.append(b*pix)

ellInd = NP.where(qSunPBEdge > NP.max(qSunPBEdge)*.5)
ellXY = NP.zeros((len(ellInd[1]),2))
ellXY[:,0] = ellInd[0]
ellXY[:,1] = ellInd[1]
sunEll.estimate(ellXY)
xc, yc, b, a, theta = sunEll.params
qSunPBRadiusViaSobel = b*pix

ellInd = NP.where(cqSunEdge > NP.max(cqSunEdge)*.5)
ellXY = NP.zeros((len(ellInd[1]),2))
ellXY[:,0] = ellInd[0]
ellXY[:,1] = ellInd[1]
sunEll.estimate(ellXY)
xc, yc, b, a, theta = sunEll.params
cqSunRadiusViaSobel = b*pix

ellInd = NP.where(cqSunPBEdge > NP.max(cqSunPBEdge)*.5)
ellXY = NP.zeros((len(ellInd[1]),2))
ellXY[:,0] = ellInd[0]
ellXY[:,1] = ellInd[1]
sunEll.estimate(ellXY)
xc, yc, b, a, theta = sunEll.params
cqSunPBRadiusViaSobel = b*pix

L = 6000
qSunVis = NP.zeros(L)
cqSunVis = NP.zeros(L)
qSunPBVis = NP.zeros(L)
cqSunPBVis = NP.zeros(L)
srhSunVis = NP.zeros((len(fitNames),L))
srhSunRadius = NP.zeros((len(fitNames)))
qSunModelsVis = NP.zeros((len(fitNames),L))
qSunModelsRadius = NP.zeros((len(fitNames)))

uv = NP.zeros(L)
scale = 1/40
for i in range(L):
    uv[i] = i*scale
    uv_arg = 2.*NP.pi*uv[i]*FOV*NP.linspace(-.5,.5,N)
    cqSunVis[i] = (NP.cos(uv_arg)*cqSun).sum()
    cqSunPBVis[i] = (NP.cos(uv_arg)*cqSunPB).sum()
    for fitNumber in range(len(fitNames)):
        qSunModelsVis[fitNumber, i] = (NP.cos(uv_arg)*qSunModels[fitNumber]).sum()
        srhSunVis[fitNumber, i] = (NP.cos(uv_arg)*srhImages[fitNumber]).sum()

cqSunVis = NP.abs(cqSunVis/cqSunVis.max())
cqSunPBVis = NP.abs(cqSunPBVis/cqSunPBVis.max())
for fitNumber in range(len(fitNames)):
    qSunModelsVis[fitNumber] = NP.abs(qSunModelsVis[fitNumber]/qSunModelsVis[fitNumber].max())
    srhSunVis[fitNumber] = NP.abs(srhSunVis[fitNumber]/srhSunVis[fitNumber].max())

_A_ = 1.22
cqSunRadius = NP.rad2deg(_A_/(NP.argmin(cqSunVis)*scale))*3600/2
cqSunPBRadius = NP.rad2deg(_A_/(NP.argmin(cqSunPBVis)*scale))*3600/2
for fitNumber in range(len(fitNames)):
    qSunModelsRadius[fitNumber] = NP.rad2deg(_A_/(NP.argmin(qSunModelsVis[fitNumber])*scale))*3600/2
    srhSunRadius[fitNumber] = NP.rad2deg(_A_/(NP.argmin(srhSunVis[fitNumber])*scale))*3600/2

theta_arg = NP.rad2deg(theta_arg)*3600
radiusColors = PL.get_cmap('rainbow')
freqColors = radiusColors(NP.linspace(1.,0.,len(fitNames)))

fig = PL.figure(figsize=(12,6));
sp0 = fig.add_subplot(2,2,1,title='qSun')
sp1 = fig.add_subplot(2,2,2,title='qSunPB')
sp2 = fig.add_subplot(2,2,3)
sp3 = fig.add_subplot(2,2,4)
sp2.set_xlim((0,L*scale))
sp3.set_xlim((0,L*scale))
sp2.set_ylim((-0.2,1.2))
sp3.set_ylim((-0.2,1.2))

sp0.plot(theta_arg,cqSun[N//2,:], color='red')
sp0.plot(theta_arg,cqSunPB[N//2,:], color='green')
sp0.plot([cqSunRadius,cqSunRadius],[0,1], color='red')
sp0.plot([cqSunPBRadius,cqSunPBRadius],[0,1], color='green')

#for fitNumber in range(len(fitNames)):
for fitNumber in range(1):
#fitNumber = 0
    sp1.plot(theta_arg,3*srhImages[fitNumber][srh_x_size//2,:]/srhImages[fitNumber].max(),color=freqColors[fitNumber])
    sp1.plot(theta_arg,3*srhImages[fitNumber][:,srh_y_size//2]/srhImages[fitNumber].max(),color=freqColors[fitNumber])
    sp1.plot([srhSunRadius[fitNumber],srhSunRadius[fitNumber]],[0,1],color=freqColors[fitNumber])
    #    sp1.plot([srhSunRadiusViaSobel[fitNumber],srhSunRadiusViaSobel[fitNumber]],[0,.5],color=freqColors[fitNumber])
        
    sp1.plot(theta_arg,cqSunPB[N//2,:], color='green')
    sp1.plot(theta_arg,cqSunPB[:,N//2], color='green')
    sp1.plot([cqSunPBRadius,cqSunPBRadius],[0,1],color='green')

sp2.plot(uv, cqSunVis, label='%.2f, %.2f'%(cqSunRadius, cqSunRadiusViaSobel), color='red')
sp2.plot(uv, cqSunPBVis, label='%.2f, %.2f'%(cqSunPBRadius, cqSunPBRadiusViaSobel), color='green')

#for fitNumber in range(len(fitNames)):
for fitNumber in range(1):
#    sp2.plot(uv, srhSunVis[fitNumber], label='%.2f, %.2f'%(srhSunRadius[fitNumber], srhSunRadiusViaSobel[fitNumber]))

#for fitNumber in range(len(fitNames)):
#for fitNumber in range(1):
    sp3.plot(uv, srhSunVis[fitNumber], label='%.2f, %.2f'%(srhSunRadius[fitNumber], srhSunRadiusViaSobel[fitNumber]),color=freqColors[fitNumber])

sp3.plot(uv, cqSunPBVis, label='%.2f, %.2f'%(cqSunPBRadius, cqSunPBRadiusViaSobel),color='green')

sp2.legend()
sp3.legend()
sp2.grid()
sp3.grid()
