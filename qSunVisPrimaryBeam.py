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

N = 256
L = 300
pix=20 # arcsec per pixel
arcsecRadius = 1080
degRadius = NP.deg2rad(arcsecRadius/3600)
radius = int(arcsecRadius/pix +0.5)
FOV = NP.deg2rad(pix * N / 3600)

qSun = NP.zeros((N, N))

theta_arg = FOV*NP.linspace(-.5,.5,N)
gx, gy = NP.meshgrid(theta_arg,theta_arg)
HBW = NP.deg2rad(2.)# half beam width [degree]
EBW = NP.deg2rad(.5)# half beam width [degree]
primaryBeam = NP.exp(-((gx/HBW)**2 + (gy/HBW)**2))
edgeBright = NP.exp(-((gx/EBW)**2 + (gy/EBW)**2))

for i in range(N):
    for j in range(N):
        x=i - N/2
        y=j - N/2
        if (NP.sqrt(x**2 + y**2) < radius and NP.sqrt(x**2 + y**2) > 0.95*radius and NP.abs(y) > NP.abs(2*x)):
            qSun[i, j] = 1.
        elif (NP.sqrt(x**2 + y**2) < radius):
            qSun[i, j] = 1.
            
qSun[int(N/2)-5:int(N/2)+5,174:182] += 10
qSunPB = qSun*primaryBeam#/edgeBright
            
dL = 1*( 15//2) + 1
kern = NP.ones((dL,dL))
cqSun = scipy.signal.fftconvolve(qSun,kern) / dL**2
cqSun = cqSun[dL//2:dL//2+N,dL//2:dL//2+N]

cqSunPB = cqSun*primaryBeam#/edgeBright

qSunEdge = filters.sobel(qSun)
cqSunEdge = filters.sobel(cqSun)
qSunPBEdge = filters.sobel(qSunPB)
cqSunPBEdge = filters.sobel(cqSunPB)

sunEll = skimage.measure.EllipseModel()

ellInd = NP.where(qSunEdge > NP.max(qSunEdge)*.5)
ellXY = NP.zeros((len(ellInd[1]),2))
ellXY[:,0] = ellInd[0]
ellXY[:,1] = ellInd[1]
sunEll.estimate(ellXY)
xc, yc, b, a, theta = sunEll.params
qSunRadiusViaSobel = b*pix

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
uv = NP.zeros(L)
scale = 1/40
for i in range(L):
    uv[i] = i*scale
    uv_arg = 2.*NP.pi*uv[i]*FOV*NP.linspace(-.5,.5,N)
    qSunVis[i] = (NP.cos(uv_arg)*qSun).sum()
    cqSunVis[i] = (NP.cos(uv_arg)*cqSun).sum()
    qSunPBVis[i] = (NP.cos(uv_arg)*qSunPB).sum()
    cqSunPBVis[i] = (NP.cos(uv_arg)*cqSunPB).sum()

qSunVis = NP.abs(qSunVis/qSunVis.max())
cqSunVis = NP.abs(cqSunVis/cqSunVis.max())
qSunPBVis = NP.abs(qSunPBVis/qSunPBVis.max())
cqSunPBVis = NP.abs(cqSunPBVis/cqSunPBVis.max())

qSunRadius = NP.rad2deg(1.22/(NP.argmin(qSunVis)*scale))*3600/2
qSunPBRadius = NP.rad2deg(1.22/(NP.argmin(qSunPBVis)*scale))*3600/2
cqSunRadius = NP.rad2deg(1.22/(NP.argmin(cqSunVis)*scale))*3600/2
cqSunPBRadius = NP.rad2deg(1.22/(NP.argmin(cqSunPBVis)*scale))*3600/2

theta_arg = NP.rad2deg(theta_arg)*3600

fig = PL.figure(figsize=(12,6));
sp0 = fig.add_subplot(2,2,1,title='qSun')
sp1 = fig.add_subplot(2,2,2,title='qSunPB')
sp2 = fig.add_subplot(2,2,3)
sp3 = fig.add_subplot(2,2,4)
sp2.set_xlim((0,L*scale))
sp3.set_xlim((0,L*scale))
sp2.set_ylim((-0.2,1.2))
sp3.set_ylim((-0.2,1.2))

#sp0.imshow(qSun)
#sp1.imshow(qSunPB)
sp0.plot(theta_arg,qSun[N//2,:])
sp0.plot(theta_arg,cqSun[N//2,:])
sp0.plot([qSunRadius,qSunRadius],[0,1])
sp1.plot(theta_arg,qSunPB[N//2,:])
sp1.plot(theta_arg,cqSunPB[N//2,:])
sp1.plot([cqSunRadius,cqSunRadius],[0,1])

sp2.plot(uv, qSunVis, label='%.2f, %.2f'%(qSunRadius, qSunRadiusViaSobel))
sp2.plot(uv, cqSunVis, label='%.2f, %.2f'%(cqSunRadius, cqSunRadiusViaSobel))
sp3.plot(uv, qSunPBVis, label='%.2f, %.2f'%(qSunPBRadius, qSunPBRadiusViaSobel))
sp3.plot(uv, cqSunPBVis, label='%.2f, %.2f'%(cqSunPBRadius, cqSunPBRadiusViaSobel))

sp2.legend()
sp3.legend()
sp2.grid()
sp3.grid()
