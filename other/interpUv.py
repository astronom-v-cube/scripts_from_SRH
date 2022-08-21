#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 07:30:55 2022

@author: sergey_lesovoi
"""

import pylab as PL
import numpy as NP

iUV = NP.zeros_like(iUV)
iiUV = phaseEdit612.srhFits.uvLcp.copy()

for i in range(128):
    ii = 2*i + 386
    for j in range(128):
        jj = 2*j + 386
        iUV[ii,jj] = 0.25*(iiUV[ii-1,jj-1] + iiUV[ii-1,jj+1]+ iiUV[ii+1,jj-1]+ iiUV[ii+1,jj+1])
        
for i in range(128):
    ii = 2*i + 386
    for j in range(128):
        jj = 2*j + 385
        iUV[ii,jj] = 0.5*(iiUV[ii,jj-1] + iiUV[ii,jj+1])

for i in range(128):
    ii = 2*i + 385
    for j in range(128):
        jj = 2*j + 386
        iUV[ii,jj] = 0.5*(iiUV[ii+1,jj]+ iiUV[ii-1,jj])

# for i in range(256):
#     ii = i + 384
#     for j in range(256):
#         jj = j + 384
#         iUV[ii,jj] = iiUV[ii,jj]
#         iUV[ii,jj] += 0.4*(iiUV[ii,  jj-1] + iiUV[ii  ,jj+1])
#         iUV[ii,jj] += 0.4*(iiUV[ii-1,jj]   + iiUV[ii+1,jj])
#         iUV[ii,jj] += 0.2*(iiUV[ii+1,jj+1] + iiUV[ii-1,jj-1])
#         iUV[ii,jj] += 0.2*(iiUV[ii+1,jj-1] + iiUV[ii-1,jj+1])

# for i in range(64):
#     ii = 385 - i
#     for j in range(256 + 2*(i+1)):
#         jj = j + 384 - (i + 1)
#         iUV[ii,jj] = 0.8*iUV[ii+1,jj]
        
for L in range(32):
    side = 256 + 2*(L+1)
    r0 = 384 - (L+1) + 1
    r1 = 384 + (L+1) + 256
    c0 = 384 - (L+1)
    c1 = 384 + (L+1) + 256
    for s in range(side):
        iUV[r0,c0 + s] = 50#0.9*iUV[r0 + 1,c0 + s]
        iUV[r1,c0 + s] = 50#0.9*iUV[r1 - 1,c0 + s]
        iUV[r0,c1 + s] = 50#0.9*iUV[r1 - 1,c0 + s]
    

iImage = NP.roll(NP.fft.fft2(NP.roll(iUV,512,axis=(0,1))).real,512,axis=(0,1))
iImage /= iImage.max()
PL.figure()
PL.imshow(iImage,cmap='gray')

iiImage = NP.roll(NP.fft.fft2(NP.roll(iiUV,512,axis=(0,1))).real,512,axis=(0,1))
iiImage /= iiImage.max()
PL.figure()
PL.imshow(iiImage,cmap='gray')

# PL.figure()
# PL.imshow(NP.abs(iUV),vmax=20)
# PL.figure()
# PL.imshow(NP.abs(iiUV[350:650,350:650]),vmax=100)

