#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 18:45:59 2021

@author: maria
"""

from srhFitsFile612 import SrhFitsFile
import pylab as PL
import time
import srh36Utils
import numpy as NP
from skimage.transform import warp, AffineTransform


start_time = time.time()
scan = 0
fileIndex = 10
freq = 0

#filePath = 'SRH0612/20210917/'
#fileNames = srh36Utils.findFits(filePath,'*.fit')
#fileNames.sort()
#file = SrhFitsFile(fileNames[fileIndex], 1024)
file = SrhFitsFile('SRH0612/20210917/srh_20210917T090258.fit', 1024)
file = SrhFitsFile('SRH0612/20210917/srh_20210917T055240.fit', 1024)
file = SrhFitsFile('SRH0612/20210917/srh_20210917T055215.fit', 1024)
file = SrhFitsFile('SRH0612/20210917/srh_20210917T055127.fit', 1024)
file = SrhFitsFile('SRH0612/20210917/srh_20210917T041734.fit', 1024)
file = SrhFitsFile('SRH0612/20210917/srh_20210917T041557.fit', 1024)
file.getHourAngle(scan)

file.setFrequencyChannel(freq)
file.vis2uv(scan)
#PL.figure()
#PL.imshow(NP.abs(file.uvRcp), vmax = 5e-3)

file.updateAntennaPhase(freq, baselinesNumber = 5)
file.vis2uv(scan,average=20)
file.uv2lmImage()
#PL.figure()
#PL.imshow(file.rcp.real)

scaling = file.getPQScale(file.sizeOfUv, NP.deg2rad(file.arcsecPerPixel * (file.sizeOfUv - 1)/3600.))
scale = AffineTransform(scale=(file.sizeOfUv/scaling[0], file.sizeOfUv/scaling[1]))
shift = AffineTransform(translation=(-file.sizeOfUv/2,-file.sizeOfUv/2))
rotate = AffineTransform(rotation = file.pAngle)
matrix = AffineTransform(matrix = file.getPQ2HDMatrix())
back_shift = AffineTransform(translation=(file.sizeOfUv/2,file.sizeOfUv/2))
O = file.sizeOfUv//2
Q = file.sizeOfUv//4
dataResult0 = warp(file.rcp.real,(shift + (scale + back_shift)).inverse)
srhData = warp(dataResult0,(shift + (matrix + back_shift)).inverse)
srhData = warp(srhData,(shift + (rotate + back_shift)).inverse)#[O-Q:O+Q,O-Q:O+Q]
srhData = NP.flip(srhData,0)
srhData = NP.flip(srhData,1)

PL.figure()
PL.imshow(srhData,vmin=0,vmax=3)
PL.title(file.dateObs + ' ' + str(file.freqList[freq]*1e-6) + ' GHz')

