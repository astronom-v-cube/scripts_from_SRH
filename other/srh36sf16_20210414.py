#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 09:55:18 2020

@author: mariagloba
"""

import sys;
import numpy as NP
from srhFitsFile36 import SrhFitsFile
from skimage.transform import warp, AffineTransform
from astropy.io import fits
from astropy.time import Time, TimeDelta
import pylab as PL
import os
from zeep import Client
import datetime as DT
import os, fnmatch;

dateName = DT.datetime.now().strftime("%Y%m%d");

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

fitsNames = findFits('/home/svlesovoi/SRH_DATA/SRH36/' + dateName ,'*.fit')
fitsNames.sort()
fitsName = fitsNames[-1]

#imSize = 256
imSize = 200
scan = 4
frequencyChannel = 0
uvSize = 2049
arcsecPerPixel = 12
baseLines = [3,2,2,1,1,1]
allImages = NP.zeros((imSize*3, imSize*2))

srhFits = SrhFitsFile(fitsName, uvSize)
srhFits.getHourAngle(scan)
srhFits.setCalibIndex(scan)
#srhFits.useNonlinearApproach = False
srhFits.useNonlinearApproach = True
srhFits.baselines = 3

#try:
#    flagFile = open('/home/serg/SRH/Flag.txt')
#    flags = flagFile.read().split(' ')
#    if flags != '':
#        for i in range(len(flags)):
#            ant = int(flags[i])
#            if ant>=176 and ant<=192:
#                srhFits.badAntsLcp[:,192-ant] = 1
#            if ant>=49 and ant<=80:
#                srhFits.badAntsLcp[:,ant-49+16] = 1
#except:
#    pass

pAngle = 0
try:
    client = Client('http://ephemeris.rao.istp.ac.ru/?wsdl')
    result = client.service.Ephemeride('SSRT','sun',srhFits.dateObs)
    pAngle = NP.deg2rad(float(result[0]['PAngle']))
except:
    pass

for frequencyChannel in range(len(srhFits.freqList)):
    srhFits.setFrequencyChannel(frequencyChannel)
#    srhFits.baselines = baseLines[frequencyChannel]
    srhFits.vis2uv(scan);
    srhFits.uv2lmImage()
    
    scaling = srhFits.getPQScale(uvSize, NP.deg2rad(arcsecPerPixel * (uvSize - 1)/3600.)*2)
    scale = AffineTransform(scale=(uvSize/scaling[0], uvSize/scaling[1]))
    shift = AffineTransform(translation=(-uvSize/2,-uvSize/2))
    rotate = AffineTransform(rotation = -pAngle)
    matrix = AffineTransform(matrix = srhFits.getPQ2HDMatrix())
    back_shift = AffineTransform(translation=(uvSize/2,uvSize/2))
    
    lcpData = srhFits.lcp.real
    dataResult0 = warp(lcpData,(shift + (scale + back_shift)).inverse)
    lcpData = warp(dataResult0,(shift + (matrix + back_shift)).inverse)
    lcpData = warp(NP.flip(lcpData, 0),(shift + (rotate + back_shift)).inverse)
    
    x = frequencyChannel//3
    y = 2 - frequencyChannel%3
    allImages[y*imSize:(y+1)*imSize, x*imSize:(x+1)*imSize] = lcpData[uvSize//2-imSize//2:uvSize//2+imSize//2, uvSize//2-imSize//2:uvSize//2+imSize//2]

allImagesHead = NP.ones((imSize*3+25, imSize*2)) * -1000
allImagesHead[:imSize*3, :] = allImages
PL.figure(figsize=(6,10))
PL.imshow(allImagesHead, origin = 'lower', cmap = 'hot', vmin = -1, vmax = 7)
for frequencyChannel in range(len(srhFits.freqList)):
    x = frequencyChannel//3
    y = 2 - frequencyChannel%3
    PL.text(x*imSize, y*imSize+5, str(int(srhFits.freqList[frequencyChannel]/1e3)), color = 'w', fontsize = 12)
PL.text(10, imSize*3+5, 'SRH ' + srhFits.dateObs.split('T')[0] + ' ' + srhFits.dateObs.split('T')[1]  + '    MHz', color = 'w', fontsize = 15)
PL.xticks([])
PL.yticks([])
dateName = DT.datetime.now().strftime("%Y%m%d")
#pngName = ' srh_'+srhFits.dateObs+'_'+str(scan)+'_nonlin3_iter50.png'
pngName = 'srh_' + dateName + '.png'
PL.savefig(pngName, bbox_inches='tight')