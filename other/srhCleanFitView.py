#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 31 09:41:28 2021

@author: sergeyvlesovoi
"""


#from BadaryRAO import BadaryRAO
#from astropy import coordinates
#from sunpy import coordinates as sunpy_coordinates
#from skimage.transform import warp, AffineTransform
#from astropy import constants
#from zeep import Client
#from skimage import measure
import numpy as NP
import pylab as PL;
from astropy.io import fits
import os, fnmatch, sys

def arcmin_format(xy, pos):
    global gArcsecPerPixel
    global gXSize
    return '%2d' % ((xy - gXSize/2) * gArcsecPerPixel / 60);

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

#fitPath = 'SRH36_20210529_60/I/'
#fitPath = 'SRH36_20210529_60/V/'
#fitPath = 'SRH36_20210604/I/'
#fitPath = 'SRH36/cleanMaps/20210606/V/'
#fitPath = 'SRH36/cleanMaps/20210503/I/'
#fitPath = 'SRH36/cleanMaps/20210503/V/'
#fitPath = 'SRH36/cleanMaps/20210504/I/'
#fitPath = 'SRH36/cleanMaps/20210505/I/'
#fitPath = 'SRH36/cleanMaps/20210505/V/'
#fitPath = 'SRH36/cleanMaps/20210522/I/'
#fitPath = 'SRH36/cleanMaps/20210522/V/'
#fitPath = 'SRH36/cleanMaps/20210610/2800/I/'
#fitPath = 'SRH36/cleanMaps/20210629/I/'
#fitPath = 'SRH36/cleanMaps/20210629/I/'
fitPath = '/home/svlesovoi/Pictures/SRH36_pngs/20210713/II/'

fitNames =  findFits(fitPath,'*.fit')
fitNames.sort()
images = []
freqList = []
gXSize = 0
gArcsecPerPixel = 0
gDateObs = []
for fileName in fitNames:
    fd = fits.open(fileName)
    images.append(NP.flipud(fd[0].data))
    freqList.append(fd[0].header['FREQUENC'])
    gXSize = int(fd[0].header['NAXIS1'])
    gArcsecPerPixel = fd[0].header['CDELT1']
    gDateObs.append(fd[0].header['DATE-OBS'].split('.')[0])
    fd.close()
      
T0 = 2e5
levels = NP.linspace(0,5,10)*2*T0 + T0
fig, pl = PL.subplots(figsize=(6,6))
fig.tight_layout()
for i in range(len(fitNames)):
    pl.clear()
    pl.xaxis.set_major_locator(PL.MultipleLocator(128));
    pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
    pl.yaxis.set_major_locator(PL.MultipleLocator(128));
    pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
#    pl.imshow(images[i],cmap='hot',vmax=1e5,vmin=-1e5)
    pl.imshow(images[i],cmap='hot',vmax=T0,vmin=0)
    pl.contour(images[i],cmap='hot',levels=levels)
    pl.set_title('SRH %s,  %s MHz'%(gDateObs[i], freqList[i]))
    pl.set_xlabel('arcmin')
    pl.set_ylabel('arcmin')
    fig.savefig(fitNames[i].split('.')[0] + '.png')

avgImage = NP.array(images)
avgImageNumber = 25
avgImageWeek = NP.zeros((avgImageNumber,avgImage.shape[1],avgImage.shape[2]))
for i in range(avgImageNumber):
    avgImageWeek[i] = avgImage[i:i+(avgImage.shape[0] - avgImageNumber)].mean(axis=0)

#sunMask = NP.zeros_like(avgImageWeek[0])
#outRadius = 220
#sunX0 = avgImage.shape[1]/2
#sunY0 = avgImage.shape[1]/2
#for i in range(avgImage.shape[1]):
#    for j in range(avgImage.shape[1]):
#        curRadius = NP.sqrt((i - sunX0)**2 + (j - sunY0)**2)
#        if curRadius < outRadius:
#            sunMask[i,j] = 1.
#
fig, pl = PL.subplots(figsize=(6,6))
fig.tight_layout()
for i in range(avgImageNumber):
    pl.clear()
    pl.xaxis.set_major_locator(PL.MultipleLocator(128));
    pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
    pl.yaxis.set_major_locator(PL.MultipleLocator(128));
    pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
    pl.imshow(avgImageWeek[i],vmin=0,vmax=500000)
    pl.set_title('SRH %s,  %s MHz'%(gDateObs[i], freqList[i]))
    pl.set_xlabel('arcmin')
    pl.set_ylabel('arcmin')
    pl.grid(ls='dotted')
    fig.savefig('srhAvg_%s_%02d'%(gDateObs[0].split('T')[0],i) + '.png')

