#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 06:58:12 2020

@author: svlesovoi
"""

import pylab as PL
import numpy as NP
from optparse import OptionParser
import re
import os, fnmatch;
from astropy.io import fits

def findFits(path, pattern, minMinutes, maxMinutes):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    timeNum = int(basename.split('T')[1].split('.')[0])
                    nameHh = timeNum // 10000
                    nameMm = (timeNum - nameHh*10000) // 100
                    nameMinutes = nameHh*60 + nameMm
                    if nameMinutes >= minMinutes and nameMinutes <= maxMinutes:
                        result.append(os.path.join(root,basename))
    return result

#fitsNames = findFits('/var/run/media/svlesovoi/SRH_DATA/SRH36/20200917/amps/','*.fit',minMinutes = 500, maxMinutes = 507)
#fitsNames = findFits('/var/run/media/svlesovoi/SRH_DATA/SRH36/20200918/amps/','*.fit',minMinutes = 500, maxMinutes = 503)
#fitsNames = findFits('/var/run/media/svlesovoi/SRH_DATA/SRH36/20200919/amps/','*.fit',minMinutes = 490, maxMinutes = 496)
#fitsNames = findFits('/var/run/media/svlesovoi/SRH_DATA/SRH36/20200920/amps/','*.fit',minMinutes = 145, maxMinutes = 152)
#fitsNames = findFits('/var/run/media/svlesovoi/SRH_DATA/SRH36/20201002/amps/','*.fit', minMinutes = 418, maxMinutes = 422)
#fitsNames = findFits('/var/run/media/svlesovoi/SRH_DATA/SRH36/20201002/amps/','*.fit', minMinutes = 490, maxMinutes = 494)
#fitsNames = findFits('/var/run/media/svlesovoi/SRH_DATA/SRH36/20201003/amps/','*.fit', minMinutes = 52, maxMinutes = 56)
fitsNames = findFits('/var/run/media/svlesovoi/SRH_DATA/SRH36/20201003/amps/','*.fit', minMinutes = 178, maxMinutes = 184)
fitsNames.sort()
#nameSuffix = '20200917_4500'
#nameSuffix = '20200920_5900'
#nameSuffix = '20201002_3000'
#nameSuffix = '20201002_4500'
nameSuffix = '20201003_5000_1'

if fitsNames[0]:
    nfits = fits.open(fitsNames[0])
    time = nfits[1].data['time']
    samplesNumber = time.shape[1]
    amplcp = nfits[1].data['amp_lcp']
    amplitudeNumber = amplcp.shape[1]//samplesNumber
    ampLcp = amplcp.reshape(1,samplesNumber,amplitudeNumber)
    amprcp = nfits[1].data['amp_rcp']
    ampRcp = amprcp.reshape(1,samplesNumber,amplitudeNumber)

    for fitsName in fitsNames[1:]:
        nfits = fits.open(fitsName)
        time = NP.concatenate((time, nfits[1].data['time']), axis = 1)
        samplesNumber = nfits[1].data['time'].shape[1]
        amplcp = nfits[1].data['amp_lcp']
        amplitudeNumber = amplcp.shape[1]//samplesNumber
        ampLcp = NP.concatenate((ampLcp, amplcp.reshape(1,samplesNumber,amplitudeNumber)), axis = 1)
        amprcp = nfits[1].data['amp_rcp']
        ampRcp = NP.concatenate((ampRcp, amprcp.reshape(1,samplesNumber,amplitudeNumber)), axis = 1)

sunWidth = 36/60
psfWidths = []
PL.figure(figsize=(10,10))
for antenna in range(48):
    antFitsInd = ((antenna // 16) * 32) + (antenna % 16)
    ampInput = ampLcp[0,:,antFitsInd]
    ampInput -= ampInput[0]
    ampInput = ampInput / ampInput.max()
    ampTime = time[0]
    ampAngleDegree = (ampTime - ampTime[0])*(9*15)/3600
    psfInd = NP.where(ampInput > .5)
    responseWidth = ampAngleDegree[psfInd[0][-1]] - ampAngleDegree[psfInd[0][0]]
    psfWidth = NP.sqrt(responseWidth**2 - sunWidth**2)
    psfWidths.append(psfWidth)
    PL.clf()
    PL.plot(ampAngleDegree,ampInput,label='antenna %d, feed %d, FrontEnd %d'%(nfits[2].data['ant_name'][antenna], nfits[2].data['ant_feed_id'][antenna], nfits[2].data['ant_fe_id'][antenna]))
    PL.plot(ampAngleDegree[psfInd],ampInput[psfInd])
    PL.plot([ampAngleDegree[psfInd[0][0]], ampAngleDegree[psfInd[0][0]]],[0,1], color='black', linewidth=.5)
    PL.plot([ampAngleDegree[psfInd[0][-1]], ampAngleDegree[psfInd[0][-1]]],[0,1], color='black', linewidth=.5)
    PL.plot([ampAngleDegree[0], ampAngleDegree[-1]],[0.5,.5], color='black', linewidth=.5)
    PL.grid()
    PL.xlabel('degree')
    PL.ylabel('arbitrary')
    PL.title('%s, %s MHz, 2.8 m, FWHM %.2f%s'%(nameSuffix.split('_')[0], nameSuffix.split('_')[1], psfWidth,'$^\circ$'))
    PL.ylim(-0.1,1.1)
    PL.legend(loc='upper left')
    PL.savefig('psf_%d_%s.png'%(nfits[2].data['ant_name'][antenna], nameSuffix))
