#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 01:58:46 2020

@author: svlesovoi
"""

import numpy as NP
import pylab as PL
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

def hhmm_format(t, pos):
  t = 0.1*t + _timeStart
  hh = int(t / 3600.)
  t -= hh*3600.
  mm = int(t / 60.)
  t -= mm*60
  ss = int(t)
  return '%02d:%02d:%02d' % (hh,mm,ss)

class VisibilityMeasurement():

    def __init__(self, parent=None):
        parser = OptionParser();
        parser.add_option("-d", "--date", dest="fitsDate",   default="2016-07-01");
        parser.add_option("-t", "--time", dest="fitsTime",   default="00:00:00");
        (cl_options,cl_args) = parser.parse_args();
        self.fitsDate = cl_options.fitsDate;
        self.fitsTime = cl_options.fitsTime;
        print('srh_' + re.sub(r'\W+','',self.fitsDate) + 'T' + re.sub(r'\W+','',self.fitsTime) + '.fit')
        listTime = self.fitsTime.split(':')
        print(int(listTime[0])*3600 + int(listTime[1])*60 + int(listTime[2]))
        self.sunFlux = 75
        self.LF_size = 13

    def show(self, antInd, ampInd, nameSuffix):
        N = self.ampLcp[0,:,antInd].shape[0]
        self.lf_filter = NP.ones(self.LF_size)
        self.ampInput = self.ampLcp[0,:,ampInd].copy()
        minValue = NP.mean(self.ampInput[0:N//4])
        maxValue = NP.mean(self.ampInput[-N//4:])
        self.ampInput = (self.ampInput - minValue)/(maxValue - minValue)*self.sunFlux
        self.ampSmoothed = NP.convolve(self.ampInput, self.lf_filter, mode='same') / self.LF_size
        self.noise = self.ampInput[self.LF_size//2:-self.LF_size//2] - self.ampSmoothed[self.LF_size//2:-self.LF_size//2]
        self.noiseSize = self.noise.shape[0]
        self.dNoiseSize = self.noiseSize//3
        if self.dNoiseSize % 2:
            self.dNoiseSize += 1
        self.skyNoise = NP.std(self.noise[0:self.noiseSize//3])
        self.sunNoise = NP.std(self.noise[-self.noiseSize//3:])
        PL.clf()
        pl0 = PL.subplot(111)
        pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
        pl0.plot(10*self.noise,label='10'+ r'$\times$'+ 'noise', color='gray')
        pl0.plot(10*self.noise[0:self.dNoiseSize], label='sky %.3f s.f.u.'%self.skyNoise, color='blue')
        pl0.plot(NP.linspace(self.noiseSize - self.dNoiseSize,self.noiseSize-1,self.dNoiseSize),10*self.noise[-self.dNoiseSize:],label='Sun %.3f s.f.u.'%self.sunNoise,color='orange')
        pl0.plot(self.ampSmoothed[self.LF_size//2:-self.LF_size//2], label = 'sky->Sun response',color='black')
        pl0.grid()
        pl0.set_xlabel('UT')
        pl0.set_ylabel('s.f.u')
        pl0.set_title(label='%s, %s MHz, 2.8 m antenna %d, feed %d, FrontEnd %d'%(nameSuffix.split('_')[0], nameSuffix.split('_')[1], self.antNames[antInd], self.antFeeds[antInd], self.antFrontEnds[antInd]))
        pl0.legend()
        PL.savefig('snr_%d_%s.png'%(self.antNames[antInd], nameSuffix))
        
    def openFits(self, path, minMinutes, maxMinutes):
        fitsNames = findFits(path,'*.fit',minMinutes = minMinutes, maxMinutes = maxMinutes)
        print(fitsNames)
        fitsNames.sort()
        
        if fitsNames[0]:
            nfits = fits.open(fitsNames[0])
            self.time = nfits[1].data['time']
            global _timeStart
            _timeStart = self.time[0][0]
            samplesNumber = self.time.shape[1]
            vislcp = nfits[1].data['vis_lcp']
            self.visListLength = nfits[1].data['vis_lcp'].size // samplesNumber
            visrcp = nfits[1].data['vis_rcp']
            self.visLcp = NP.reshape(nfits[1].data['vis_lcp'],(samplesNumber, self.visListLength));
            self.visRcp = NP.reshape(nfits[1].data['vis_rcp'],(samplesNumber, self.visListLength));
            amplcp = nfits[1].data['amp_lcp']
            amplitudeNumber = amplcp.shape[1]//samplesNumber
            self.ampLcp = amplcp.reshape(1,samplesNumber,amplitudeNumber)
            amprcp = nfits[1].data['amp_rcp']
            self.ampRcp = amprcp.reshape(1,samplesNumber,amplitudeNumber)
            self.antNames = nfits[2].data['ant_name']
            self.antFeeds = nfits[2].data['ant_feed_id']
            self.antFrontEnds = nfits[2].data['ant_fe_id']
        
            for fitsName in fitsNames[1:]:
                nfits = fits.open(fitsName)
                self.time = NP.concatenate((self.time, nfits[1].data['time']), axis = 1)
                samplesNumber = nfits[1].data['time'].shape[1]
                vislcp = NP.reshape(nfits[1].data['vis_lcp'],(samplesNumber, self.visListLength));
                visrcp = NP.reshape(nfits[1].data['vis_rcp'],(samplesNumber, self.visListLength));
                self.visLcp = NP.concatenate((self.visLcp, vislcp), axis = 1)
                self.visRcp = NP.concatenate((self.visRcp, visrcp), axis = 1)
                amplcp = nfits[1].data['amp_lcp']
                amplitudeNumber = amplcp.shape[1]//samplesNumber
                self.ampLcp = NP.concatenate((self.ampLcp, amplcp.reshape(1,samplesNumber,amplitudeNumber)), axis = 1)
                amprcp = nfits[1].data['amp_rcp']
                self.ampRcp = NP.concatenate((self.ampRcp, amprcp.reshape(1,samplesNumber,amplitudeNumber)), axis = 1)

visPsf = VisibilityMeasurement()
#visPsf.openFits('/var/run/media/svlesovoi/SRH_DATA/SRH36/20200919/', minMinutes = 33, maxMinutes = 48)
visPsf.openFits('/var/run/media/svlesovoi/SRH_DATA/SRH36/20201002/amps/', minMinutes = 490, maxMinutes = 494)

PL.figure(figsize=(10,10))
antenna = 16
antFitsInd = ((antenna // 16) * 32) + (antenna % 16)
visPsf.show(antenna, antFitsInd, '20200918_4500')
