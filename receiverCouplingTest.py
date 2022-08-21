#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 11:29:40 2021

@author: sergey_lesovoi
"""

import pylab as PL
import numpy as NP
import matplotlib

def smoothSig(signal):
    lf_filter = NP.ones(5)
    return NP.convolve(signal, lf_filter, mode = 'same') / lf_filter.shape[0]

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

times = phaseEdit.srhFits.freqTime[0]
cmap = matplotlib.cm.get_cmap('prism')
colorScale = 5
freq = 0
#N17-N32 20211015_rec8
#ampInd0 = 15        #reciever7 channel0
#visInd0 = 3022      #N16-N17

#E48-E63 20211016_rec6
#ampInd0 = 111       #reciever5 channel0
#visInd0 = 3552      #E48-E49

#E32-E47 20211016_rec5
ampInd0 = 95       #reciever4 channel0
visInd0 = 3536      #E32-E33

fig = PL.figure(figsize=(20,12))
sub = fig.subplots(nrows=2, ncols=1)
fig.suptitle('%s, %d MHz'%(phaseEdit.srhFits.dateObs, phaseEdit.srhFits.freqList[freq]/1000))

sub[0].xaxis.set_major_locator(PL.MultipleLocator(120));
sub[0].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
sub[0].xaxis.set_minor_locator(PL.MultipleLocator(30));
sub[0].set_xlabel('UT')
sub[0].set_ylim(0,1e4)
sub[0].grid()
for i in range(16): 
    sub[0].plot(times, (NP.abs(phaseEdit.srhFits.ampLcp[freq,:,i+ampInd0]) + NP.abs(phaseEdit.srhFits.ampRcp[freq,:,i+ampInd0])) + 0*i, 
                '-', color=cmap(i*colorScale), label='%s'%phaseEdit.srhFits.antennaNames[i + ampInd0])
sub[0].legend()

sub[1].xaxis.set_major_locator(PL.MultipleLocator(120));
sub[1].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
sub[1].xaxis.set_minor_locator(PL.MultipleLocator(30));
sub[1].set_xlabel('UT')
sub[1].set_yscale('log')
sub[1].grid()
for i in NP.arange(13,14):
    antA = str(phaseEdit.srhFits.antennaA[i+visInd0])
    antB = str(phaseEdit.srhFits.antennaB[i+visInd0])
    antAname = phaseEdit.srhFits.antennaNames[ NP.where(phaseEdit.srhFits.antennaNumbers == antA)[0][0]]
    antBname = phaseEdit.srhFits.antennaNames[ NP.where(phaseEdit.srhFits.antennaNumbers == antB)[0][0]]
    sub[1].plot(times, (smoothSig(NP.abs(128/2*phaseEdit.srhFits.visLcp[freq,:,i+visInd0])) + smoothSig(NP.abs(128/2*phaseEdit.srhFits.visRcp[freq,:,i+visInd0])))*1**(i+1),
                '-', color=cmap(i*colorScale),label = antAname + '-' + antBname)
sub[1].legend()
