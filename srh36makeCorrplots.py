#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 01:06:57 2021

@author: svlesovoi
"""
import srh36Utils
import datetime as DT;
from srhFitsFile36_amp import SrhFitsFile
from astropy.io import fits
import numpy as NP

currentDate = DT.datetime.now().date().strftime("%Y%m%d")
fitPath = '/home/svlesovoi/SRH_0306/' + currentDate + '/'

fitNames =  srh36Utils.findFits(fitPath,'*.fit')
fitNames.sort()
newestFits = SrhFitsFile(fitNames[-1],256)
freqAmount  = newestFits.freqListLength;
freqList = newestFits.freqList;
startTime= DT.datetime.strptime(newestFits.hduList[0].header['DATE-OBS'] + ' ' + newestFits.hduList[0].header['TIME-OBS'], '%Y-%m-%d %H:%M:%S');
newestTimes = newestFits.freqTime;
r_real = abs(newestFits.visRcp.real[:,:,0:3007])
r_imag = abs(newestFits.visRcp.imag[:,:,0:3007])
corrFluxRcp = NP.mean(NP.sqrt(r_real**2. + r_imag**2.),axis=2);
r_real = abs(newestFits.visLcp.real[:,:,0:3007])
r_imag = abs(newestFits.visLcp.imag[:,:,0:3007])
corrFluxLcp = NP.mean(NP.sqrt(r_real**2. + r_imag**2.),axis=2);
r_real = 0;
r_imag = 0;
ampFluxRcp = NP.mean(newestFits.ampRcp, axis = 2)
ampFluxLcp = NP.mean(newestFits.ampLcp, axis = 2)
newestFits.close()

dateName = DT.datetime.now().strftime("%Y%m%d");
try:
    currentCpFits = fits.open('srh_cp_' + dateName + '.fits')
    currentCpFreqList = currentCpFits[1].data['frequencies']
    currentCpTime = currentCpFits[2].data['time']
    currentCpCorrI = currentCpFits[2].data['I']
    currentCpCorrV = currentCpFits[2].data['V']
    currentCpFluxI = currentCpFits[2].data['flux_I']
    currentCpFluxV = currentCpFits[2].data['flux_I']
    currentCpFits.close()
except IOError:
    print('There is not cp fits')
