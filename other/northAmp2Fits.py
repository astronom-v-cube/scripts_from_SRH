#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 09:43:52 2020

@author: mariagloba
"""
from astropy.io import fits
import numpy as NP
import pylab as PL
import struct
import os
import fnmatch

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

def visIndex(hh, vv):
    h = hh if hh>vv else vv
    v = vv if hh>vv else hh
    ind = (h//4 - v//4) * 16 + (2 * 52 - v//4 + 1) * (v//4) * 16 // 2 + v%4 * 4 + h%4
    return ind

def antennaIndex(ant):
    receiver1 = NP.array([81,83,85,87,89,91,95,98,101,104,107,112,116,120,124,128])
    receiver3 = NP.array([48,46,44,42,40,38,34,31,28,25,22,17,13,9,5,1])
    receiver5 = NP.array([229,231,233,235,237,239,241,243,245,247,249,251,253,255])
    receiver7 = NP.array([197,199,201,203,205,207,209,211,213,215,217,219,221,223,225,227])
    ind = []
    if ant >= 81 and ant <=128:
        ind = NP.where(receiver1 == ant)[0][0]
    elif ant >= 1 and ant <=48:
        ind = NP.where(receiver3 == ant)[0][0] + 32
    elif ant >= 229 and ant <=255:
        ind = NP.where(receiver5 == ant)[0][0] + 64
    elif ant >= 197 and ant <=227:
        ind = NP.where(receiver7 == ant)[0][0] + 96
    return ind


#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201029/vis2'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201030/vis0'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201030/vis1'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201030/vis4'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201031/vis1'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201031/vis2'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201130/vis1'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201130/vis3'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201201/vis1'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201201/vis2'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201201/vis3'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201201/vis4'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201201/vis5'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201201/vis7'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201202/vis3'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201202/vis5'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201203/vis1'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201204/vis5'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201204/vis8'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201205/vis0'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201205/vis1'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201205/vis2'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201230/vis0'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201230/vis1'
path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201230/vis2'

out_path = ''
#ant = NP.linspace(0,29,30,dtype='int')*2 + 197
ant = NP.linspace(0,29,30,dtype='int')*2 + 197

fits_list = findFits(path,'*.fit')
fits_list.sort()
dateName = fits_list[0].split('/')[-1].split('_')[1].split('.')[0]

for i in range(30):
    nfits = fits.open(fits_list[0])
    time = nfits[1].data['time']
    freqList = nfits[1].data['frequency']
    samplesNumber = time.shape[1]
    freqeuncyNumber = freqList.shape[0]
    amplcp1 = nfits[1].data['amp_lcp']
    amplitudeNumber = amplcp1.shape[1]//samplesNumber
    amplcp = amplcp1.reshape(freqeuncyNumber,samplesNumber,amplitudeNumber)[:, :, antennaIndex(ant[i])]
    
    amprcp1 = nfits[1].data['amp_rcp']
    amplitudeNumber = amprcp1.shape[1]//samplesNumber
    amprcp = amprcp1.reshape(freqeuncyNumber,samplesNumber,amplitudeNumber)[:, :, antennaIndex(ant[i])]

    for file in fits_list[1:]:
        nfits = fits.open(file)
        time = NP.concatenate((time, nfits[1].data['time']), axis = 1)
        samplesNumber = nfits[1].data['time'].shape[1]
        amplcp1 = nfits[1].data['amp_lcp']
        amplitudeNumber = amplcp1.shape[1]//samplesNumber
        amplcp = NP.concatenate((amplcp, amplcp1.reshape(freqeuncyNumber,samplesNumber,amplitudeNumber)[:, :, antennaIndex(ant[i])]),axis=1)

        amprcp1 = nfits[1].data['amp_rcp']
        amplitudeNumber = amprcp1.shape[1]//samplesNumber
        amprcp = NP.concatenate((amprcp, amprcp1.reshape(freqeuncyNumber,samplesNumber,amplitudeNumber)[:, :, antennaIndex(ant[i])]),axis=1)
        

    freqFormat = 'D';
    dataFormat = str(time.shape[1]) + 'D';
    freqColumn = fits.Column(name='frequency',format=freqFormat,array=freqList);
    timeColumn = fits.Column(name='time',format=dataFormat,array=time);
    rColumn = fits.Column(name='amp_rcp',format=dataFormat,array=amprcp);
    lColumn = fits.Column(name='amp_lcp',format=dataFormat,array=amplcp);
    
    pHeader = fits.Header();
    pHdu = fits.PrimaryHDU(header=pHeader);
    dTableHdu = fits.BinTableHDU.from_columns([freqColumn, timeColumn, rColumn, lColumn]);
    hduList = fits.HDUList([pHdu, dTableHdu]);
    fName = 'srh_amp_' + str(ant[i]) + '_' + dateName + '.fits'
    hduList.writeto(os.path.join(out_path, fName),clobber=True);
    hduList.close();
