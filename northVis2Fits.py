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


#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20200917/amps'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20200919/vis'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201003/vis'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201010/vis'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201011/vis'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201012'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201012/amps'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201016/vis'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201017/vis'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201018'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201028/amps'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201028/amps_service_point'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201029'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201029/vis0'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201029/vis2'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201030/vis4'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201126/vis1'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201127/vis2'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201127/vis3'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201128/vis1'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201128/vis5'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201201/vis1'
#path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201201/vis7'
path = '/var/run/media/svlesovoi/SRH_DATA/SRH36/20201204/vis5'

out_path = ''
#visList = [ [211,213], [215,217], [217, 219], [219,221] ]
visList = [ [211,213], [213,215], [215,217], [217,219], [211,215], [213,217], [215,219] ]
#visList = [ [197,199] ]

fits_list = findFits(path,'*.fit')
fits_list.sort()
nfits = fits.open(fits_list[0])
dateName = fits_list[0].split('/')[-1].split('_')[1].split('.')[0]

for i in range(len(visList)):
    print(visList[i])
    time = nfits[1].data['time']
    freqList = nfits[1].data['frequency']
    samplesNumber = time.shape[1]
    vis_ind = visIndex(antennaIndex(visList[i][0]), antennaIndex(visList[i][1]))
    rawlcp = nfits[1].data['vis_lcp']
    visibilityNumber = rawlcp.shape[1]//samplesNumber
    unpackFormat = ">%dl" % (rawlcp.shape[1]*2)
    ff = NP.array(struct.unpack(unpackFormat,rawlcp))
    fff = ff.reshape(rawlcp.shape[1],2)
    fff_real_lcp = NP.array(fff[:,0].reshape(1,samplesNumber,visibilityNumber))[0, :, vis_ind]
    fff_imag_lcp = NP.array(fff[:,1].reshape(1,samplesNumber,visibilityNumber))[0, :, vis_ind]
    
    rawrcp = nfits[1].data['vis_rcp']
    visibilityNumber = rawrcp.shape[1]//samplesNumber
    unpackFormat = ">%dl" % (rawrcp.shape[1]*2)
    ff = NP.array(struct.unpack(unpackFormat,rawrcp))
    fff = ff.reshape(rawrcp.shape[1],2)
    fff_real_rcp = NP.array(fff[:,0].reshape(1,samplesNumber,visibilityNumber))[0, :, vis_ind]
    fff_imag_rcp = NP.array(fff[:,1].reshape(1,samplesNumber,visibilityNumber))[0, :, vis_ind]
    
    for file in fits_list[1:]:
        nfits = fits.open(file)
        time = NP.concatenate((time, nfits[1].data['time']), axis = 1)
        samplesNumber = nfits[1].data['time'].shape[1]
        rawlcp = nfits[1].data['vis_lcp']
        visibilityNumber = rawlcp.shape[1]//samplesNumber
        unpackFormat = ">%dl" % (rawlcp.shape[1]*2)
        ff = NP.array(struct.unpack(unpackFormat,rawlcp))
        fff = ff.reshape(rawlcp.shape[1],2)
        fff_real_lcp = NP.concatenate((fff_real_lcp, fff[:,0].reshape(1,samplesNumber,visibilityNumber)[0, :, vis_ind]))
        fff_imag_lcp = NP.concatenate((fff_imag_lcp, fff[:,1].reshape(1,samplesNumber,visibilityNumber)[0, :, vis_ind]))
        rawrcp = nfits[1].data['vis_rcp']
        visibilityNumber = rawrcp.shape[1]//samplesNumber
        unpackFormat = ">%dl" % (rawrcp.shape[1]*2)
        ff = NP.array(struct.unpack(unpackFormat,rawrcp))
        fff = ff.reshape(rawrcp.shape[1],2)
        fff_real_rcp = NP.concatenate((fff_real_rcp, fff[:,0].reshape(1,samplesNumber,visibilityNumber)[0, :, vis_ind]))
        fff_imag_rcp = NP.concatenate((fff_imag_rcp, fff[:,1].reshape(1,samplesNumber,visibilityNumber)[0, :, vis_ind]))

    vis_lcp = fff_real_lcp + fff_imag_lcp * 1j
    vis_rcp = fff_real_rcp + fff_imag_rcp * 1j

    dataFormat = 'D';
    timeColumn = fits.Column(name='time',format=dataFormat,array=time[0]);
    freqColumn = fits.Column(name='frequency',format=dataFormat,array=freqList);
    dataFormat = 'C';
    rColumn = fits.Column(name='vis_rcp',format=dataFormat,array=vis_rcp);
    lColumn = fits.Column(name='vis_lcp',format=dataFormat,array=vis_lcp);

    pHeader = fits.Header();
    pHdu = fits.PrimaryHDU(header=pHeader);
    dTableHdu = fits.BinTableHDU.from_columns([freqColumn, timeColumn, rColumn, lColumn]);
    hduList = fits.HDUList([pHdu, dTableHdu]);
    fName = 'srh_vis_' + str(visList[i][0]) + '_' + str(visList[i][1]) + '_' + dateName + '.fits'
    hduList.writeto(os.path.join(out_path, fName),clobber=True);
    hduList.close();
    