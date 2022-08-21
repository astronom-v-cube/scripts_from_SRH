#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  1 11:44:53 2022

@author: sergey_lesovoi
"""

from pyuvdata import UVData, utils
from srhFitsFile612 import SrhFitsFile
import numpy as NP
from astropy.time import Time
import astropy.coordinates as COORD
import base2uvw_612 as UVW

def hhmmssFromSec(sec):
  hh = int(sec / 3600.)
  sec -= hh*3600.
  mm = int(sec / 60.)
  sec -= mm*60
  ss = int(sec)
  return '%02d:%02d:%02d' % (hh,mm,ss)

#------------------------------------------------------------------------------
lastVisibility = 8192
srhFits = SrhFitsFile('srh_20211116T063224.fit',1024)
freq = 0
scan = 0
srhFits.getHourAngle(0)
srhFits.setFrequencyChannel(freq)
print('calibrating ...')
srhFits.calibrate(scan)
srhFits.vis2uv(scan)
print('centering ...')
srhFits.centerDisk()

visRcp = srhFits.visRcp[freq,scan,0:lastVisibility].copy()
visLcp = srhFits.visLcp[freq,scan,0:lastVisibility].copy()
print('correcting ...')
for vis in range(lastVisibility):
    i = vis // 128
    j = vis % 128 
    visLcp[vis] *= NP.exp(1j*(-(srhFits.ewAntPhaLcp[freq, j]+srhFits.ewLcpPhaseCorrection[freq, j]) + (srhFits.sAntPhaLcp[freq, i] + srhFits.sLcpPhaseCorrection[freq, i])))
    visRcp[vis] *= NP.exp(1j*(-(srhFits.ewAntPhaRcp[freq, j]+srhFits.ewRcpPhaseCorrection[freq, j]) + (srhFits.sAntPhaRcp[freq, i] + srhFits.sRcpPhaseCorrection[freq, i])))
    visLcp[vis] /= (srhFits.ewAntAmpLcp[freq, j] * srhFits.sAntAmpLcp[freq, i])
    visLcp[vis] /= (srhFits.ewAntAmpRcp[freq, j] * srhFits.sAntAmpRcp[freq, i])

a = 6378137.0000
b = 6356752.3141
e2 = (a*a - b*b) / (a*a)
ssrtLat = NP.deg2rad(51.759)
ssrtLon = NP.deg2rad(102.217)
ssrtHeight = 799.
v = a / (NP.sqrt(1. - e2*(NP.sin(ssrtLat) * NP.sin(ssrtLat))))
srhX = (v + ssrtHeight)*NP.cos(ssrtLat)*NP.cos(ssrtLon)
srhY = (v + ssrtHeight)*NP.cos(ssrtLat)*NP.sin(ssrtLon)
srhZ = (((1. - e2)*v) + ssrtHeight)*NP.sin(ssrtLat)

#------------------------------------------------------------------------------
UV = UVData()
UV.ant_1_array = srhFits.antennaA[0:lastVisibility]
UV.ant_2_array = srhFits.antennaB[0:lastVisibility]
UV.antenna_names = srhFits.antennaNames
UV.antenna_numbers = NP.array(list(map(int,srhFits.antennaNumbers)))
UV.Nants_data = NP.union1d(srhFits.antennaA[0:lastVisibility],srhFits.antennaB[0:lastVisibility]).size
UV.Nants_telescope = srhFits.antennaNames.shape[0]
UV.Nfreqs = 1
UV.Npols = 2
UV.Ntimes = 1
UV.Nspws = 1
UV.Nbls = lastVisibility
UV.Nblts = lastVisibility * UV.Ntimes
UV.phase_type = 'phased'
coords = COORD.get_sun(Time(srhFits.dateObs))
UV.phase_center_ra = coords.ra.rad
UV.phase_center_dec = coords.dec.rad
UV.phase_center_app_ra = NP.full(UV.Nblts,coords.ra.rad)
UV.phase_center_app_dec = NP.full(UV.Nblts,coords.dec.rad)
UV.phase_center_epoch = 2000.0
UV.channel_width = 1e7
UV.freq_array = NP.zeros((1,1))
UV.freq_array[0] = 5.8e9
UV.history = 'SRH'
UV.instrument = 'SRH'
UV.integration_time = NP.full(UV.Nblts,0.1)
UV.antenna_diameters = NP.full(UV.Nants_telescope,2.)
UV.lst_array = NP.full(UV.Nblts,0.)
UV.object_name = 'Sun'
UV.polarization_array = NP.array([-1,-2])
UV.spw_array = [1]
UV.telescope_location = [srhX, srhY, srhZ]
UV.telescope_name = 'SRH'
UV.time_array = NP.full(UV.Nblts,Time(srhFits.dateObs).jd)
UV.data_array = NP.zeros((UV.Nblts,1,UV.Nfreqs,UV.Npols),dtype='complex')
UV.data_array[:,0,0,0] = visRcp
UV.data_array[:,0,0,1] = visLcp
    
UV.flag_array = NP.full((UV.Nblts,1,UV.Nfreqs,UV.Npols),False,dtype='bool')
UV.nsample_array = NP.full((UV.Nblts,1,UV.Nfreqs,UV.Npols),1,dtype='float')
UV.vis_units = 'uncalib'

UV.antenna_positions = NP.zeros((UV.Nants_telescope,3))
for ant in NP.arange(0, 128):
#    UV.antenna_positions[ant] = [(ant - 63.5) * 4.9, 0, 0]
    UV.antenna_positions[ant] = [0, (ant - 63.5) * 4.9, 0]
#    UV.antenna_positions[ant] = [0, -(ant - 63.5) * 4.9, 0]

for ant in NP.arange(128, 192):
#    UV.antenna_positions[ant] =  [0,-(ant - 127.5) * 4.9, 0]
    UV.antenna_positions[ant] =  [-(ant - 127.5) * 4.9, 0, 0]

UV.baseline_array = 2048 * (UV.ant_2_array + 1) + UV.ant_1_array + 1 + 2**16

uvw_array = UV.antenna_positions[UV.ant_2_array] - UV.antenna_positions[UV.ant_1_array]

hA = srhFits.getHourAngle(scan)

#UV.uvw_array = utils.phase_uvw(coords.ra.rad, coords.dec.rad,uvw_array)

for vis in range(lastVisibility):
    uvw_array[vis] = UVW.base2uvw(hA,coords.dec.rad,UV.ant_2_array[vis] + 1, UV.ant_1_array[vis] + 1)

UV.uvw_array = uvw_array
#------------------------------------------------------------------------------
UV.write_uvfits('srh.uvfits',write_lst=False,spoof_nonessential=True)
