#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  1 11:44:53 2022

@author: sergey_lesovoi
"""

from pyuvdata import UVData
from srhFitsFile36 import SrhFitsFile
import numpy as NP
from astropy.time import Time
import astropy.coordinates as COORD
from astropy.coordinates import EarthLocation
from astropy import units as u
from datetime import timedelta

class Srh0306UVData(UVData):
    """
    """
    def __init__(self):
        self.srhLat = 51.759
        self.srhLon = 102.217
        self.srhAlt = 799.
        self.srhEarthLocation = EarthLocation(lon=self.srhLon * u.deg, lat=self.srhLat * u.deg, height=self.srhAlt * u.m)
        self.srhLocation = self.srhEarthLocation.T.to_value()
        self.srhLat = NP.deg2rad(self.srhLat)
        self.srhLon = NP.deg2rad(self.srhLon)
        self.phi_operator = NP.array([
            [-NP.sin(self.srhLat), 0., NP.cos(self.srhLat)],
            [0., 1., 0.],
            [NP.cos(self.srhLat), 0., NP.sin(self.srhLat)]
            ])
        super().__init__()
        
    def baseline2uvw(self, hourAngle, declination, antenna0, antenna1):
        if ((antenna0 >= 1 and antenna0 <= 97) and (antenna1 >= 98 and antenna1 <= 128)):
            base = NP.array([-antenna1 + 97, antenna0 - 33, 0.])
        elif ((antenna1 >= 1 and antenna1 <= 97) and (antenna0 >= 98 and antenna0 <= 128)):
            base = NP.array([-antenna0 + 97, antenna1 - 33, 0.])
        elif ((antenna0 >= 1 and antenna0 <= 97) and (antenna1 >= 1 and antenna1 <= 97)):
            base = NP.array([0.,antenna0 - antenna1,0.])
        elif ((antenna0 >= 98 and antenna0 <= 128) and (antenna1 >= 98 and antenna1 <= 128)):
            base = NP.array([-antenna0 + antenna1,0.,0.])
        
        base *= 9.8;

        uvw_operator = NP.array([
            [ NP.sin(hourAngle),		 NP.cos(hourAngle),		0.	  ],
            [-NP.sin(declination)*NP.cos(hourAngle),  NP.sin(declination)*NP.sin(hourAngle), NP.cos(declination)], 
            [ NP.cos(declination)*NP.cos(hourAngle), -NP.cos(declination)*NP.sin(hourAngle), NP.sin(declination)]  
            ])
    

        return NP.dot(uvw_operator, NP.dot(self.phi_operator, base))

    def write_uvfits(self, srcName, dstName, lastVisibility=3007, frequency=0, scan=0, calibrating=True, uvSize=1024, averaging=0):
        srhFits = SrhFitsFile(srcName, uvSize)
        srhFits.setFrequencyChannel(frequency)
        srhFits.getHourAngle(scan)
        if (calibrating):
            srhFits.calibrate(frequency)
        srhFits.vis2uv(scan, average=averaging)
        visRcp = srhFits.visRcp[frequency,scan,0:lastVisibility].copy()
        visLcp = srhFits.visLcp[frequency,scan,0:lastVisibility].copy()
        scanTime = Time(srhFits.dateObs.split('T')[0] + 'T' + str(timedelta(seconds=srhFits.freqTime[frequency,scan])))
        
        coords = COORD.get_sun(scanTime)
        
        for vis in range(lastVisibility):
            i = vis // 97
            j = vis % 97 
            visLcp[vis] *= NP.exp(1j*(-(srhFits.ewAntPhaLcp[frequency, j]+srhFits.ewLcpPhaseCorrection[frequency, j]) + (srhFits.nAntPhaLcp[frequency, i] + srhFits.nLcpPhaseCorrection[frequency, i])))
            visRcp[vis] *= NP.exp(1j*(-(srhFits.ewAntPhaRcp[frequency, j]+srhFits.ewRcpPhaseCorrection[frequency, j]) + (srhFits.nAntPhaRcp[frequency, i] + srhFits.nRcpPhaseCorrection[frequency, i])))
            visLcp[vis] /= (srhFits.ewAntAmpLcp[frequency, j] * srhFits.nAntAmpLcp[frequency, i])
            visRcp[vis] /= (srhFits.ewAntAmpRcp[frequency, j] * srhFits.nAntAmpRcp[frequency, i])
        
        self.ant_1_array = srhFits.antennaA[0:lastVisibility]
        self.ant_2_array = srhFits.antennaB[0:lastVisibility]
        self.antenna_names = srhFits.antennaNames
        self.antenna_numbers = NP.array(list(map(int,srhFits.antennaNumbers)))
        self.Nants_data = NP.union1d(srhFits.antennaA[0:lastVisibility],srhFits.antennaB[0:lastVisibility]).size
        self.Nants_telescope = srhFits.antennaNames.shape[0]
        self.Nfreqs = 1
        self.Npols = 2
        self.Ntimes = 1
        self.Nspws = 1
        self.Nbls = lastVisibility
        self.Nblts = lastVisibility * self.Ntimes
        self.phase_type = 'phased'
        self.phase_center_ra = coords.ra.rad
        self.phase_center_dec = coords.dec.rad
        self.phase_center_app_ra = NP.full(self.Nblts,coords.ra.rad)
        self.phase_center_app_dec = NP.full(self.Nblts,coords.dec.rad)
        self.phase_center_epoch = 2000.0
        self.channel_width = 1e7
        self.freq_array = NP.zeros((1,1))
        self.freq_array[0] = srhFits.freqList[frequency]*1e3
        self.history = 'SRH'
        self.instrument = 'SRH0306'
        self.integration_time = NP.full(self.Nblts,0.1)
        self.antenna_diameters = NP.full(self.Nants_telescope,3.)
        self.lst_array = NP.full(self.Nblts,0.)
        self.object_name = 'Sun'
        self.polarization_array = NP.array([-1,-2])
        self.spw_array = [1]
        self.telescope_location = list(self.srhLocation)
        self.telescope_name = 'SRH'
        self.time_array = NP.full(self.Nblts,scanTime.jd)
        self.data_array = NP.zeros((self.Nblts,1,self.Nfreqs,self.Npols),dtype='complex')
        self.data_array[:,0,0,0] = visRcp
        self.data_array[:,0,0,1] = visLcp
        
        self.flag_array = NP.full((self.Nblts,1,self.Nfreqs,self.Npols),False,dtype='bool')
        self.nsample_array = NP.full((self.Nblts,1,self.Nfreqs,self.Npols),1,dtype='float')
        self.vis_units = 'uncalib'
        
        flags_ew_lcp = NP.where(srhFits.ewAntAmpLcp[frequency] == 1e6)[0]
        flags_ew_rcp = NP.where(srhFits.ewAntAmpRcp[frequency] == 1e6)[0]
        flags_ew = NP.unique(NP.append(flags_ew_lcp, flags_ew_rcp))
        flags_n_lcp = NP.where(srhFits.nAntAmpLcp[frequency] == 1e6)[0]
        flags_n_rcp = NP.where(srhFits.nAntAmpRcp[frequency] == 1e6)[0]
        flags_n = NP.unique(NP.append(flags_n_lcp, flags_n_rcp))
        
        flags_arr = NP.zeros((31,97), dtype = 'bool')
        flags_arr[flags_n,:] = True
        flags_arr[:,flags_ew] = True
        flags_arr = NP.reshape(flags_arr, (31*97))
        
        self.flag_array[:,0,0,0] = flags_arr
        self.flag_array[:,0,0,1] = flags_arr

        self.antenna_positions = NP.zeros((self.Nants_telescope,3))
        for ant in NP.arange(0, 97):
            self.antenna_positions[ant] = [0, (ant - 32) * 9.8, 0]
        for ant in NP.arange(97, 128):
            self.antenna_positions[ant] =  [-(ant - 96) * 9.8, 0, 0]
        self.baseline_array = 2048 * (self.ant_2_array + 1) + self.ant_1_array + 1 + 2**16
#        self.uvw_array = self.antenna_positions[self.ant_1_array] - self.antenna_positions[self.ant_2_array]
        self.uvw_array = NP.zeros((lastVisibility,3))

        lst = Time(scanTime,scale='utc', location=self.srhEarthLocation).sidereal_time('mean')
        hourAngle = lst.to('rad').value - self.phase_center_ra
        for vis in range(lastVisibility):
            self.uvw_array[vis] = self.baseline2uvw(hourAngle,coords.dec.rad,self.ant_2_array[vis], self.ant_1_array[vis])
        
        super().write_uvfits(dstName,write_lst=False,spoof_nonessential=True,run_check=False)
