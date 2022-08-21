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
from astropy.coordinates import EarthLocation
from astropy import units as u

class UVSrhUvFits(UVData):
    """
    """
    def __init__(self):
        self.srhEarthLocation = EarthLocation(lon=102.217 * u.deg, lat=51.759 * u.deg, height=799.0 * u.m)
        self.srhLocation = srhEarthLocation.T.to_value()
        self.phi_operator = NP.array([
            [-NP.sin(self.srhLat), 0., NP.cos(self.srhLat)],
            [0., 1., 0.],
            [NP.cos(self.srhLat), 0., NP.sin(self.srhLat)]
            ])
        
        def baseline2uvw(hourAngle, declination, antenna0, antenna1):
            if ((antenna0 >= 1 and antenna0 <= 128) and (antenna1 >= 129 and antenna1 <= 192)):
                base = np.array([antenna1 - 128.5,antenna0 - 64.5,0.])
            elif ((antenna1 >= 1 and antenna1 <= 128) and (antenna0 >= 129 and antenna0 <= 192)):
                base = np.array([antenna0 - 128.5,antenna1 - 64.5,0.])
            elif ((antenna0 >= 1 and antenna0 <= 128) and (antenna1 >= 1 and antenna1 <= 128)):
                base = np.array([0.,antenna0 - antenna1,0.])
            elif ((antenna0 >= 129 and antenna0 <= 192) and (antenna1 >= 129 and antenna1 <= 192)):
                base = np.array([antenna0 - antenna1,0.,0.])
            
            base *= 4.9;
            
            uvw_operator = NP.array([
                [ NP.sin(hourAngle),		 NP.cos(hourAngle),		0.	  ],
                [-NP.sin(declination)*NP.cos(hourAngle),  NP.sin(declination)*NP.sin(hourAngle), NP.cos(declination)], 
                [ NP.cos(declination)*NP.cos(hourAngle), -NP.cos(declination)*NP.sin(hourAngle), NP.sin(declination)]  
                ])
        
            return NP.dot(uvw_operator, NP.dot(self.phi_operator, base))

    def write_uvfits(self, srcName, dstName, lastVisibilities=8192, frequeny=0, scan=0, calibrating=True, uvSize=1024, averaging=0):
        srhFits = SrhFitsFile(srcName, uvSize)
        srhFits.setFrequencyChannel(frequency)
        if (calibrating):
            srhFits.calibrate(frequency)
        srhFits.vis2uv(scan, average=averaging)
        visRcp = srhFits.visRcp[frequency,scan,0:lastVisibility].copy()
        visLcp = srhFits.visLcp[frequency,scan,0:lastVisibility].copy()
        coords = COORD.get_sun(Time(srhFits.dateObs))
        for vis in range(lastVisibility):
            i = vis // 128
            j = vis % 128 
            visLcp[vis] *= NP.exp(1j*(-(srhFits.ewAntPhaLcp[freq, j]+srhFits.ewLcpPhaseCorrection[freq, j]) + (srhFits.sAntPhaLcp[freq, i] + srhFits.sLcpPhaseCorrection[freq, i])))
            visRcp[vis] *= NP.exp(1j*(-(srhFits.ewAntPhaRcp[freq, j]+srhFits.ewRcpPhaseCorrection[freq, j]) + (srhFits.sAntPhaRcp[freq, i] + srhFits.sRcpPhaseCorrection[freq, i])))
            visLcp[vis] /= (srhFits.ewAntAmpLcp[freq, j] * srhFits.sAntAmpLcp[freq, i])
            visRcp[vis] /= (srhFits.ewAntAmpRcp[freq, j] * srhFits.sAntAmpRcp[freq, i])
        
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
        self.freq_array[0] = 5.8e9
        self.history = 'SRH'
        self.instrument = 'SRH'
        self.integration_time = NP.full(self.Nblts,0.1)
        self.antenna_diameters = NP.full(self.Nants_telescope,2.)
        self.lst_array = NP.full(self.Nblts,0.)
        self.object_name = 'Sun'
        self.polarization_array = NP.array([-1,-2])
        self.spw_array = [1]
        self.telescope_location = srhLocation
        self.telescope_name = 'SRH'
        self.time_array = NP.full(self.Nblts,Time(srhFits.dateObs).jd)
        self.data_array = NP.zeros((self.Nblts,1,self.Nfreqs,self.Npols),dtype='complex')
        self.data_array[:,0,0,0] = visRcp
        self.data_array[:,0,0,1] = visLcp
            
        self.flag_array = NP.full((self.Nblts,1,self.Nfreqs,self.Npols),False,dtype='bool')
        self.nsample_array = NP.full((self.Nblts,1,self.Nfreqs,self.Npols),1,dtype='float')
        self.vis_units = 'uncalib'

        self.antenna_positions = NP.zeros((self.Nants_telescope,3))
        for ant in NP.arange(0, 128):
            self.antenna_positions[ant] = [0, (ant - 63.5) * 4.9, 0]
        for ant in NP.arange(128, 192):
            self.antenna_positions[ant] =  [(ant - 127.5) * 4.9, 0, 0]
        self.baseline_array = 2048 * (self.ant_2_array + 1) + self.ant_1_array + 1 + 2**16
        
        self.uvw_array = self.antenna_positions[self.ant_1_array] - self.antenna_positions[self.ant_2_array]
        obs_time1 = Time(Time(srhFits.dateObs),scale='utc', location=self.srhEarthLocation)
        lst = obs_time1.sidereal_time('mean')
        hA = lst.to('rad').value - self.phase_center_ra
        for vis in range(self.lastVisibility):
            self.uvw_array[vis] = self.baseline2uvw(hA,coords.dec.rad,self.ant_2_array[vis] + 1, self.ant_1_array[vis] + 1)
        
        self.write_uvfits(dstName,write_lst=False,spoof_nonessential=True)

#------------------------------------------------------------------------------
# UV = UVData()
# UV.ant_1_array = srhFits.antennaA[0:lastVisibility]
# UV.ant_2_array = srhFits.antennaB[0:lastVisibility]
# UV.antenna_names = srhFits.antennaNames
# UV.antenna_numbers = NP.array(list(map(int,srhFits.antennaNumbers)))
# UV.Nants_data = NP.union1d(srhFits.antennaA[0:lastVisibility],srhFits.antennaB[0:lastVisibility]).size
# UV.Nants_telescope = srhFits.antennaNames.shape[0]
# UV.Nfreqs = 1
# UV.Npols = 2
# UV.Ntimes = 1
# UV.Nspws = 1
# UV.Nbls = lastVisibility
# UV.Nblts = lastVisibility * UV.Ntimes
# UV.phase_type = 'phased'
# coords = COORD.get_sun(Time(srhFits.dateObs))
# UV.phase_center_ra = coords.ra.rad
# UV.phase_center_dec = coords.dec.rad
# UV.phase_center_app_ra = NP.full(UV.Nblts,coords.ra.rad)
# UV.phase_center_app_dec = NP.full(UV.Nblts,coords.dec.rad)
# UV.phase_center_epoch = 2000.0
# UV.channel_width = 1e7
# UV.freq_array = NP.zeros((1,1))
# UV.freq_array[0] = 5.8e9
# UV.history = 'SRH'
# UV.instrument = 'SRH'
# UV.integration_time = NP.full(UV.Nblts,0.1)
# UV.antenna_diameters = NP.full(UV.Nants_telescope,2.)
# UV.lst_array = NP.full(UV.Nblts,0.)
# UV.object_name = 'Sun'
# UV.polarization_array = NP.array([-1,-2])
# UV.spw_array = [1]
# UV.telescope_location = srhLocation
# UV.telescope_name = 'SRH'
# UV.time_array = NP.full(UV.Nblts,Time(srhFits.dateObs).jd)
# UV.data_array = NP.zeros((UV.Nblts,1,UV.Nfreqs,UV.Npols),dtype='complex')
# UV.data_array[:,0,0,0] = visRcp
# UV.data_array[:,0,0,1] = visLcp
    
# UV.flag_array = NP.full((UV.Nblts,1,UV.Nfreqs,UV.Npols),False,dtype='bool')
# UV.nsample_array = NP.full((UV.Nblts,1,UV.Nfreqs,UV.Npols),1,dtype='float')
# UV.vis_units = 'uncalib'

# phi = 0.903338787600965
# phi_operator = NP.array([
#     [-NP.sin(phi), 0., NP.cos(phi)],
#     [0., 1., 0.],
#     [NP.cos(phi), 0., NP.sin(phi)]
#     ])


# UV.antenna_positions = NP.zeros((UV.Nants_telescope,3))
# for ant in NP.arange(0, 128):
#     UV.antenna_positions[ant] = [0, (ant - 63.5) * 4.9, 0]

# for ant in NP.arange(128, 192):
#     UV.antenna_positions[ant] =  [(ant - 127.5) * 4.9, 0, 0]




# for ant in NP.arange(0, 128):
#     UV.antenna_positions[ant] = utils.ECEF_from_ENU(UV.antenna_positions[ant],ssrtLat,ssrtLon,ssrtHeight) - UV.telescope_location

# for ant in NP.arange(128, 192):
#     UV.antenna_positions[ant] =  utils.ECEF_from_ENU(UV.antenna_positions[ant],ssrtLat,ssrtLon,ssrtHeight) - UV.telescope_location





# UV.baseline_array = 2048 * (UV.ant_2_array + 1) + UV.ant_1_array + 1 + 2**16

# UV.uvw_array = UV.antenna_positions[UV.ant_1_array] - UV.antenna_positions[UV.ant_2_array]

# for vis in range(lastVisibility):
#     UV.uvw_array[vis] = NP.dot(phi_operator, UV.uvw_array[vis])

# uvw_array_phased = utils.phase_uvw(coords.ra.rad, coords.dec.rad,UV.uvw_array)

# PL.figure()
# PL.plot(UV.uvw_array.T[0],UV.uvw_array.T[1],'.')
# PL.plot(uvw_array_phased.T[0],uvw_array_phased.T[1],'.')

# location = EarthLocation(lon=102.217 * u.deg, lat=51.759 * u.deg, height=799.0 * u.m)
# obs_time1 = Time(Time(srhFits.dateObs),scale='utc', location=location)
# lst = obs_time1.sidereal_time('mean')
# hA = lst.to('rad').value - UV.phase_center_ra

# uvw_array = NP.zeros_like(UV.uvw_array)
# for vis in range(lastVisibility):
#     uvw_array[vis] = UVW.base2uvw(hA,coords.dec.rad,UV.ant_2_array[vis] + 1, UV.ant_1_array[vis] + 1)

# UV.uvw_array = uvw_array
# PL.plot(UV.uvw_array.T[0],UV.uvw_array.T[1],'.')
#------------------------------------------------------------------------------
#UV.write_uvfits('srh.uvfits',write_lst=False,spoof_nonessential=True)
