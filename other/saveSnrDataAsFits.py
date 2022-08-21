#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 06:20:03 2022

@author: sergey_lesovoi
"""

from astropy.io import fits

freqs0306 = phaseEdit.srhFits.freqList
time0306 = phaseEdit.srhFits.freqTime
lcpAmps0306 = phaseEdit.srhFits.ampLcp
rcpAmps0306 = phaseEdit.srhFits.ampRcp
lcpRVisNorth0306 = phaseEdit.srhFits.visLcp[:,:,3007:3007+30]
rcpRVisNorth0306 = phaseEdit.srhFits.visLcp[:,:,3007:3007+30]
lcpRVisWestEast0306 = phaseEdit.srhFits.visLcp[:,:,3472:3472+96]
rcpRVisWestEast0306 = phaseEdit.srhFits.visLcp[:,:,3472:3472+96]

freqs0612 = phaseEdit_0612.srhFits.freqList
time0612 = phaseEdit_0612.srhFits.freqTime
lcpAmps0612 = phaseEdit_0612.srhFits.ampLcp
rcpAmps0612 = phaseEdit_0612.srhFits.ampRcp
lcpRVisSouth0612 = phaseEdit_0612.srhFits.visLcp[:,:,8192:8192+63]
rcpRVisSouth0612 = phaseEdit_0612.srhFits.visRcp[:,:,8192:8192+63]
lcpRVisWestEast0612 = phaseEdit_0612.srhFits.visLcp[:,:,10208:10208+127]
rcpRVisWestEast0612 = phaseEdit_0612.srhFits.visRcp[:,:,10208:10208+127]

freqs1224 = phaseEdit1224.srhFits.freqList
time1224 = phaseEdit1224.srhFits.freqTime
lcpAmps1224 = phaseEdit1224.srhFits.ampLcp
rcpAmps1224 = phaseEdit1224.srhFits.ampRcp
lcpRVisSouth1224 = phaseEdit1224.srhFits.visLcp[:,:,9452:9452+67]
rcpRVisSouth1224 = phaseEdit1224.srhFits.visLcp[:,:,9452:9452+67]
lcpRVisWestEast1224 = phaseEdit1224.srhFits.visLcp[:,:,11730:11730+138]
rcpRVisWestEast1224 = phaseEdit1224.srhFits.visRcp[:,:,11730:11730+138]

dataFormat = str(time0306.shape[1]) + 'D'
freqsColumn0306 = fits.Column(name='freqs0306',format='D',array=freqs0306)
freqsColumn0612 = fits.Column(name='freqs0612',format='D',array=freqs0612)
freqsColumn1224 = fits.Column(name='freqs1224',format='D',array=freqs1224)

timeColumn0306 = fits.Column(name='time0306',format=dataFormat,array=time0306)
timeColumn0612 = fits.Column(name='time0612',format=dataFormat,array=time0612)
timeColumn1224 = fits.Column(name='time1224',format=dataFormat,array=time1224)

dataFormat = str(lcpAmps0306.shape[1]*lcpAmps0306.shape[2]) + 'D'
lcpColumn0306 = fits.Column(name='lcpAmp0306',format=dataFormat,array=lcpAmps0306)
rcpColumn0306 = fits.Column(name='rcpAmp0306',format=dataFormat,array=rcpAmps0306)

dataFormat = str(lcpRVisNorth0306.shape[1]*lcpRVisNorth0306.shape[2]) + 'D'
lcpRVisNorthColumn0306 = fits.Column(name='lcpRVisNorth0306',format=dataFormat,array=lcpRVisNorth0306)
rcpRVisNorthColumn0306 = fits.Column(name='rcpRVisNorth0306',format=dataFormat,array=rcpRVisNorth0306)
dataFormat = str(lcpRVisWestEast0306.shape[1]*lcpRVisWestEast0306.shape[2]) + 'D'
lcpRVisWestEastColumn0306 = fits.Column(name='lcpRVisNorth0306',format=dataFormat,array=lcpRVisWestEast0306)
rcpRVisWestEastColumn0306 = fits.Column(name='rcpRVisNorth0306',format=dataFormat,array=rcpRVisWestEast0306)

dataFormat = str(lcpAmps0612.shape[1]*lcpAmps0612.shape[2]) + 'D'
lcpColumn0612 = fits.Column(name='lcpAmp0612',format=dataFormat,array=lcpAmps0612)
rcpColumn0612 = fits.Column(name='rcpAmp0612',format=dataFormat,array=rcpAmps0612)

dataFormat = str(lcpRVisSouth0612.shape[1]*lcpRVisSouth0612.shape[2]) + 'D'
lcpRVisSouthColumn0612 = fits.Column(name='lcpRVisSouth0612',format=dataFormat,array=lcpRVisSouth0612)
rcpRVisSouthColumn0612 = fits.Column(name='rcpRVisSouth0612',format=dataFormat,array=rcpRVisSouth0612)
dataFormat = str(lcpRVisWestEast0612.shape[1]*lcpRVisWestEast0612.shape[2]) + 'D'
lcpRVisWestEastColumn0612 = fits.Column(name='lcpRVisWestEast0612',format=dataFormat,array=lcpRVisWestEast0612)
rcpRVisWestEastColumn0612 = fits.Column(name='rcpRVisWestEast0612',format=dataFormat,array=rcpRVisWestEast0612)

dataFormat = str(lcpAmps1224.shape[1]*lcpAmps1224.shape[2]) + 'D'
lcpColumn1224 = fits.Column(name='lcpAmp1224',format=dataFormat,array=lcpAmps1224)
rcpColumn1224 = fits.Column(name='rcpAmp1224',format=dataFormat,array=rcpAmps1224)

dataFormat = str(lcpRVisSouth1224.shape[1]*lcpRVisSouth1224.shape[2]) + 'D'
lcpRVisSouthColumn1224 = fits.Column(name='lcpRVisSouth1224',format=dataFormat,array=lcpRVisSouth1224)
rcpRVisSouthColumn1224 = fits.Column(name='rcpRVisSouth1224',format=dataFormat,array=rcpRVisSouth1224)
dataFormat = str(lcpRVisWestEast1224.shape[1]*lcpRVisWestEast1224.shape[2]) + 'D'
lcpRVisWestEastColumn1224 = fits.Column(name='lcpRVisWestEast1224',format=dataFormat,array=lcpRVisWestEast1224)
rcpRVisWestEastColumn1224 = fits.Column(name='rcpRVisWestEast1224',format=dataFormat,array=rcpRVisWestEast1224)

fTableHdu = fits.BinTableHDU.from_columns([freqsColumn0306, freqsColumn0612, freqsColumn1224]);
tTableHdu = fits.BinTableHDU.from_columns([timeColumn0306, timeColumn0612, timeColumn1224]);

dTableHdu0306 = fits.BinTableHDU.from_columns([lcpColumn0306, rcpColumn0306])
dTableHduNorthVis0306 = fits.BinTableHDU.from_columns([lcpRVisNorthColumn0306, rcpRVisNorthColumn0306])
dTableHduWestEastVis0306 = fits.BinTableHDU.from_columns([lcpRVisWestEastColumn0306, rcpRVisWestEastColumn0306])

dTableHdu0612 = fits.BinTableHDU.from_columns([lcpColumn0612, rcpColumn0612])
dTableHduSouthVis0612 = fits.BinTableHDU.from_columns([lcpRVisSouthColumn0612, rcpRVisSouthColumn0612])
dTableHduWestEastVis0612 = fits.BinTableHDU.from_columns([lcpRVisWestEastColumn0612, rcpRVisWestEastColumn0612])

dTableHdu1224 = fits.BinTableHDU.from_columns([lcpColumn1224, rcpColumn1224])
dTableHduSouthVis1224 = fits.BinTableHDU.from_columns([lcpRVisSouthColumn1224, rcpRVisSouthColumn1224])
dTableHduWestEastVis1224 = fits.BinTableHDU.from_columns([lcpRVisWestEastColumn1224, rcpRVisWestEastColumn1224])

pHeader = fits.Header()
pHeader['DATE-OBS']     = '2022-04-25'
pHeader['TIME-OBS']     = '03:40:00'
pHeader['INSTRUME']     = 'SRH'
pHeader['ORIGIN']       = 'ISTP'
pHeader['OBS-LAT']      = '51.759'
pHeader['OBS-LONG']     = '102.217'
pHeader['OBS-ALT']      = '799'
pHeader['ANT_0306']      = '128'
pHeader['ANT_0612']      = '192'
pHeader['ANT_1224']      = '207'

pHdu = fits.PrimaryHDU(header=pHeader)
hduList = fits.HDUList([pHdu, fTableHdu, tTableHdu, dTableHdu0306, dTableHdu0612, dTableHdu1224])
hduList.writeto('srh_amp_snr_20220425' + '.fits',clobber=True)
hduList.close()

pHeader = fits.Header()
pHeader['DATE-OBS']     = '2022-04-25'
pHeader['TIME-OBS']     = '03:40:00'
pHeader['INSTRUME']     = 'SRH'
pHeader['ORIGIN']       = 'ISTP'
pHeader['OBS-LAT']      = '51.759'
pHeader['OBS-LONG']     = '102.217'
pHeader['OBS-ALT']      = '799'
pHeader['VISN0306']      = '30'
pHeader['VISW0306']      = '96'
pHeader['VISS0612']      = '63'
pHeader['VISW0612']      = '127'
pHeader['VISS1224']      = '67'
pHeader['VISW1224']      = '138'

pHdu = fits.PrimaryHDU(header=pHeader)
hduList = fits.HDUList([pHdu, fTableHdu, tTableHdu, 
                        dTableHduNorthVis0306, dTableHduWestEastVis0306, 
                        dTableHduSouthVis0612, dTableHduWestEastVis0612,
                        dTableHduSouthVis1224, dTableHduWestEastVis1224])

hduList.writeto('srh_vis_snr_20220425' + '.fits',clobber=True)
hduList.close()
