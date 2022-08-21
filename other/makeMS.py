#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 18:45:59 2021

@author: maria
"""

from srhFitsFile36 import SrhFitsFile
import pylab as PL
import srh36MS2

dateName = '20210405'
#fitsName = 'srh_20210329T024125'
#fitsName = 'srh_20210323T053828'
#fitsName = 'srh_20210323T030312'
#fitsName = 'srh_20210402T020941'
fitsName = 'srh_20210405T034842'

frequency = 0
scan = 4
fileName = '/home/sergeyvlesovoi/SRH36/' + dateName + '/' + fitsName + '.fit'

file = SrhFitsFile(fileName, 2049)
file.baselines = 22
file.useNonlinearApproach = True
file.setFrequencyChannel(frequency)
file.vis2uv(scan)
file.uv2lmImage()
PL.clf()
PL.imshow(file.lcp.real)

saveName = './' + fitsName + ('_%d'%frequency) + '.ms'
ms2Table = srh36MS2.SrhMs2(saveName)
ms2Table.initDataTable(file, file.frequencyChannel, phaseCorrect = True, amplitudeCorrect = False)
ms2Table.initAntennaTable(file)
ms2Table.initSpectralWindowTable(file, file.frequencyChannel)
ms2Table.initDataDescriptionTable()
ms2Table.initPolarizationTable()
ms2Table.initSourceTable()
ms2Table.initFieldTable()
ms2Table.initFeedTable(file, file.frequencyChannel)
ms2Table.initObservationTable(file)