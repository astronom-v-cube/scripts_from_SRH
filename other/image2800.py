#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  1 09:41:54 2021

@author: mariagloba
"""

from srhFitsFile36 import SrhFitsFile
import pylab as PL
import srh36MS2
import numpy as NP
import os

#fileName = '/home/mariagloba/Work/fits/3-6/20210419/srh_20210419T050330.fit'
#fileName = 'SRH36_temp_20210501_4/srh_20210501T083414.fit'
#fileName = 'SRH36_temp_20210502_1/srh_20210502T030649.fit'
#fileName = 'SRH36_temp_20210502_2/srh_20210502T082025.fit'
#fileName = 'SRH36_temp_20210502_3/srh_20210502T084834.fit'
#fileName = 'SRH36_temp_20210429_1/srh_20210429T083051.fit'
#fileName = 'SRH36_temp_20210502_4/srh_20210429T020017.fit'
#fileName = 'SRH36_temp_20210503_1/srh_20210503T013622.fit'
#fileName = 'SRH36_temp_20210503_1/srh_20210503T013715.fit'
#fileName = 'SRH36_temp_20210503_1/srh_20210503T013808.fit'
#fileName = 'SRH36_temp_20210503_1/srh_20210503T013834.fit'
#fileName = 'SRH36_temp_20210503_1/srh_20210503T013900.fit'
#fileName = 'SRH36_temp_20210503_1/srh_20210503T013927.fit'
#fileName = 'SRH36_temp_20210503_1/srh_20210503T013953.fit'
#fileName = 'SRH36_temp_20210503_1/srh_20210503T014020.fit'
#fileName = 'SRH36_temp_20210509_1/srh_20210509T014748.fit'
#fileName = 'SRH36_temp_20210509_1/srh_20210509T014814.fit'
#fileName = 'SRH36_temp_20210509_1/srh_20210509T014841.fit'
#fileName = 'SRH36_temp_20210509_1/srh_20210509T014907.fit'
#fileName = 'SRH36_temp_20210509_1/srh_20210509T014933.fit'
#fileName = 'SRH36_temp_20210512_1/srh_20210512T055401.fit'
#fileName = 'SRH36_temp_20210512_1/srh_20210512T055454.fit'
#fileName = 'SRH36_temp_20210512_2/srh_20210512T014516.fit'
#fileName = 'SRH36_temp_20210513_1/srh_20210513T034149.fit'
#fileName = 'SRH36_temp_20210513_2/srh_20210513T072235.fit'
#fileName = 'SRH36_temp_20210513_2/srh_20210513T072209.fit'
#fileName = 'SRH36_temp_20210513_2/srh_20210513T072142.fit'
#fileName = 'SRH36_temp_20210513_2/srh_20210513T072116.fit'
#fileName = 'SRH36_temp_20210511_1/srh_20210511T024852.fit'
#fileName = 'SRH36_temp_20210510_1/srh_20210510T024749.fit'
#fileName = 'SRH36_temp_20210516_1/srh_20210516T031151.fit'
#fileName = 'SRH36_temp_20210515_1/srh_20210515T031232.fit'
#fileName = 'SRH36_temp_20210520_1/srh_20210520T030715.fit'
fileName = 'SRH36_temp_20210522_1/srh_20210522T015224.fit'

file = SrhFitsFile(fileName, 2049)
freq = 0
file.baselines = 5
file.useNonlinearApproach = True
file.getHourAngle(0)
file.solarPhase(freq)
file.updateAntennaPhase(freq, baselinesNumber = 5)
file.setFrequencyChannel(freq)
file.vis2uv(0, average = 20)
file.uv2lmImage()
#PL.clf()
#PL.imshow(file.lcp.real, vmax = 3)

saveName = fileName.split('.')[0] + '.ms'
ms2Table = srh36MS2.SrhMs2(saveName)
ms2Table.createMS(file, frequencyChannel = [0], phaseCorrect = True, amplitudeCorrect = True)

flags_ew = NP.where(file.ewAntAmpLcp[0] == 1e6)[0] + 1
flags_n = NP.where(file.nAntAmpLcp[0] == 1e6)[0]+98
flags = ','.join(map(str, NP.append(flags_ew, flags_n)))

#clean_script = 'casaclean.py'
#with open(clean_script, 'w') as f:
#    f.write('vis_name = ' + saveName)
#    f.write('flags = ' + flags)
#    f.write('tclean()')
#    f.write('\nimagename = ' + self.cleanParams['imagename = '][:-1] + '_dirty\'\n')
#    f.write('niter = 0\ntclean()')
#os.system('casa -c ' + clean_script.replace(' ', '\ '))

command = 'casa -c casaclean.py \'' + saveName + '\'  \'' + flags + '\''
os.system(command)
