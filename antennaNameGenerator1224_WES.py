#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 09:17:13 2020

@author: svlesovoi
"""
import numpy as NP
import scipy.constants
import base2uvw_1224

Cf = scipy.constants.c / 1.467

cSectionName = 'CenterAntennaDescriptor'
eSectionName = 'EastAntennaDescriptor'
wSectionName = 'WestAntennaDescriptor'
sSectionName = 'SouthAntennaDescriptor'
loFrequencyName = 'loFrequency'

diameter = 1
spacing = 2450

#print('Name Receiver Channel X Y Z FrontEndID FeedID\n')

srh1224WestAntennaNumber = 69
srh1224EastAntennaNumber = 69
srh1224SouthAntennaNumber = 68
srh1224CenterAntennaNumber = 1

srh1224AntennaNumber = srh1224WestAntennaNumber + srh1224CenterAntennaNumber + srh1224EastAntennaNumber + srh1224SouthAntennaNumber

optEastAntDelay = [\
350.75, 350.64, 350.43, 350.86, 350.55, 350.96, 350.96, 351.17, 351.17, 351.08, 350.74, 350.96, 350.98, 350.96, 350.74, 350.96, \
351.17, 351.3,  351,    350.66, 351,    350.86, 350.74, 350.96, 350.85, 351.17, 350.74, 350.74, 350.34, 351.01, 350.77, 351.19, \
351.01, 350.98, 350.17, 350.77, 350.98, 351.01, 350.98, 351.19, 350.98, 351.01, 350.76, 350.77, 350.54, 350.84, 351.06, 350.87, \
345.45, 350.77, 350.84, 350.96, 350.21, 350.81, 351.02, 350.17, 351.17, 350.96, 350.74, 350.96, 350.11, 351.38, 351.62, 356 ]

optWestAntDelay = [\
349.57, 351.19, 351.4,  351.5,  351.29, 350.65, 348.85, 351.29, 350.79, 351.36, 350.55, 351.3,  351.25, 351.86, 350.03, 350.76, \
352.26, 351.57, 351.4,  351.68, 351.23, 351.22, 351.72, 351.72, 351.72, 351.82, 351.5,  351.72, 351.71, 351.72, 351.86, 351.43, \
351.5,  351.82, 351.08, 351.93, 351.23, 351.96, 351.01, 351.5,  351.07, 351.72, 351.1,  351.07, 350.86, 350.78, 350.7,  351, \
348.65, 350.98, 350.4,  351.07, 351.1,  351.01, 350.86, 350.77, 351.02, 350.97, 351.45, 351.65, 351.4,  351.06, 351.01, 356.65 ]

optSouthAntDelay = [\
351.06, 350.42, 350.63, 351.38, 351.27, 351.27, 351.06, 351.27, 350.42, 351.7,  351.87, 350.81, 351.02, 350.59, 350.81, 350.59, \
350.81, 349.49, 351.06, 350.63, 349.68, 350.63, 351.7,  350.84, 351.27, 350.27, 350.3,  351.08, 352.57, 350.34, 350.06, 350.06, \
350.06, 350.06, 350.07, 350.88, 351.07, 350.64, 349.7,  350.32, 350.96, 350.11, 351.81, 350.74, 350.32, 350.11, 350.11, 351.2, \
351.17, 350.53, 350.53, 350.8,  351.65, 351.22, 351.22, 351.01, 351.2,  351.01, 350.79, 350.58, 351.22, 349.94, 352.23, 356.11 ]

westAntDelay = NP.zeros(srh1224WestAntennaNumber, dtype='int')
eastAntDelay = NP.zeros(srh1224EastAntennaNumber, dtype='int')
southAntDelay = NP.zeros(srh1224SouthAntennaNumber, dtype='int')
westFrontEndID = NP.zeros(srh1224WestAntennaNumber, dtype='int')
eastFrontEndID = NP.zeros(srh1224EastAntennaNumber, dtype='int')
southFrontEndID = NP.zeros(srh1224SouthAntennaNumber, dtype='int')
westFeedID = NP.zeros(srh1224WestAntennaNumber, dtype='int')
eastFeedID = NP.zeros(srh1224EastAntennaNumber, dtype='int')
southFeedID = NP.zeros(srh1224SouthAntennaNumber, dtype='int')

receivers = NP.array([1,2,3,4, 5,6,7,8, 9,10,11,12,13])

# for i in range(srh1224SouthAntennaNumber):
#     southAntDelay[i] = NP.ceil((optSouthAntDelay[i]/Cf)*1e12)
# for i in range(srh1224WestAntennaNumber):
#     westAntDelay[i] = NP.ceil((optWestAntDelay[i]/Cf)*1e12)
# for i in range(srh1224EastAntennaNumber):
#     eastAntDelay[i] = NP.ceil((optEastAntDelay[i]/Cf)*1e12)

maxs = NP.array([southAntDelay.max(), westAntDelay.max(), eastAntDelay.max()])
commonMaxCableLength = maxs.max()

southAntDelay = commonMaxCableLength - southAntDelay
westAntDelay  = commonMaxCableLength - westAntDelay
eastAntDelay  = commonMaxCableLength - eastAntDelay

for i in range(srh1224AntennaNumber):
    rec = receivers[i//16]
    chan = i%16 + 1
    if (i == 0):
        print('[%s]'%wSectionName)
    elif (i == srh1224WestAntennaNumber):
        print('\n[%s]'%cSectionName)
    elif (i == srh1224WestAntennaNumber + 1):
        print('\n[%s]'%eSectionName)
    elif (i == srh1224WestAntennaNumber + srh1224EastAntennaNumber + 1):
        print('\n[%s]'%sSectionName)
        
    if (i < srh1224WestAntennaNumber):
        print(r'%d\%s="W10%02d, %d, %d, 0, %d, 0, %d, %d, %d, %d"'%(i+1,\
                                                                    wSectionName,\
                                                                    srh1224WestAntennaNumber - i,\
                                                                    rec, \
                                                                    chan,\
                                                                    base2uvw_1224.distFromCenter(i + 1)*spacing,\
                                                                    westFrontEndID[srh1224WestAntennaNumber - i - 1], \
                                                                    westFeedID[srh1224WestAntennaNumber - i - 1], \
                                                                    diameter, \
                                                                    westAntDelay[srh1224WestAntennaNumber - i - 1]))
    elif  (i < srh1224WestAntennaNumber + 1):
        print(r'%d\%s="C10%02d, %d, %d, 0, %d, 0, %d, %d, %d, %d"'%(i - srh1224WestAntennaNumber + 1,\
                                                                              cSectionName,\
                                                                              i - srh1224WestAntennaNumber + 1, \
                                                                              rec, \
                                                                              chan, \
                                                                              base2uvw_1224.distFromCenter(i + 1)*spacing,\
                                                                              eastFrontEndID[i - srh1224WestAntennaNumber], \
                                                                              eastFeedID[i - srh1224WestAntennaNumber], \
                                                                              diameter,\
                                                                              eastAntDelay[i - srh1224WestAntennaNumber]))
    elif  (i < srh1224WestAntennaNumber + srh1224EastAntennaNumber + 1):
        print(r'%d\%s="E10%02d, %d, %d, 0, %d, 0, %d, %d, %d, %d"'%(i-srh1224WestAntennaNumber,\
                                                                              eSectionName,\
                                                                              i - srh1224WestAntennaNumber, \
                                                                              rec, \
                                                                              chan, \
                                                                              base2uvw_1224.distFromCenter(i + 1)*spacing,\
                                                                              eastFrontEndID[i - srh1224WestAntennaNumber - 1], \
                                                                              eastFeedID[i - srh1224WestAntennaNumber - 1], \
                                                                              diameter,\
                                                                              eastAntDelay[i - srh1224WestAntennaNumber - 1]))
    else:
        print(r'{0}\{1}="S10{2:02d}, {3}, {4}, {5}, 0, 0, {6}, {7}, {8}, {9}"'.format(\
                                                                              i-srh1224WestAntennaNumber-srh1224EastAntennaNumber,\
                                                                              sSectionName,\
                                                                              i - srh1224WestAntennaNumber - srh1224EastAntennaNumber,\
                                                                              rec,\
                                                                              chan,\
                                                                              -base2uvw_1224.distFromCenter(i + 1)*spacing,\
                                                                              southFrontEndID[i - srh1224WestAntennaNumber - srh1224EastAntennaNumber - 1], \
                                                                              southFeedID[i - srh1224WestAntennaNumber - srh1224EastAntennaNumber - 1], \
                                                                              diameter,\
                                                                              southAntDelay[i - srh1224WestAntennaNumber - srh1224EastAntennaNumber - 1]))
    
print('\n[array]')
print('antennaNumber={0}'.format(srh1224AntennaNumber))
print('westAntennaNumber={0}'.format(srh1224WestAntennaNumber))
print('eastAntennaNumber={0}'.format(srh1224EastAntennaNumber))
print('southAntennaNumber={0}'.format(srh1224SouthAntennaNumber))
print('centerAntennaNumber={0}'.format(1))

print('\n[loFrequencies]')
print('1\loFrequency={0}'.format(12000000))
print('2\loFrequency={0}'.format(12800000))
print('3\loFrequency={0}'.format(13600000))
print('4\loFrequency={0}'.format(14400000))
print('5\loFrequency={0}'.format(15200000))
print('6\loFrequency={0}'.format(16000000))
print('7\loFrequency={0}'.format(16800000))
print('8\loFrequency={0}'.format(17600000))
print('9\loFrequency={0}'.format(18400000))
print('10\loFrequency={0}'.format(19200000))
print('11\loFrequency={0}'.format(20000000))
print('12\loFrequency={0}'.format(20800000))
print('13\loFrequency={0}'.format(21600000))
print('14\loFrequency={0}'.format(22400000))
print('15\loFrequency={0}'.format(23200000))
print('16\loFrequency={0}'.format(24000000))
print('size=16')

print('\n[network]')
print('SyncDriverIP=10.1.10.3')
print('SyncDriverPort=56765')
print('correlatorIP=10.1.11.5')
print('correlatorPort=56765')
print('localOscillatorIP=10.1.11.6')

print('\n[receiver]')
print('autoStart=true')
print('dataDelay=10000')
print('dataDuration=2000000')
print('delayTracking=true')
print('fringeStopping=true')
print('internalSync=false')
print('localOscillatorStartStop=false')
print('oneBitCorrelation=0')

print('\n[SyncDriver]')
print('frequencyPulseDuration=100000')
print('frequencySetTime=100000')
print('frequencySwitchTime=200000')
print('frequencybandMask=7')
print('leftPolarizationDuration=10000000')
print('rightPolarizationDuration=10000000')
print('polarizationCyclesPerFrequency=1')

print('\n[FITS]')
print('fitsPath=/../../SRH1224/currentData')
print('frequencyListSize=16')
print('fullPacketsInFits=320')
