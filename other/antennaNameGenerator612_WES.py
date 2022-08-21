#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 09:17:13 2020

@author: svlesovoi
"""
import numpy as NP
import scipy.constants

Cf = scipy.constants.c / 1.467

eSectionName = 'EastAntennaDescriptor'
wSectionName = 'WestAntennaDescriptor'
sSectionName = 'SouthAntennaDescriptor'
loFrequencyName = 'loFrequency'

westReceiver = [1,2,3,4]
eastReceiver = [5,6,7,8]
southReceiver = [9,10,11,12]
diameter = 2
spacing = 4900

print('Name Receiver Channel X Y Z FrontEndID FeedID\n')

print('[network]')
print('SyncDriverIP=10.1.x.x')
print('SyncDriverPort=56765')
print('correlatorIP=10.1.12.2')
print('correlatorPort=56765')
print('localOscillatorIP=10.1.1.x')
print('')

srh0612WestAntennaNumber = 64
srh0612EastAntennaNumber = 64
srh0612SouthAntennaNumber = 64
srh0612AntennaNumber = srh0612SouthAntennaNumber + srh0612WestAntennaNumber + srh0612EastAntennaNumber

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

westAntDelay = NP.zeros(srh0612WestAntennaNumber, dtype='int')
eastAntDelay = NP.zeros(srh0612EastAntennaNumber, dtype='int')
southAntDelay = NP.zeros(srh0612SouthAntennaNumber, dtype='int')
westFrontEndID = NP.zeros(srh0612WestAntennaNumber, dtype='int')
eastFrontEndID = NP.zeros(srh0612EastAntennaNumber, dtype='int')
southFrontEndID = NP.zeros(srh0612SouthAntennaNumber, dtype='int')
westFeedID = NP.zeros(srh0612WestAntennaNumber, dtype='int')
eastFeedID = NP.zeros(srh0612EastAntennaNumber, dtype='int')
southFeedID = NP.zeros(srh0612SouthAntennaNumber, dtype='int')

receivers = NP.array([1,2,3,4, 5,6,7,8, 9,10,11,12])

for i in range(srh0612SouthAntennaNumber):
    southAntDelay[i] = NP.ceil((optSouthAntDelay[i]/Cf)*1e12)
for i in range(srh0612WestAntennaNumber):
    westAntDelay[i] = NP.ceil((optWestAntDelay[i]/Cf)*1e12)
for i in range(srh0612EastAntennaNumber):
    eastAntDelay[i] = NP.ceil((optEastAntDelay[i]/Cf)*1e12)

maxs = NP.array([southAntDelay.max(), westAntDelay.max(), eastAntDelay.max()])
commonMaxCableLength = maxs.max()

southAntDelay = commonMaxCableLength - southAntDelay
westAntDelay  = commonMaxCableLength - westAntDelay
eastAntDelay  = commonMaxCableLength - eastAntDelay

for i in range(srh0612AntennaNumber):
    rec = receivers[i//16]
    chan = i%16 + 1
    if (i == 0):
        print('[%s]'%wSectionName)
    elif (i == srh0612WestAntennaNumber):
        print('\n[%s]'%eSectionName)
    elif (i == srh0612WestAntennaNumber + srh0612EastAntennaNumber):
        print('\n[%s]'%sSectionName)
        
    if (i < srh0612WestAntennaNumber):
        print(r'%d\%s="W20%02d, %d, %d, 0, %d, 0, %d, %d, %d, %d"'%(i+1,\
                                                                    wSectionName,\
                                                                    srh0612WestAntennaNumber - i,\
                                                                    rec, \
                                                                    chan,\
                                                                    -(srh0612WestAntennaNumber - i - .5)*spacing,\
                                                                    westFrontEndID[srh0612WestAntennaNumber - i - 1], \
                                                                    westFeedID[srh0612WestAntennaNumber - i - 1], \
                                                                    diameter, \
                                                                    westAntDelay[srh0612WestAntennaNumber - i - 1]))
    elif  (i < srh0612WestAntennaNumber + srh0612EastAntennaNumber):
        print(r'%d\%s="E20%02d, %d, %d, 0, %d, 0, %d, %d, %d, %d"'%(i-srh0612WestAntennaNumber + 1,\
                                                                              eSectionName,\
                                                                              i - srh0612WestAntennaNumber + 1, \
                                                                              rec, \
                                                                              chan, \
                                                                              (i - srh0612WestAntennaNumber + .5)*spacing,\
                                                                              eastFrontEndID[i - srh0612WestAntennaNumber], \
                                                                              eastFeedID[i - srh0612WestAntennaNumber], \
                                                                              diameter,\
                                                                              eastAntDelay[i - srh0612WestAntennaNumber]))
    else:
#        print(r'%d\%s="S20%02d, %d, %d, %d, 0, 0, %d, %d, %d, %d"'%(i-srh0612WestAntennaNumber-srh0612EastAntennaNumber,\
#                                                                              sSectionName,\
#                                                                              i - srh0612WestAntennaNumber - srh0612EastAntennaNumber + 1,\
#                                                                              rec,\
#                                                                              chan,\
#                                                                              -(i - srh0612WestAntennaNumber - srh0612EastAntennaNumber + .5)*spacing,\
#                                                                              southFrontEndID[i - srh0612WestAntennaNumber - srh0612EastAntennaNumber], \
#                                                                              southFeedID[i - srh0612WestAntennaNumber - srh0612EastAntennaNumber], \
#                                                                              diameter,\
#                                                                              southAntDelay[i - srh0612WestAntennaNumber - srh0612EastAntennaNumber]))
        print(r'{0}\{1}="S20{2:02d}, {3}, {4}, {5}, 0, 0, {6}, {7}, {8}, {9}"'.format(\
                                                                              i-srh0612WestAntennaNumber-srh0612EastAntennaNumber + 1,\
                                                                              sSectionName,\
                                                                              i - srh0612WestAntennaNumber - srh0612EastAntennaNumber + 1,\
                                                                              rec,\
                                                                              chan,\
                                                                              int(-(i - srh0612WestAntennaNumber - srh0612EastAntennaNumber + .5)*spacing),\
                                                                              southFrontEndID[i - srh0612WestAntennaNumber - srh0612EastAntennaNumber], \
                                                                              southFeedID[i - srh0612WestAntennaNumber - srh0612EastAntennaNumber], \
                                                                              diameter,\
                                                                              southAntDelay[i - srh0612WestAntennaNumber - srh0612EastAntennaNumber]))
