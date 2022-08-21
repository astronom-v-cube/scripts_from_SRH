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
nSectionName = 'NorthAntennaDescriptor'
nSectionName = 'NorthAntennaDescriptor'
loFrequencyName = 'loFrequency'

westReceiver = [1,2]
northReceiver = [3,4]
eastReceiver = [5,6,7,8]
centerReceiver = 9
diameter = 3

print('Name Receiver Channel X Y Z FrontEndID FeedID\n')

print('[%s]'%eSectionName)
Enames = 64
for i in range(Enames):
    if eastReceiver[i//16]:
        print(r'%d\%s="E30%02d, %d, %d, 0, %d, 0, 1909000, 1911000, %d"'%(i+1,eSectionName,i+1,eastReceiver[i//16],i%16+1,(i+1)*9800,diameter))
    else:
        print(r'%d\%s="E30%02d, %d, %d, 0, 0, 0, 1909000, 1911000, %d"'%(i+1,eSectionName,i+1,eastReceiver[i//16],0,diameter))
print('size=%d\n'%Enames)
    
print('[CenterAntennaDescriptor]')
if centerReceiver:
    print(r'%d\%s="C30, %d, %d, 0, 0, 0, 1909000, 1911000, %d"'%(1,'CenterAntennaDescriptor',centerReceiver,1,diameter))
else:
    print(r'%d\%s="C30, %d, %d, 0, 0, 0, 1909000, 1911000, %d"'%(1,'CenterAntennaDescriptor',centerReceiver,0,diameter))
print('size=1\n')

print('[%s]'%wSectionName)
Wnames = 32
for i in range(Wnames):
    if westReceiver[i//16]:
        print(r'%d\%s="W30%02d, %d, %d, 0, %d, 0, 1909000, 1911000, %d"'%(i+1,wSectionName,i+1,westReceiver[i//16],i%16+1,-(i+1)*9800,diameter))
    else:
        print(r'%d\%s="W30%02d, %d, %d, 0, 0, 0, 1909000, 1911000, %d"'%(i+1,wSectionName,i+1,westReceiver[i//16],0,diameter))
print('size=%d\n'%Wnames)

print('[%s]'%nSectionName)
Nnames = 32
for i in range(Nnames):
    if northReceiver[i//16]:
        print(r'%d\%s="N30%02d, %d, %d, %d, 0, 0, 1909000, 1911000, %d"'%(i+1,nSectionName,i+1,northReceiver[i//16],i%16+1,(i+1)*9800,diameter))
    else:
        print(r'%d\%s="N30%02d, %d, %d, 0, 0, 0, 1909000, 1911000, %d"'%(i+1,nSectionName,i+1,northReceiver[i//16],0,diameter))
print('size=%d'%Nnames)

srh0306NorthAntennaNumber = 31
srh0306WestAntennaNumber = 32
srh0306EastAntennaNumber = 64
srh0306AntennaNumber = srh0306NorthAntennaNumber + srh0306WestAntennaNumber + 1 + srh0306EastAntennaNumber

optCenterAntLength = 800.

optEastAntLength = [\
                   799.38, 799.59, 800.23, 800.44, 799.80, 800.12, 800.01, 799.91, 799.80, 799.59, 800.23, 800.44, 799.80, 800.02, 800.00, 800.01, \
                   800.33, 799.80, 800.01, 799.59, 800.76, 799.80, 800.01, 800.01, 799.80, 800.44, 799.59, 800.23, 800.01, 800.23, 800.01, 799.80, \
                   799.34, 799.77, 799.56, 799.34, 800.19, 799.98, 800.41, 800.41, 800.41, 799.66, 800.41, 799.98, 800.41, 800.19, 800.09, 800.41, \
                   800.62, 800.73, 799.34, 800.51, 800.62, 800.83, 799.56, 800.62, 800.62, 800.30, 800.41, 799.98, 800.41, 800.62, 800.51, 800.41]

optWestAntLength = [\
                   800.44, 800.65, 799.80, 800.12, 800.66, 799.80, 799.91, 799.91, 799.38, 800.33, 800.23, 799.80, 793.48, 799.80, 799.80, 799.80, \
                   800.65, 800.02, 800.01, 799.80, 799.59, 799.59, 800.01, 799.80, 799.48, 800.01, 799.59, 800.44, 800.01, 799.80, 800.01, 799.80]

optNorthAntLength = [\
                   800.66, 800.66, 800.45, 800.45, 800.45, 800.45, 800.66, 800.23, 800.02, 800.00, 800.45, 800.24, 800.45, 800.24, 800.24, 800.77, \
                   800.45, 800.87, 800.66, 801.94, 800.66, 800.65, 800.45, 800.66, 800.00, 800.87, 800.00, 800.66, 802.15, 800.88, 808.30, 802.57] 


westAntDelay = NP.zeros(srh0306WestAntennaNumber, dtype='int')
eastAntDelay = NP.zeros(srh0306EastAntennaNumber, dtype='int')
northAntDelay = NP.zeros(srh0306NorthAntennaNumber, dtype='int')

centerAntDelay = NP.ceil((optCenterAntLength[i]/Cf)*1e12)
for i in range(srh0306WestAntennaNumber):
    westAntDelay[i] = NP.ceil((optWestAntLength[i]/Cf)*1e12)
for i in range(srh0306EastAntennaNumber):
    eastAntDelay[i] = NP.ceil((optEastAntLength[i]/Cf)*1e12)
for i in range(srh0306NorthAntennaNumber):
    northAntDelay[i] = NP.ceil((optNorthAntLength[i]/Cf)*1e12)

maxs = NP.array([centerAntDelay, northAntDelay.max(), westAntDelay.max(), eastAntDelay.max()])
commonMaxCableLength = maxs.max()

northAntDelay = commonMaxCableLength - northAntDelay
westAntDelay  = commonMaxCableLength - westAntDelay
eastAntDelay  = commonMaxCableLength - eastAntDelay

#centerAntDelay = 7640.
#northAntDelay = NP.array([7572.,  9432.,  7078.,  8049.,  8400.,  6765.,  5798.,
#                          6796.,  2668.,  5593.,  3881.,  4493.,  3020.,  2241.,  3250.,
#                          3601.,  1286.,  1883.,  1947.,  6133.,  2313.,  1796.,  1343.,
#                          900.,   803.,  1459.,  1026.,   538.,  9637.,     0., 36697.])
#westAntDelay = NP.array([32714., 32550., 33959., 34356., 34450., 34761., 34236.,
#                         33600., 40501., 44335., 42747., 43969., 42838., 43711., 43109.,
#                         46946., 43980., 42890., 41852., 43138., 43475., 42824., 42679.,
#                         44180., 44150., 42495., 41634., 41543., 43724., 42257., 40674.,
#                         42335.])
#eastAntDelay = NP.array([11360., 12356., 16668., 43506., 45131., 43642., 45268.,
#                         18893., 14622., 11307., 17551., 42778., 16379., 16204., 39515.,
#                         44141., 16300., 45462., 43362., 44233., 45112., 46281., 44502.,
#                         45415., 45248., 44212., 41435., 42418., 40050., 41116., 41441.,
#                         39318., 41910., 41346., 42030., 19122., 41470., 41783., 41901.,
#                         40903., 42476., 42414., 43417., 41416., 43766., 41709., 41315.,
#                         41528., 40886., 41778., 37450., 39819., 41900., 40881., 43879.,
#                         37711., 43216., 41605., 40908., 43330., 46986., 41835., 41894.,
#                         42690.])

print('\nWCEN order---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')

for i in range(srh0306AntennaNumber):
    rec = i//16 + 1
    chan = i%16 + 1
    if (i == 0):
        print('[%s]'%wSectionName)
    elif (i == srh0306WestAntennaNumber):
        print('\n[CenterAntennaDescriptor]')
    elif (i == srh0306WestAntennaNumber + 1):
        print('\n[%s]'%eSectionName)
    elif (i == srh0306WestAntennaNumber + 1 + srh0306EastAntennaNumber):
        print('\n[%s]'%nSectionName)
        
    if (i < srh0306WestAntennaNumber):
        print(r'%d\%s="W30%02d, %d, %d, 0, %d, 0, 1909000, 1911000, %d, %d"'%(i+1,wSectionName,srh0306WestAntennaNumber - i, rec, chan, -(srh0306WestAntennaNumber - i)*9800,diameter,westAntDelay[srh0306WestAntennaNumber - 1 - i]))
    elif  (i == srh0306WestAntennaNumber):
        print(r'%d\%s="C30, %d, %d, 0, 0, 0, 1909000, 1911000, %d, %d"'%(1,'CenterAntennaDescriptor', rec, chan,diameter,centerAntDelay))
    elif  (i < srh0306WestAntennaNumber + srh0306EastAntennaNumber + 1):
        print(r'%d\%s="E30%02d, %d, %d, 0, %d, 0, 1909000, 1911000, %d, %d"'%(i-srh0306WestAntennaNumber,eSectionName,i - srh0306WestAntennaNumber, rec, chan, (i - srh0306WestAntennaNumber)*9800,diameter,eastAntDelay[i - srh0306WestAntennaNumber - 1]))
    else:
        print(r'%d\%s="N30%02d, %d, %d, %d, 0, 0, 1909000, 1911000, %d, %d"'%(i-srh0306WestAntennaNumber-srh0306EastAntennaNumber,nSectionName,i - srh0306WestAntennaNumber - srh0306EastAntennaNumber, rec, chan, (i - srh0306WestAntennaNumber - srh0306EastAntennaNumber)*9800,diameter,northAntDelay[i - srh0306WestAntennaNumber - srh0306EastAntennaNumber-1]))


dF = .1e9
f0 = 2.8e9
Nf = 30
scpiStr = 'LIST:FREQ '
for f in range(Nf):
        print(r'%d\%s=%2d'%(f+1,loFrequencyName,int((f0 + f*dF)*1e-3)))
        scpiStr += '%d MHz'%(int((f0 + f*dF)*1e-6))
        if f < Nf- 1:
            scpiStr += ','
        else:
            scpiStr += '\n'
            

