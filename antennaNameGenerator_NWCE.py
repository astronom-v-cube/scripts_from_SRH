#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 09:17:13 2020

@author: svlesovoi
"""

eSectionName = 'EastAntennaDescriptor'
wSectionName = 'WestAntennaDescriptor'
nSectionName = 'NorthAntennaDescriptor'

westReceiver = [1,2]
northReceiver = [3,4]
eastReceiver = [5,6,7,8]
centerReceiver = 9

print('Name Receiver Channel X Y Z FrontEndID FeedID\n')

print('[%s]'%eSectionName)
Enames = 64
for i in range(Enames):
    if eastReceiver[i//16]:
        print(r'%d\%s="E30%02d, %d, %d, 0, %d, 0, 1909000, 1911000"'%(i+1,eSectionName,i+1,eastReceiver[i//16],i%16+1,(i+1)*9800))
    else:
        print(r'%d\%s="E30%02d, %d, %d, 0, 0, 0, 1909000, 1911000"'%(i+1,eSectionName,i+1,eastReceiver[i//16],0))
print('size=%d\n'%Enames)
    
print('[CenterAntennaDescriptor]')
if centerReceiver:
    print(r'%d\%s="C30, %d, %d, 0, 0, 0, 1909000, 1911000"'%(1,'CenterAntennaDescriptor',centerReceiver,1))
else:
    print(r'%d\%s="C30, %d, %d, 0, 0, 0, 1909000, 1911000"'%(1,'CenterAntennaDescriptor',centerReceiver,0))
print('size=1\n')

print('[%s]'%wSectionName)
Wnames = 32
for i in range(Wnames):
    if westReceiver[i//16]:
        print(r'%d\%s="W30%02d, %d, %d, 0, %d, 0, 1909000, 1911000"'%(i+1,wSectionName,i+1,westReceiver[i//16],i%16+1,-(i+1)*9800))
    else:
        print(r'%d\%s="W30%02d, %d, %d, 0, 0, 0, 1909000, 1911000"'%(i+1,wSectionName,i+1,westReceiver[i//16],0))
print('size=%d\n'%Wnames)

print('[%s]'%nSectionName)
Nnames = 32
for i in range(Nnames):
    if northReceiver[i//16]:
        print(r'%d\%s="N30%02d, %d, %d, %d, 0, 0, 1909000, 1911000"'%(i+1,nSectionName,i+1,northReceiver[i//16],i%16+1,(i+1)*9800))
    else:
        print(r'%d\%s="N30%02d, %d, %d, 0, 0, 0, 1909000, 1911000"'%(i+1,nSectionName,i+1,northReceiver[i//16],0))
print('size=%d'%Nnames)

srh0306NorthAntennaNumber = 31
srh0306WestAntennaNumber = 32
srh0306EastAntennaNumber = 64
srh0306AntennaNumber = srh0306NorthAntennaNumber + srh0306WestAntennaNumber + 1 + srh0306EastAntennaNumber

print('Another order')

for i in range(srh0306AntennaNumber):
    if (i == 0):
        print('[%s]'%wSectionName)
    elif (i == srh0306WestAntennaNumber):
        print('\n[CenterAntennaDescriptor]')
    elif (i == srh0306WestAntennaNumber + 1):
        print('\n[%s]'%eSectionName)
    elif (i == srh0306WestAntennaNumber + 1 + srh0306EastAntennaNumber):
        print('\n[%s]'%nSectionName)
        
    if (i < srh0306WestAntennaNumber):
        print(r'%d\%s="W30%02d, %d, %d, 0, %d, 0, 1909000, 1911000"'%(i+1,wSectionName,srh0306WestAntennaNumber - i,i//16,i%16+1,-(srh0306WestAntennaNumber - i)*9800))
    elif  (i == srh0306WestAntennaNumber):
        print(r'%d\%s="C30, %d, %d, 0, 0, 0, 1909000, 1911000"'%(i+1,'CenterAntennaDescriptor',i//16,i%16+1))
    elif  (i < srh0306WestAntennaNumber + srh0306EastAntennaNumber + 1):
        print(r'%d\%s="E30%02d, %d, %d, 0, %d, 0, 1909000, 1911000"'%(i+1,eSectionName,i - srh0306WestAntennaNumber,i//16,i%16+1,(i - srh0306WestAntennaNumber)*9800))
    else:
        print(r'%d\%s="N30%02d, %d, %d, %d, 0, 0, 1909000, 1911000"'%(i+1,nSectionName,i - srh0306WestAntennaNumber - srh0306EastAntennaNumber,i//16,i%16+1,(i - srh0306WestAntennaNumber - srh0306EastAntennaNumber)*9800))
