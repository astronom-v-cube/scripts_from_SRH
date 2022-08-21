#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 16:29:03 2020

@author: mariagloba
"""

import ftplib
import datetime
import sys
import os
import fnmatch
from makeSrhImage import srhImage

if len(sys.argv) != 5:
    print('usage: makeImage [date: YYYYMMDD] [start time: HHMM] [stop time: HHMM] [number of frequency channel]')
    sys.exit()
    
date = sys.argv[1] 
timeStart = int(sys.argv[2])*100
timeStop = int(sys.argv[3])*100
freq = int(sys.argv[4])

path = 'SRH/'+date[:4]+'/'+date[4:6]+'/'+date[6:]+'/'
#path = 'SRH/2020/03/14/'

result = []
fd = ftplib.FTP('ftp.rao.istp.ac.ru', 'anonymous', 'anonymous')
filenames_list = fd.nlst(path)
filenames_list.sort()
for i in range(len(filenames_list)):
    basename = filenames_list[i].split('/')[-1].split('.')[0]
    timeNum = int(basename[-6:])
    if timeNum>timeStart and timeNum<timeStop:
        result.append(filenames_list[i])
        tempFile = open('temp.fits', 'wb')
        fd.retrbinary("RETR " + filenames_list[i] ,tempFile.write)
        tempFile.close()
        images = srhImage('./temp.fits', freq)
        images.makeImages()
