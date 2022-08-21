#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 09:26:44 2019

@author: svlesovoi & pbdedyukhin
"""

import os, fnmatch;
import ftplib;
import datetime;
  

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

today = datetime.datetime.today();
fitNames = findFits('/home/serg/SRH/data/heliograph48/','mf_' + today.strftime('%Y%m%d') + '_*.fit')

fitNames.sort();
fitNames = fitNames[1:];
"""
today = datetime.datetime.today();
"""

store = ftplib.FTP('84.237.21.39','ftpwriter','cn5HeuA4u')
for fitName in fitNames:
    fi = open(fitName,'rb')
    store.storbinary('STOR /SRH' + today.strftime('/%Y/%m/%d/') + fitName.split('/')[-1], fi)
    fi.close()

store.close()
