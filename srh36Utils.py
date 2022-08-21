#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 01:10:02 2021

@author: svlesovoi
"""
import os, fnmatch

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

def number2name():
    srh36AntNumber2AntName = {}
    for ant in range(32):
        srh36AntNumber2AntName[ant + 1] = 'W%d'%(32 - ant)
    srh36AntNumber2AntName[33] = 'C0'
    for ant in range(64):
        srh36AntNumber2AntName[ant + 34] = 'E%d'%(ant + 1)
    for ant in range(32):
        srh36AntNumber2AntName[ant + 98] = 'N%d'%(ant + 1)
    return srh36AntNumber2AntName
