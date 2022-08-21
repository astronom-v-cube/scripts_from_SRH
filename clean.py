#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 14:02:30 2021

@author: maria
"""

#fitsName = 'srh_20210329T024125'
#fitsName = 'srh_20210323T053828'
#fitsName = 'srh_20210402T020941_0'
fitsName = 'srh_20210405T034842_0'
visName = '/home/sergeyvlesovoi/Python Scripts/' + fitsName + '.ms'
flagdata(vis = visName, antenna = '23, 45, 57, 114, 119, 126')
#flagdata(vis = '/home/sergeyvlesovoi/Python Scripts/srh_.ms', antenna = '23, 45, 57, 114, 119, 126')

tclean(vis=visName,
      imagename='/home/sergeyvlesovoi/Pictures/images/' + fitsName,
      cell=3.0,
      niter=10000,
      gain = 0.1,
      threshold = '0.1mJy',
      datacolumn = 'corrected',
      imsize=[1024,1024],
      stokes = 'RR',
      deconvolver='multiscale',
      weighting='uniform')
