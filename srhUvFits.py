#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 31 10:08:14 2021

@author: sergey_lesovoi
"""

import numpy as NP
from astropy.time import Time, TimeDelta
import base2uvw_612 as UVW
from astropy.io import fits

class srhUvFits():
    def __init__(self, fileName):
        self.fitsPath = fileName
        self.primaryHeader = fits.Header()
        self.primaryHeader['NAXIS'] = 2
        self.primaryHeader['NAXIS1'] = 777777701
        self.primaryHeader['NAXIS2'] = 0
        self.primaryHeader['EXTEND'] = ('T', 'All data in tebles')
        self.primaryHDU = fits.PrimaryHDU(header = self.primaryHeader)

    def flush(self):
        hduList = fits.HDUList(self.primaryHDU)
        hduList.writeto(self.fitsPath)
    
    
arr = NP.zeros(16)
col = fits.Column(name='C',format='D',array=arr)
tab = fits.BinTableHDU.from_columns([col])
hdu = fits.PrimaryHDU()
hdu.header.remove('NAXIS')
hdu.header.insert('EXTEND', ('NAXIS',2))
hdu.header.insert('EXTEND', ('NAXIS1',777777701))
hdu.header.insert('EXTEND', ('NAXIS2',0))
hdu.header['EXTEND'] = False
hdu.data = arr

print(hdu.header)
hduList = fits.HDUList([hdu])
hduList.writeto('qqqq.fits',output_verify='ignore',overwrite=True)
#hduList.writeto('qqqq.fits',overwrite=True)
print(hdu.header)

fd = fits.open('qqqq.fits',mode='update',output_verify='ignore')
fd[0].header.remove('NAXIS')
# fd[0].header.insert('EXTEND', ('NAXIS',2))
# fd[0].header.insert('EXTEND', ('NAXIS1',777777701))
# fd[0].header.insert('EXTEND', ('NAXIS2',0))
fd.writeto('qqqq1.fits',output_verify='ignore')
