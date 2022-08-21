# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 10:45:18 2018

@author: Sergey
"""

from scipy import signal
import numpy as np
import pylab as pl

N = 6000
t = np.linspace(0,1,N)
sig = np.cos(500*np.pi*t**3)
fsig = signal.stft(sig,nperseg=250)
#fig = pl.figure()
#sp1 = fig.subplot(111)
#sp2 = fig.subplot(212)
#sp1.plot(sig)
#sp2.imshow(np.abs(fsig[2]))

pl.imshow(np.abs(fsig[2]))
pl.show()