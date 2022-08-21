#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 03:40:03 2019

@author: svlesovoi
"""

import pandas as PD
import pylab as PL
from matplotlib.ticker import MultipleLocator
import numpy as NP

def num_format(t, pos):
    if t > 0. and t < 49.:
#        return '%03d' % (int(t));
        return '%03d' % (ant_number[int(t)]);
    else:
        return('')


adc_delta = 0.5
ant_table = PD.ExcelFile('SRH_antenna_inputs_0-11.xlsx')
ant_sheet = ant_table.parse(0)

ant_number = ant_sheet['Unnamed: 1'][1:49]
ant_zero = ant_sheet['Unnamed: 2'][1:49]
ant_sun = ant_sheet['Unnamed: 3'][1:49]
ant_sky = ant_sheet['Unnamed: 4'][1:49]
ant_pow = ant_sheet['Unnamed: 5'][1:49]

fig1 = PL.figure(1,figsize=(20,12));
fig1.suptitle('SRH antenna inputs',fontsize='large');

sp1 = fig1.add_subplot(1,1,1);
sp1.xaxis.set_major_locator(MultipleLocator(1));
sp1.xaxis.set_major_formatter(PL.FuncFormatter(num_format));

sp1.plot(ant_sun / adc_delta)
sp1.plot(ant_sun / adc_delta * 10)

sp1.grid()
sp1.set_yscale('log', basey=2)

