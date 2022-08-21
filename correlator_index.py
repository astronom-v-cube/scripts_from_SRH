#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 09:01:37 2020

@author: mariagloba
"""

h = 2 # horizontal antenna index
v = 1 # vertical antenna index

ind = (h//4 - v//4) * 16 + (2 * 52 - v//4 + 1) * (v//4) *16 / 2 + v%4 * 4 + h%4
print(h, v, ind)