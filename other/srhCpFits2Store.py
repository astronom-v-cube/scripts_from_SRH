#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
last mod 03.31.2021 by Pavel
@author: pavel
"""
import ftplib;
import datetime as DT;
  

dateName = DT.datetime.now().strftime("%Y%m%d");
datePName = DT.datetime.now().strftime("/%Y/%m/");

store = ftplib.FTP('10.1.1.35','sergeylesovoi','jill21ik');

fi = open('srh_cp_'+ dateName + '.fits','rb');
store.storbinary('STOR /mnt/FTP_DISK/SRH/corrPlot'+datePName+'srh_cp_'+ dateName + '.fits', fi);
fi.close();
store.close();


