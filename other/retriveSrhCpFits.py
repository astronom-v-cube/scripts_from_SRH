# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 01:47:01 2016

@author: Sergey
"""
from optparse import OptionParser;
import urllib.request
import datetime

parser = OptionParser();
parser.add_option("-d", "--date", dest="inDate",   default="20200111");
(cl_options,cl_args) = parser.parse_args();
requestedDate = cl_options.inDate;



date = datetime.datetime(2020,2,11)
for i in range(1): 
    date.year
    requestedDate = ('%04d%02d%02d') % (date.year, date.month, date.day)
    srhCpName = 'srh_cp_' + requestedDate + '.fits'
    srhURL = 'http://archive.rao.istp.ac.ru/SRH/' + srhCpName
    urllib.request.urlretrieve(srhURL, srhCpName)
    print(srhCpName)
    date += datetime.timedelta(days=1)

