#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 07:51:28 2019

@author: sergey
"""

import pymysql
import requests

db = pymysql.connect('10.1.1.69','lesovoi','wol6GhlN1', 'datacollector')
cursor = db.cursor()
#cursor.execute("select filename from files where bucket='SRH' and filename like 'mf_20190531%' order by filename desc limit 3;")
cursor.execute("select filename from files where filename like 'mf_20190530%' order by filename desc limit 30;")
fits = cursor.fetchall()
db.close()

req = requests.get('http://10.2.1.2:7480/SRH/%s' % fits[0])
downloadFits = open('%s' % fits[0], 'wb')
downloadFits.write(req.content)
downloadFits.close()

