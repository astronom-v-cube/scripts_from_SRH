#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 06:55:11 2021

@author: svlesovoi
"""

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
mail_content = '''Hello,world!'''

sender_address = 'svlesovoi@gmail.com'
sender_pass = 'Jill21_iK'
receiver_address = 'globa@iszf.irk.ru'

message = MIMEMultipart()
message['From'] = sender_address
message['To'] = receiver_address
message['Subject'] = 'A test mail sent by Python. It has an attachment.'   #The subject line
message.attach(MIMEText(mail_content, 'plain'))

session = smtplib.SMTP('smtp.gmail.com', 587) #use gmail with port
session.starttls() #enable security
session.login(sender_address, sender_pass) #login with mail_id and password
text = message.as_string()
session.sendmail(sender_address, receiver_address, text)
session.quit()
print('Mail Sent')
