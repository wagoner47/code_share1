#!/usr/bin/env python
# encoding: utf-8

import smtplib
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from datetime import datetime
import os, sys


def noticeEMail(starttime, usr, psw, fromaddr, toaddr, exit=0, nohup=None, plots=None):
	"""
	Sends an email message through GMail once the script is completed.  
	Developed to be used with AWS so that instances can be terminated 
	once a long job is done. Only works for those with GMail accounts.

	starttime : a datetime() object for when to start run time clock

	usr : the GMail username, as a string

	psw : the GMail password, as a string 

	fromaddr : the email address the message will be from, as a string

	toaddr : a email address, or a list of addresses, to send the 
			 message to
	"""
	
	
	# Calculate run time
	runtime=datetime.now() - starttime

	# Initialize SMTP server
	server=smtplib.SMTP('smtp.gmail.com:587')
	server.starttls()
	server.login(usr,psw)

	# Send email
	senddate=datetime.strftime(datetime.now(), '%Y-%m-%d')
	subject="Your job has completed"
	msg = MIMEMultipart()
	msg['From']=fromaddr
	msg['To']=toaddr
	msg['Date']=senddate
	msg['Subject']=subject
	# m="Date: %s\r\nFrom: %s\r\nTo: %s\r\nSubject: %s\r\nX-Mailer: My-Mail\r\n\r\n" % (senddate, fromaddr, toaddr, subject)
	text='''

	Exit status = {}
	Job runtime = {}\n\n'''.format(exit, runtime)
	msg.attach(MIMEText(text))
	if nohup is not None :
		part = MIMEApplication(open(nohup, 'rb').read(), Content_Disposition='attachement; filename="%s"' % os.path.basename(nohup), Name=os.path.basename(nohup))
		msg.attach(part)
	
	msg.attach(MIMEText('\n\n'))

	if plots is None :
		plots = []
	
	for p in plots :
		part = MIMEApplication(open(p, 'rb').read(), Content_Disposition='attachment; filename="%s"' % os.path.basename(p), Name=os.path.basename(p))
		msg.attach(part)
		msg.attach(MIMEText('\n\n'))

	server.sendmail(fromaddr, toaddr, msg.as_string())
	server.quit()


if __name__ == '__main__':
    # Start time of script
    starttime=datetime.now()

    # Some long job...we'll count to 100,000,000
    count=0
    for i in xrange(100000000):
        count+=1
        
    # Fill these in with the appropriate info...
    usr=''
    psw=''
    fromaddr=''
    toaddr=''

    # Send notification email
    noticeEMail(starttime, usr, psw, fromaddr, toaddr)
