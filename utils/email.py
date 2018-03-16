import smtplib
import socket


def notify(msg, subj=None, machine_name=socket.gethostname()):
    if subj is None:
        subj = msg

    sender = 'david.schmidt@rwth-aachen.de'
    receivers = ['davidschmidt314@gmail.com']

    message = """From: %s
    To: David Schmidt <davidschmidt314@gmail.com>
    Subject: %s

    %s
    """ % (machine_name, subj, msg)

    try:
        smtpObj = smtplib.SMTP('localhost')
        smtpObj.sendmail(sender, receivers, message)
        print("Successfully sent email")
    except smtplib.SMTPException:
        print("Error: unable to send email")
