# Convenience script to send an email after a workflow has finished.
import smtplib
import ssl
import sys

data = sys.stdin.read()

port = 465  # For SSL
smtp_server = "smtp.gmail.com"
sender_email = ""  # Enter your address
receiver_email = ""  # Enter receiver address
password = ""
message = f"""\
Subject: Snakemake Notification!

Hi there,
This message is sent from Python.
Your recent snakemake job either failed or finished running.
See the following log.

Have a nice day!

{data}

"""

context = ssl.create_default_context()
with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
    server.login(sender_email, password)
    server.sendmail(sender_email, receiver_email, message)
