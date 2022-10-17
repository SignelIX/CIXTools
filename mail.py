import smtplib
from os.path import basename
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.utils import  formatdate

from email.message import EmailMessage
def send_mail(to_email, subject, message, server='mail.signelix.com',
              from_email='eric@signelix.com'):

    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = from_email
    msg['To'] = ', '.join(to_email)
    msg.set_content(message)
    server = smtplib.SMTP(server, 26)
    server.set_debuglevel(1)
    pswd = input ('password')
    server.login('eric@signelix.com', pswd)  # user & password
    server.send_message(msg)
    server.quit()
    print('successfully sent the mail.')


def send_mail_ext(to_email, subject, message, server='mail.signelix.com',
              from_email='eric@signelix.com', files=None,):
    assert isinstance(to_email, list)

    msg = MIMEMultipart()
    msg['From'] = from_email
    msg['To'] = ', '.join(to_email)
    msg['Date'] = formatdate(localtime=True)
    msg['Subject'] = subject

    msg.attach(MIMEText(message))

    for f in files or []:
        with open(f, "rb") as fil:
            part = MIMEApplication(
                fil.read(),
                Name=basename(f)
            )
        # After the file is closed
        part['Content-Disposition'] = 'attachment; filename="%s"' % basename(f)
        msg.attach(part)


    smtp = smtplib.SMTP(server)
    smtp.set_debuglevel(1)
    pswd = input ('password')
    smtp.login('eric@signelix.com', pswd)
    smtp.sendmail(from_email, to_email, msg.as_string())
    smtp.close()

if __name__== "__main__":
    send_mail_ext(to_email=['eric@citadeldiscovery.io'],
          subject='hello', message='test e-mail!', files = ['/Users/eric/testcomponent/LICENSE'])