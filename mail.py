import smtplib
from email.message import EmailMessage
def send_mail(to_email, subject, message, server='mail.signelix.com',
              from_email='<>'):

    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = from_email
    msg['To'] = ', '.join(to_email)
    msg.set_content(message)
    server = smtplib.SMTP(server, 26)
    server.set_debuglevel(1)
    pswd = input ('password')
    server.login('<>', pswd)  # user & password
    server.send_message(msg)
    server.quit()
    print('successfully sent the mail.')

send_mail(to_email=['<>'],
          subject='hello', message='test e-mail!')