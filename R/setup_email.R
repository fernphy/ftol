library(keyring)
library(blastula)

# Store SMTP credentials as a file
# with the filename "gmail_creds"
create_smtp_creds_file(
  file = "gmail_creds",
  user = "jnitta.no.reply@gmail.com",
  host = "smtp.gmail.com",
  port = 465,
  use_ssl = TRUE
)

# Send email
email <- compose_email("This is a test")

smtp_send(
    email,
    from = "jnitta.no.reply@gmail.com",
    to = "joelnitta@gmail.com",
    subject = "Testing the `smtp_send()` function",
    credentials = creds_file(file = "gmail_creds")
  )
