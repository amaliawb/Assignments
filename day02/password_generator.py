import secrets
import string

characters = string.ascii_letters + string.digits
password = "".join(secrets.choice(characters) for _ in range(4))

print(f"Your password is: {password}")
