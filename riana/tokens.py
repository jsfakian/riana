from django.contrib.auth.tokens import PasswordResetTokenGenerator

class AccountActivationTokenGenerator(PasswordResetTokenGenerator):
    """
    Generates a one-time-use token for email activation.
    """
    def _make_hash_value(self, user, timestamp):
        # Combine user PK, timestamp, and active status
        return f"{user.pk}{timestamp}{user.is_active}"

# Instantiate the generator
activation_token = AccountActivationTokenGenerator()