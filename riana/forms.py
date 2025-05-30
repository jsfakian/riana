from django import forms
from django.contrib.auth.models import User
from django.contrib.auth.forms import AuthenticationForm, PasswordResetForm
from django.contrib.auth.password_validation import validate_password
from django.contrib.auth import authenticate, get_user_model

class ContactForm(forms.Form):
    name    = forms.CharField(max_length=100, label="Your Name")
    email   = forms.EmailField(label="Your Email")
    subject = forms.CharField(max_length=150)
    message = forms.CharField(widget=forms.Textarea, label="Message")

class EmailOrUsernameAuthenticationForm(AuthenticationForm):
    username = forms.CharField(label="Username or Email", max_length=254)

    def clean(self):
        username_or_email = self.cleaned_data.get("username")
        password = self.cleaned_data.get("password")
        if username_or_email and password:
            UserModel = get_user_model()
            # Try to fetch by email; fall back to treating it as username
            try:
                user = UserModel.objects.get(email__iexact=username_or_email)
                username = user.get_username()
            except UserModel.DoesNotExist:
                username = username_or_email
            self.user_cache = authenticate(
                self.request, username=username, password=password
            )
            if self.user_cache is None:
                raise self.get_invalid_login_error()
            self.confirm_login_allowed(self.user_cache)
        return self.cleaned_data

class RegisterForm(forms.ModelForm):
    username = forms.CharField(label='Username')
    email    = forms.EmailField(label='Email')
    password1 = forms.CharField(
        label='Password',
        widget=forms.PasswordInput,
        help_text='Must be at least 8 characters, not too common, and not entirely numeric.'
    )
    password2 = forms.CharField(
        label='Confirm Password',
        widget=forms.PasswordInput
    )

    class Meta:
        model = User
        fields = ['username', 'email']

    def clean_password1(self):
        p1 = self.cleaned_data.get('password1')
        if p1:
            # this will raise a ValidationError if it fails any of your
            # AUTH_PASSWORD_VALIDATORS (length, common, numeric, etc.)
            validate_password(p1, self.instance)
        return p1

    def clean(self):
        cleaned = super().clean()
        p1, p2 = cleaned.get('password1'), cleaned.get('password2')
        if p1 and p2 and p1 != p2:
            self.add_error('password2', 'Passwords must match')
        return cleaned

class CalcForm(forms.Form):
    Material = forms.ChoiceField(label='Material', choices=[
        ('Cr', 'Chromium'),
        ('Ti', 'Titanium'),
        ('Ni', 'Nickel'),
        ('Au', 'Gold'),
        ('Ag', 'Silver'),
        ('Al', 'Aluminum'),
        ('Cu', 'Copper'),
        ('W', 'Tungsten'),
        ('Pt', 'Platinum'),
        ('Steel', 'Steel'),
        ('Mo', 'Molybdenum'),
    ])
    Substrate = forms.ChoiceField(label='Substrate', choices=[
        ('SiO2', 'SiO2'),
        ('Soda lime glass', 'Soda lime glass'),
        ('Si', 'Si'),
        ('Air', 'Air'),
    ])
    thickness = forms.FloatField(
        label='Thickness (nm)',
        min_value=10, max_value=1000,
        help_text='Enter a value between 10 and 1000'
    )
    fluence = forms.FloatField(
        label='Fluence (J/cm²)',
        min_value=0.001, max_value=2,
        help_text='Recommended value range between 0.001 and 2',
        widget=forms.NumberInput(attrs={
            'class': 'mt-1 block w-full bg-gray-100 border border-gray-300 rounded p-2',
            'min': '0.001',
            'step': '0.001',
        })
    )
    wavelength = forms.FloatField(
        label='Wavelength (nm)',
        min_value=248, max_value=15000,
        help_text='Enter a value between 248 and 15000'
    )
    pulse_dur = forms.FloatField(
        label='Pulse Duration (fs)',
        min_value=100, max_value=500000,
        help_text='Enter a value between 100 and 500000'
    )
    pulse_sep = forms.FloatField(
        label='Pulse Separation (fs)',
        min_value=0, max_value=500000,
        help_text='Enter a value between 0 and 500000'
    )
    max_time = forms.FloatField(
        label='Maximum Time (ps)',
        min_value=0,
        help_text='Max allowed end interval four times the Recommended Value',
    )
    # Additional optical constants, default to 1
    n1 = forms.FloatField(label='n1', initial=1)
    k1 = forms.FloatField(label='k1', initial=1)
    n2 = forms.FloatField(label='n2', initial=1)
    k2 = forms.FloatField(label='k2', initial=1)

