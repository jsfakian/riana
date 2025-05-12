import os, shutil, tarfile
from datetime import datetime

from django.shortcuts import render, redirect
from django.contrib.auth import login, logout
from django.contrib.auth.models import User
from django.contrib.sites.shortcuts import get_current_site
from django.template.loader import render_to_string
from django.utils.http import urlsafe_base64_encode, urlsafe_base64_decode
from django.http import FileResponse
from django.utils.encoding import force_bytes, force_str
from django.core.mail import EmailMessage
from django.contrib.auth.views import PasswordResetView, PasswordResetDoneView, PasswordResetConfirmView, PasswordResetCompleteView
from django.contrib.auth.decorators import login_required
from .forms import EmailOrUsernameAuthenticationForm, RegisterForm, CalcForm, ContactForm
from .tokens import activation_token
from . import settings
from django.contrib.auth.views import LoginView
import matlab.engine

# Register
def register(request):
    if request.method == 'POST':
        form = RegisterForm(request.POST)
        if form.is_valid():
            user = form.save(commit=False)
            user.is_active = False
            user.set_password(form.cleaned_data['password1'])
            user.save()
            # send activation email
            current_site = get_current_site(request)
            mail_subject = 'Activate your account'
            message = render_to_string('activation_email.html', {
                'user': user,
                'domain': current_site.domain,
                'uid': urlsafe_base64_encode(force_bytes(user.pk)),
                'token': activation_token.make_token(user),
            })
            EmailMessage(mail_subject, message, to=[user.email]).send()
            return render(request, 'activation_complete.html')
    else:
        form = RegisterForm()
    return render(request, 'register.html', {'form': form})

# Activate
def activate(request, uidb64, token):
    try:
        uid = force_str(urlsafe_base64_decode(uidb64))
        user = User.objects.get(pk=uid)
    except Exception:
        user = None
    if user and activation_token.check_token(user, token):
        user.is_active = True
        user.save()
        login(request, user)
        return redirect('login')
    return render(request, 'activation_invalid.html')

# Login & logout
def home_login(request):
    if request.user.is_authenticated:
        return redirect("calc")
    return LoginView.as_view(
        template_name="home_login.html",
        authentication_form=EmailOrUsernameAuthenticationForm
    )(request)

def logout_view(request):
    logout(request)
    return redirect('login')

# Forgot Password uses built‑ins
class ForgotPassword(PasswordResetView):
    template_name = 'forgot_password.html'
    email_template_name = 'password_reset_email.html'
    success_url = '/password-reset-done/'

class FPWDone(PasswordResetDoneView):
    template_name = 'password_reset_done.html'

class FPRConfirm(PasswordResetConfirmView):
    template_name = 'password_reset_confirm.html'
    success_url = '/password-reset-complete/'

class FPWComplete(PasswordResetCompleteView):
    template_name = 'password_reset_complete.html'

@login_required
def manual(request):
    """
    Show the PDF user manual (embedded).
    """
    # Relative path to the manual PDF in media
    manual_path = settings.MEDIA_URL + 'manual_tsibidis_RIANA_final.pdf'
    # Build absolute URI including host and port
    manual_url = request.build_absolute_uri(manual_path)
    return render(request, 'manual.html', {
        'manual_url': manual_url
    })

@login_required
def manual_url(request):
    """
    Streams the PDF file itself.
    """
    manual_path = os.path.join(
        settings.MEDIA_ROOT,
        'manual_tsibidis_RIANA_final.pdf'
    )
    return FileResponse(
        open(manual_path, 'rb'),
        content_type='application/pdf'
    )

# MATLAB Calc (protected)
@login_required

def calc(request):
    # Parse all your inputs from request.POST or request.GET
    Ep1               = float(request.POST['fluence'])
    wavelength1       = float(request.POST['wavelength'])
    tp1               = float(request.POST['pulse_dur'])
    t_delay1          = float(request.POST['pulse_sep'])
    t_max1            = float(request.POST['max_time'])
    L1                = float(request.POST['thickness'])
    material          = request.POST['material']
    material_substrate= request.POST['substrate']
    n1, k1            = map(float, request.POST['n1']), map(float, request.POST['k1'])
    n2, k2            = map(float, request.POST['n2']), map(float, request.POST['k2'])

    user_id = request.user.id
    user_email = request.user.email
    username = request.user.username

    def background_job():
        # --- Call MATLAB function ---
        eng = matlab.engine.start_matlab()
        tf_dir = os.path.join(settings.BASE_DIR, 'Thin_films_ode15s')
        eng.cd(tf_dir, nargout=0)

        # no outputs, we just let calc write its files
        eng.calc(
            Ep1, wavelength1, tp1, t_delay1, t_max1, L1,
            material, material_substrate, n1, k1, n2, k2,
            nargout=0
        )
        eng.quit()

        # assume MATLAB wrote to BASE_DIR/<material>/
        src_folder = os.path.join(settings.BASE_DIR, material)

        # prepare user‐specific directory
        user_dir = os.path.join(settings.MEDIA_ROOT, 'simulations', str(user_id))
        os.makedirs(user_dir, exist_ok=True)

        # build a unique name
        ts     = datetime.now().strftime('%Y%m%d_%H%M%S')
        params = f"t{L1}_f{Ep1}_wl{wavelength1}_pd{tp1}_ps{t_delay1}_mt{t_max1}"
        new_name = f"{ts}_{material}_{params}"

        dest_folder = os.path.join(user_dir, new_name)
        shutil.move(src_folder, dest_folder)

        # --- Create tar.gz ---
        tar_name = new_name + '.tar.gz'
        tar_path = os.path.join(user_dir, tar_name)
        with tarfile.open(tar_path, 'w:gz') as tz:
            tz.add(dest_folder, arcname=new_name)

        # --- Build download URL ---
        download_url = request.build_absolute_uri(
            settings.MEDIA_URL + f"simulations/{user_id}/{tar_name}"
        )

        # --- Email the user ---
        subject = 'Your Riana Simulation Results'
        body = (
            f"Hello {username},\n\n"
            "Your simulation has completed with these parameters:\n"
            f"  • Material: {material}\n"
            f"  • Substrate: {material_substrate}\n"
            f"  • Thickness: {L1} nm\n"
            f"  • Fluence: {Ep1} J/cm²\n"
            f"  • Wavelength: {wavelength1} nm\n"
            f"  • Pulse Duration: {tp1} fs\n"
            f"  • Pulse Separation: {t_delay1} fs\n"
            f"  • Max Time: {t_max1} ps\n\n"
            f"Download your full results here:\n{download_url}\n\n"
            "Thank you for using Riana."
        )
        EmailMessage(subject, body, to=[user_email]).send()
    else:
        form = CalcForm()
    
    # kick off the thread (daemon=True won't block Django shutdown)
    t = threading.Thread(target=background_job, daemon=True)
    t.start()

    # immediately return a :202 Accepted" or similar
    return JsonResponse({'status': 'started'}, status=202)

def contact(request):
    if request.method == 'POST':
        form = ContactForm(request.POST)
        if form.is_valid():
            cd = form.cleaned_data

            # Build the message body
            body = (
                f"From: {cd['name']} <{cd['email']}>\n\n"
                f"{cd['message']}"
            )

            # Create the email
            email = EmailMessage(
                subject=f"[Contact] {cd['subject']}",
                body=body,
                from_email=settings.DEFAULT_FROM_EMAIL,
                to=[settings.ADMIN_EMAIL],
                reply_to=[cd['email']],                # user's address in Reply-To
            )
            # Override the "From:"" header so it reads as from the user
            email.extra_headers = {
                'From': f"{cd['name']} <{cd['email']}>"
            }

            email.send(fail_silently=False)
            return render(request, 'contact_sent.html', {'name': cd['name']})
    else:
        form = ContactForm()
    return render(request, 'contact.html', {'form': form})
