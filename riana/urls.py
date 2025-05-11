"""
URL configuration for riana project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path
from django.conf.urls.static import static
from rest_framework_simplejwt.views import TokenObtainPairView, TokenRefreshView
from . import views, settings

urlpatterns = [
    path('admin/', admin.site.urls),
    path('api/token/', TokenObtainPairView.as_view(), name='token_obtain_pair'),
    path('api/token/refresh/', TokenRefreshView.as_view(), name='token_refresh'),

    # Web UI routes
    path('', views.home_login, name='home'),
    path('register/', views.register, name='register'),
    path('activate/<uidb64>/<token>/', views.activate, name='activate'),
    path('login/', views.home_login, name='login'),
    path('logout/', views.logout_view, name='logout'),
    path('contact/', views.contact, name='contact'),
    path('manual/', views.manual, name='manual'),
    path('manual/url/', views.manual_url,  name='manual_url'),  # streams the PDF

    # Password reset routes
    path('forgot-password/', views.ForgotPassword.as_view(), name='forgot_password'),
    path('password-reset-done/', views.FPWDone.as_view(), name='password_reset_done'),
    path('reset/<uidb64>/<token>/', views.FPRConfirm.as_view(), name='password_reset_confirm'),
    path('password-reset-complete/', views.FPWComplete.as_view(), name='password_reset_complete'),

    # Protected MATLAB calc page
    path('calc/', views.calc, name='calc'),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
