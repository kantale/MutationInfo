from django.conf.urls import patterns, url

from MutationInfoApp import views

urlpatterns = patterns('',
    url(r'', views.do_MutationInfo, name='do_MutationInfo'),
)


