from django.shortcuts import render

# Create your views here.

def do_MutationInfo(request):

    return render(request, 'MutationInfoApp/MutationInfoApp.html', {
#        'upload_form':upload_form, 
    })



