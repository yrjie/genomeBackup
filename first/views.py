'''
All methods that provide JSON to Track components go here.

@author: Fabi
'''

from django.contrib.auth.decorators import login_required
from django.db import connections
from django.shortcuts import render_to_response
from django.template.context import RequestContext
import json, config
from django.http import HttpResponseForbidden, HttpResponse, HttpResponseRedirect
import os, datetime
from modules.uploader.models import FileForm, File
from django.views.decorators.csrf import csrf_exempt
import requests
import time
from browser.models import Track
import settings
from util.http import cookie, get_params
from django.db.models import Max


@login_required
def visit_home(request):
    return render_to_response('home.html', request)

def serve_static(request, path):
    response = HttpResponse(mimetype='text/{0}'.format(path.split('.')[-1]))
    with open(os.path.join('/home/browser/BASIC/basic2/modules/uploader/static', path)) as f:
        response.write(f.read())
    return response

def login_pgsadmin(request):
    if request.POST:
        form = LoginForm(request.POST)
        if form.is_valid():
            username = request.POST['username']
            password = request.POST['password']
            user = authenticate(username=username, password=password)
            if user is not None and user.is_active:
                login(request, user)
                return list_libraries(request)
            else:
                form._errors["password"] = ErrorList([u"Your username and password didn\'t match. Please try again."])
    else:
        form = LoginForm()
    context = RequestContext(request, {
        'form': form,
    })
    return render_to_response('user/account.login.html', context)

def logout_pgsadmin(request):
    logout(request)
    return visit_home(request)


def list_source(request):
    lst=get_params(request, 'menu')
    if lst=='all_list':
        context = RequestContext(request, {
            'subject_list': File.objects.order_by('-id')
        })
    else:
        context = RequestContext(request, {
            'subject_list': File.objects.filter(submitted_by=request.user).order_by('-id')
        })
    return render_to_response('subject.list.html', context)

def list_track(request):
    lst=get_params(request, 'menu')
    if lst=='all_list':
        context = RequestContext(request, {
            'subject_list': File.objects.order_by('-id')
        })
    else:
        context = RequestContext(request, {
            'subject_list': File.objects.filter(submitted_by=request.user).order_by('-id')
        })
    return render_to_response('subject.list.html', context)

def show_file(request):
    sub_id = request.REQUEST['id']
    #sub = utils.get_or_none(File, pk=sub_id)
    file = File.objects.filter(id=sub_id)[0]
    try:
        track= Track.objects.filter(name=file.track)[0]
    except:
        return HttpResponseRedirect("/basic2/uploader/list/")
    variables = RequestContext(request, {
        'track': track,
        'fileID': sub_id,
        'sourceID': file.source,
    })
    return render_to_response('subject.detail.html', variables)

@csrf_exempt
def update_track(request):
    if request.method == 'POST':
        file=File.objects.filter(id=request.POST['file_id'])[0]
        track= Track.objects.filter(id=request.POST['track_id'])[0]
        if 'save' in request.POST:
            track.shortLabel=request.POST['shortLabel']
            track.longLabel=request.POST['longLabel']
            track.save()
        elif 'delete' in request.POST:
            file.delete()
            track.delete()
            cmd='curl '+settings.BCS_URL+'/delete/'+request.POST['source_id']+'/'
            print cmd
            os.system(cmd)
    return HttpResponseRedirect("/basic2/uploader/list/")

@csrf_exempt
def submit_file(request):
    try:
        menu=get_params(request, 'menu')
        menu_file='file' in menu
    except:
        menu_file='files[]' in request.FILES
    form = FileForm(request.POST,username=request.user)
    if request.method == 'POST' and request.POST['type']!='' and form.is_valid():
        if len(request.POST['library'])<1:
            lib='user_'+str(request.user)
        else:
            lib=request.POST['library']
        asm=request.POST['asm']
        name=''
        if menu_file:
            name=str(request.FILES['files[]'])
            request.POST['input']='upload://files[]'
        else:
            url=request.POST['url']
            request.POST['type']='bigWig'
            name=url.split('/')[-1]
            request.POST['importmode']='remote'
            request.POST['input']=url
        shortLabel=request.POST['shortLabel']
        longLabel=request.POST['longLabel'] 
        fid=int(File.objects.all().aggregate(Max('id'))['id__max'])+1
        request.POST['tags']='{"assembly":"'+asm+'","project":"test","UCSC_type":'+'\"'+request.POST['type']+'\"}'
        request.POST['hook']='http://biogpu.d1.comp.nus.edu.sg:8001/reqbin/?asm='+asm+'&lib='+lib+'&track='+name.replace('.','_')+'&shortLabel='+shortLabel+'&longLabel='+longLabel+'&fid='+str(fid)
        request.POST['format']=request.POST['type']
        driver='mysql_table'
        if request.POST['type']=='bigWig':
            driver='bigWig'
        if request.POST['type']=='bigBed':
            driver='bigBed'
        if request.POST['type']=='bed':
            fileL=max(3,len(request.FILES['files[]'].readline().split('\t')))
            request.POST['format']+=str(fileL)
        if not menu_file:
            res=requests.post(config.BCS_URL+'/create/'+driver+'/', data=request.POST)
        else:
            res=requests.post(config.BCS_URL+'/create/'+driver+'/', data=request.POST, files=request.FILES)
        res_dict=json.loads(res.text)
        new_subject = form.save(commit=False)
        new_subject.status="Processing"
        if res_dict['status']=='OK':
            new_subject.name=name
            new_subject.type=request.POST['type']
            new_subject.source=res_dict['result']['source']['id']
            new_subject.track=name.replace('.','_')
            new_subject.asm=asm
            new_subject.library=lib
            new_subject.submitted_by = request.user
            new_subject.pub_date = datetime.datetime.now().replace(microsecond=0)
            new_subject.save()
        prop = [{'name':name, 'type':"", 'library': "", 'path': "", 'asm': "", 'status': res_dict['status']}]
        response_dict={"files":prop}        
        if 'submit' not in request.POST:
            return HttpResponse(json.dumps(response_dict), content_type='application/json')
        else:
            return HttpResponseRedirect("/basic2/uploader/list/")
        #return HttpResponse("{\"files\":[{\"path\": \"\", \"type\": \"\", \"name\": \""+str(name)+"\", \"library\": \""+"\", \"asm\": \""+"\", \"status\": \""+new_subject.status+"\"}]}", content_type='application/json')
    else:
        form = FileForm(username=request.user)
    if menu_file:
        form.fields.pop('url')
    variables = RequestContext(request, {
        'form': form,
        'menu_file': menu_file
    })
    return render_to_response('subject.submit.html', variables)

@csrf_exempt
def reqbin(request,term):
    print request.POST
    if request.method == 'POST' and request.POST['status']=='OK':
        temp=request.POST['source']
        beg=temp.find('id')+5
        end=temp.find(',',beg)
        id=temp[beg:end]
        cmd='/home/browser/BASIC/basic2-zzz/py /home/browser/BASIC/basic2-zzz/console/track_util.py bcs -l 10005 -n \'test_upload_'+id+'\' bed '+id
        ret=os.system(cmd)
        if ret==0:
            print 'uploaded successfully'
    return None

