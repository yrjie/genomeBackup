'''
Created on Jul 27, 2012

@author: mulawadifh
'''
from os import path
from django.conf.urls import patterns

import pymongo

NAME = 'Uploader' # displayed on navbar

AUTH_GROUPS = ['uploader',]

DATABASE = {
    'ENGINE': 'django.db.backends.mysql',
    'NAME': 'plap08_basic2',
    'USER': 'browser',
    'PASSWORD': 'basic',
    'HOST': 'localhost',
    'PORT': '',
}

_mongoconn_pgs = pymongo.Connection(DATABASE['HOST'] or 'localhost')
MONGODB = _mongoconn_pgs[DATABASE['NAME']]

_whereami = path.dirname(path.abspath(__file__))
urlpatterns = patterns('',
    ('^static/css/(?P<path>.*)$', 'django.views.static.serve', { 
        'document_root': path.join(_whereami, 'static', 'css'),
    }),
    ('^static/js/(?P<path>.*)$', 'django.views.static.serve', { 
        'document_root': path.join(_whereami, 'static', 'js'),
    }),
    
    ('^static/bootstrap/css/(?P<path>.*)$', 'django.views.static.serve', { 
        'document_root': path.join(_whereami, 'static', 'bootstrap', 'css'),
    }),
    ('^static/bootstrap/js/(?P<path>.*)$', 'django.views.static.serve', { 
        'document_root': path.join(_whereami, 'static', 'bootstrap', 'js'),
    }),
)
urlpatterns += patterns('modules.uploader.views',
    ('^$', 'visit_home'),
    ('^list_source$', 'list_source'),
    ('^list_track$', 'list_track'),
    ('^submit$', 'submit_file'),
)