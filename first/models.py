from django.db import models
from django.contrib.auth.models import User
from django.forms import ModelForm, Textarea, TextInput
from django import forms
from django.forms.widgets import CheckboxSelectMultiple
from django.shortcuts import render_to_response
import datetime

class File(models.Model):
    def __unicode__(self):
        return self.subject_id

    subject_id = models.CharField(max_length=50, unique=True)
    gender = models.CharField(max_length=10, choices=(('M', 'Male'), ('F', 'Female'), ('U', 'Unknown')))
    ethnicity = models.CharField(max_length=50, blank=True)
    relatives = models.CharField(max_length=50, blank=True)
    human_pedigree = models.CharField(max_length=50, blank=True)
    serum_monitoring = models.CharField(max_length=50, blank=True)
    treatment = models.CharField(max_length=50, blank=True)
    status = models.CharField('Subject status', max_length=1024, default='', blank=True, null=True)
    pub_date = models.DateTimeField('Date submitted')
    submitted_by = models.ForeignKey(User, null=True, blank=True)
    
class FileForm(ModelForm):
    subject_id = forms.CharField(widget=forms.TextInput({ "placeholder": "Any unique name" }), help_text='Any unique name', error_messages={'required': 'Please enter a unique subject name'})
    gender = forms.TypedChoiceField(choices=(('U', 'Unknown'), ('M', 'Male'), ('F', 'Female')), empty_value='Gender cannot be empty')
    ethnicity = forms.CharField(widget=forms.TextInput({ "placeholder": "Chinese, Singapore..." }), required=False)
    class Meta:
        model = File
        exclude = ('status','submitted_by','pub_date',)
        widgets = {
            'relatives': Textarea(attrs={'cols': 40, 'rows': 5}),
            'treatment': Textarea(attrs={'cols': 40, 'rows': 5}),
        }