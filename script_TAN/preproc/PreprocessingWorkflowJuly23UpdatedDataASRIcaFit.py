#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================
Get ica component and decide rejecting ones on Enrico data using MNE and - local computer version
===============================================
We firstly import subject list from sbatch 
we define session list (1, 2), state list (VD, FA, OP), reject dict, 
then we import eeglab format Raw data of one state during one session for single subject with MNE package. We apply:
1) make sure that there is event 254 and event 255, and crop the raw data between 254 and 255
2) a notch filter to remove powerline artifact (50 Hz)
3) a 1Hz-100Hz band-pass filter

Then concatenate the data of the same session with annotation engineering, detail utils please turn to 
utils_preProcessingWorkflowJuly05.py
3) ASR and ICA fitting:

====> output = subj+session_ica fif file that save ica mixing matrix for one session for one subject
Note: 
1. exception : subject 36, some subject can have several 254,255 events
------ please refer to excel Enrico recording summary.xlsx
2. events code: state + condition + session 
    1. state: 1:VD 2:FA 3:OP
    2. condition: 1:baseline 2:safe 3:threat
    3. session: 1:session1 2:session2
3. we fix sampling rate at 512 = for those file whose sfreq = 2048, we do a downsampling

Suggestions:
1) decide infomation storage format
2) 

Updated on July 2019

@author: Gansheng TAN aegean0045@outlook.com    based on Manu's codes
"""

##############################################################  Set-up header ###########################################
import mne
import importlib
import numpy as np
import numpy.matlib
from mne.report import Report
from autoreject import AutoReject
from autoreject import compute_thresholds
from autoreject import get_rejection_threshold 
import matplotlib
import matplotlib.pyplot as plt  # noqa
import matplotlib.patches as patches  # noqa
from autoreject import set_matplotlib_defaults  # noqa
from utils_ASR import *
from utils_PreprocessingWorkflowJuly23UpdatedData import *
from scipy.linalg import toeplitz
from scipy import signal
import sys
import encodings
import os

matplotlib.use('Agg')
mne.set_log_level('WARNING')

##################### OS path in INSERM computer #####################################################################
# raw_data_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/raw_data/'
# montage_fname = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/raw_data/Biosemi64_MAS_EOG.locs'
# # report_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/report/'
# full_epochs_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/'
#
##################### OS path in cluster ######################################################################
raw_data_path = '/mnt/data/gansheng/raw_data/'
montage_fname = '/mnt/data/gansheng/raw_data/Biosemi64_MAS_EOG.locs'
preProc_ica_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/preProc_ica/'
report_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/report/'
# full_epochs_path = '/mnt/data/gansheng/preClean_data/'


########################################## Algorithme parameter ############################################
cutoff = 10
pca_n_comp = 0.98
decim = 2

########################################## Initialization parameter##########################################
subj = sys.argv[1]
# subj = '25'
session_list=['1','2']
#state list defines the concatenating order
# state_list = ['VD','FA','OP']


################################ step00: cut and filter data and concatenate 3 recording in one session ############

###### set up montage
montage_biosemi=mne.channels.read_montage(montage_fname)

###### preproc for each raw file

psd_figs=[]
psd_captions=[]
ASR_figs=[]
ASR_captions=[]

session2conctn_list=[]
############### single subject report ###########################
rep = Report(image_format = 'png', subject = 'subj0'+subj)

for session in session_list:

    reject_state=[]
    conctn_list = []
    conctn_anno_list=[]
    
    epochs_ASR_clean,psd_figs2conctn,psd_captions2conctn,ASR_figs2conctn,ASR_captions2conctn = get_epochs_ASR_clean(subj,session)
    psd_figs=psd_figs+psd_figs2conctn
    psd_captions = psd_captions+psd_captions2conctn
    ASR_figs = ASR_figs + ASR_figs2conctn
    ASR_captions = ASR_captions+ASR_captions2conctn
    ############### step02 ICA components check ##########################
    if epochs_ASR_clean == False:
        continue
    else:
        ica = mne.preprocessing.ICA(n_components=pca_n_comp, method='fastica', random_state=11, max_iter=100)
        ica.fit(epochs_ASR_clean,decim=decim)
        preProc_ica_fname = preProc_ica_path+'subj0'+subj+'session'+session+'preProc_ica.fif'
        ica.save(preProc_ica_fname)

rep.add_figs_to_section(figs=psd_figs, captions=psd_captions, section = 'preprocessed full epochs plot',
                            image_format = 'svg')
rep.add_figs_to_section(figs=ASR_figs, captions=ASR_captions, section = 'preprocessed full epochs plot',
                            image_format = 'svg')
        
f_report = report_path + 'subj0 ' + subj+' asr report.html'
rep.save(f_report, open_browser = False, overwrite = True)
            
        

