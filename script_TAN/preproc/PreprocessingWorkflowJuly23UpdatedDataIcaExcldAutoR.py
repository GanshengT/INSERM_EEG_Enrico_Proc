### cluster version
#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================
Preprocessing on Enrico data using MNE and ASR - cluster version
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
    we replicate ASR and load ica file
        
4) Autoreject and concatenate two sessions
====> output = full_epoch fif file that save the the full recording for one subject
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
# #
##################### OS path in cluster ######################################################################
raw_data_path = '/mnt/data/gansheng/raw_data/'
montage_fname = '/mnt/data/gansheng/raw_data/Biosemi64_MAS_EOG.locs'
report_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/report/'
preProc_ica_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/preProc_ica/'
full_epochs_path = '/mnt/data/gansheng/preClean_data/'


########################################## Algorithme parameter ############################################
cutoff = 4
pca_n_comp = 0.98
decim = 2

########################################## Initialization parameter##########################################
subj_list = [sys.argv[1]]
# subj_list = ['02']
session_list=['1','2']
#state list defines the concatenating order
# state_list = ['VD','FA','OP']
state_list = ['VD','FA','OP']
power_freq_array = [50]
reject_raw_data_session1 = {'74':['FA','OP','VD'],'62':['FA','OP','VD'],'75':['FA','OP','VD']}
reject_raw_data_session2 = {'74':['FA','OP','VD'],'62':['FA','OP','VD'],'75':['FA','OP','VD']}

# bad channel rejection is not apllied in the preproc, bad channels will be defined by eyes later
bad_channels={'02':{'1':['P2','FC5'],'2':['P2','FC5']},
              '04':{'2':['FC6']},
              '07':{'1':['Iz'],'2':['F8','T7','TP8']},
              '10':{'1':['F8','AF8','F6','Fp1','AF7','FC6','FT8'],'2':['Fp1','Fpz','Fp2','AF8','F8']},
              '11':{'2':['T8','F4']},
              '12':{'1':['P2']},
              '14':{'1':['P2'],'2':['P2']},
              '16':{'1':['Fp1','AF7','Fp2','AF8'],'2':['Fp1','AF7','Fp2','AF8','F6']},
              '19':{'1':['P2','T7'],'2':['P2']},
              '21':{'2':['Iz']},
              '22':{'1':['T7','TP7','Iz','P1','POz','Pz'],'2':['T7','TP7','Iz']},
              '25':{'1':['T8','Fp1'],'2':['FC1','C3','PO4','F2','Pz']},
              '26':{'1':['T8'],'2':['CPz']},
              '28':{'2':['CP2','PO7','Oz','POz']},
              '29':{'1':['F3','F5','F7','AF8','F4','F6','F8'],'2':['F6','T7','F3','F5','F7','AF8','F8']},
              '32':{'1':['P2'],'2':['P2','TP7']},
              '34':{'1':['P2'],'2':['P2']},
              '35':{'1':['T7','T8'],'2':['T8','PO8']},
              '36':{'1':['P2','PO4'],'2':['P2','PO4']},
              '37':{'1':['Iz']},
              '38':{'2':['C5','FC5','TP8']},
              '39':{'1':['P2','F8','AF8','Fp1','AF7'],'2':['P2','FT8','AF8','T8','P10']},
              '40':{'1':['P2','TP7'],'2':['P2','TP7']},
              '42':{'1':['P2'],'2':['P2']},
              '50':{'1':['T7']},
              '51':{'1':['P2'],'2':['P2']},
              '53':{'1':['T7'],'2':['T7']},
              '54':{'1':['P2'],'2':['P2']},
              '55':{'1':['Iz']},
              '56':{'1':['P2','T7'],'2':['P2','TP7']},
              '57':{'1':['P2'],'2':['P2']},
              '58':{'1':['P2','T8'],'2':['PO4']},
              '59':{'1':['P2','PO4']},
              '60':{'1':['P2'],'2':['P2']},
              '63':{'2':['PO8']},
              '64':{'1':['C1']},
              '65':{'1':['P2'],'2':['P2']},
              '67':{'1':['FT8']},
              '68':{'1':['P2'],'2':['P2']},
              '70':{'1':['PO4','O2','FC3','FC5','F4','F6'],'2':['PO4','O2','FC5','FC3']},
              '71':{'1':['P2','Iz'],'2':['P2','Iz','C1','Cz']},
              '73':{'2':['FCz','FT8']},
            #'75':{'1':['C6','P2','FT8','AF8','CP1','P9','PO4','O2']},
              '76':{'2':['T7']},
              '77':{'1':['F6'],'2':['O1','Oz','F6','O2']},
              '78':{'1':['P2'],'2':['P2']},
              '79':{'1':['P2','POz'],'2':['P2','POz','T7','Fp1','AF7']},
              '81':{'1':['Iz','Oz','Pz','CPz','PO4','P2','POz'],'2':['Iz','Oz','POz','CPz','P2','PO4','FC1','C1','Pz']},
              '82':{'1':['P2'],'2':['AFz']},
              '83':{'1':['T7'],'2':['T7']},
              '87':{'2':['P2']},
              '88':{'1':['FC2','T8'],'1':['F4','P8','CP4']},
              '90':{'1':['T7','P2'],'2':['P2']},
              '91':{'1':['P2'],'2':['P2']},
              '93':{'1':['FC5','Fp1','F3','PO4'],'2':['Fp1','F3','FC3','PO4']},
              '94':{'1':['Fp1','F6','AF8','Fp2','T7','T8'],'2':['Fp1','F6','AF8','Fp2','T7']},
              '95':{'1':['P2'],'2':['P2']}
             }


compnts2exclude_ica =     {'02':{'1':[0,17,18,20,22],'2':[0,17,11,15,23]},
                           '04':{'1':[0,1,6,12,13,14,17],'2':[0,19,16,14,15]},
                           '07':{'1':[1,4,0,20,31,33,36],'2':[28,26,1]},
                           '10':{'1':[0,4],'2':[0]},
                           '11':{'1':[0,11,12,13],'2':[0,12,13]},
                           '12':{'1':[0,1,3,19,20,21,22,23],'2':[0,4,2,19,14,15,16,22,23]},
                           '14':{'1':[2,13,17,8,16,18,20,31,30],'2':[0,1,7,10,12,25,26,28,22]},
                           '16':{'1':[1,3,6,5,9,14,15,22,23],'2':[1,7,12,14]},
                           '18':{'1':[1,7,3,14,17,18,19,22,23,26,24],'2':[0,4,9,18,19,20,25,27,28]},
                           '19':{'1':[0,1,9,11,12,23,21],'2':[9,8,16,20,28]},
                           '21':{'1':[0,8,13,15,17,11],'2':[0,10,18,19,20,24]},
                           '22':{'1':[0,6,8,11],'2':[1,7,9,12,18,14,19,15,16,17]},
                           '25':{'1':[0,9,23,15,16,17],'2':[0,18,17,12,21]},
                           '26':{'1':[0,5,14,15,17],'2':[6,8,11,12,21]},
                           '28':{'1':[0,10,12,16,17,24,23,32],'2':[0,6,13,15,21,23]},
                           '29':{'1':[1,4,8,13],'2':[0,8,20]},
                           '30':{'1':[0,6,12,13,18,19,21,22,24],'2':[0,10,13,14,18,19,20,21,23,24,25,29]},
                           '32':{'1':[14,12,13,20,21,24,25],'2':[1,6,7,19,18,17,21]},
                           '34':{'1':[0,7,11,17,15,21,22,32],'2':[0,4,11,16,24,30,31,35]},
                           '35':{'1':[0,2,18,23,34,37,38,41],'2':[10,17,43]},
                           '36':{'1':[0,4,9,16,20,25,23,27],'2':[0,3,4,8,9,15,17,19,21]},
                           '37':{'1':[0,27,22,31,41,42],'2':[0,7,21,22,38,33]},
                           '38':{'1':[5,4,8,11,18,17,16,20,26,27,29],'2':[0,2,5,19,22,27,36,35,34]},
                           '39':{'1':[0,18,21],'2':[2,27]},
                           '40':{'1':[0,5,4,13,15,12,29,25,23,30],'2':[0,9,14,10,16,19,29,23,24,20,32]},
                           '42':{'1':[0,5,10,18,14,15,21,24],'2':[0,10,15,17]},
                           '50':{'1':[1,3,5,7,16,10,19,30],'2':[3,4,5,6,7,13]},
                           '51':{'1':[0,3,4,5,7,8,10,11,12,14,16],'2':[1,13,11,16]},
                           '52':{'1':[0,12,13,22],'2':[1,31,34]},
                           '53':{'1':[5,6,13,16,21,23,20,22,33,34],'2':[5,8,14,18,19,22,23,21]},
                           '54':{'1':[0,20,23,25,31,34],'2':[0,15,25,30,31]},
                           '55':{'1':[0,10,16,14,15,18],'2':[2,8,7,11,14,13,15,16]},
                           '56':{'1':[1,10,24,25],'2':[3,6,15,14,21,28]},
                           '57':{'1':[],'2':[18]},
                           '58':{'1':[0,10,14,22,23],'2':[0,16,21,26,27]},
                           '59':{'1':[0,1,19,12,28,35],'2':[1,10,13,17,21,20,28,32,33]},
                           '60':{'1':[2,6,9,11,17,19,20,27,34],'2':[5,6,8,17,10,11,19,20,28,29]},
                           '63':{'1':[4,14,19,15,20,29],'2':[6,11,16,21]},
                           '64':{'1':[0,17,11,19],'2':[1]},
                           '65':{'1':[1,9,15,16,17,18,14],'2':[2,3,11,12,17,20,24,26,28,29,39]},
                           '67':{'1':[5,9,16,20],'2':[0,5,34,26,38]},
                           '68':{'1':[13,14,15,17,1],'2':[0,12,14,15,17]},
                           '69':{'1':[3,15,16,14],'2':[8,12,13]},
                           '70':{'1':[5,4,6,9,11,0,15,20,25,31,30],'2':[0,8,17,25]},
                           '71':{'1':[1,9,26,29,37],'2':[1,15]},
                           '72':{'1':[1,4,11,15,17],'2':[0,10,14,17,18,19,20,21,22]},
                           '73':{'1':[2,8,9,10,13,14,15,16,18,19,26,28,29,22,30,32],'2':[9,17,23,24,26,31]},
                           '78':{'1':[0,5,19,11,14],'2':[0,7,8,31]},
                           '79':{'1':[0,1,2,6,7,8,10,14,15,16,27,25,28],'2':[3,0,2,7,9,11,13]},
                           '80':{'1':[7],'2':[1,7,18,20,21,23,26,28]},
                           '81':{'1':[0,10,11],'2':[0,13,11]},
                           '82':{'1':[0,15,13,16,17,18],'2':[3,12,15,16,17]},
                           '83':{'1':[0,9,17,7,4,10,14,19,22,28],'2':[0,4,14,12,13,7,15]},
                           '87':{'1':[1,13,15,19,25,33,31,4,7,9],'2':[1,2,8,10,14,13,15,19,24,20]},
                           '88':{'1':[1,4,9,14,16,18,20],'2':[0,4,6,9]},
                           '90':{'1':[1,9,12,19,25],'2':[0,7,12,21]},
                           '91':{'1':[2,17,15,20,21,22,23],'2':[1,13]},
                           '93':{'1':[0,9,6,10,15],'2':[0,8,16]},
                           '94':{'1':[0,6,9,18,14,16,19,25],'2':[0,5,10,13]},
                           '95':{'1':[0,14,15,17,18,19,20,21],'2':[1,3,15]},
                           '96':{'1':[0,17,20],'2':[0,10,21,22]},
                           '76':{'1':[0,14,15,21],'2':[1,10,11,22,27,28,29,32,33]},
                           '77':{'1':[0,1,14,28,45],'2':[0]},
                    
                          }

# normal ones 110 75 (48subjs)
# subjs=['02', '04', '07', '11', '12', '14', '16', '18', '19', '21', '22', '26', '28', '30',
#        '32', '34', '36', '37', '38', '40', '42', '50', '51', '52', '53', '54', '55', '56', '58', '59',
#        '60', '63', '65', '67', '68', '70', '72', '73', '78', '83', '87', '88', '90', '91', '93', '94', '95', '96']

#subject acquire specific parameters 130 75(10 subjs)
#SUBJECTS=(10 25 29 39 57 64 69 80 81 82)

#subject acquire specific parameters 140 75(4 subjs)
#SUBJECTS=(35 71 76 77)

#subject acquire specific parameters 160 75(1subj)
#SUBJECTS=(79)

# set-aside subjs 74,75,62

################################ step00: cut and filter data and concatenate 3 recording in one session ############

###### set up montage
montage_biosemi=mne.channels.read_montage(montage_fname)

###### preproc for each raw file
for subj in subj_list:
    psd_full_figs=[]
    psd_full_caption=[]
    session2conctn_list=[]
    ############### single subject report ###########################
    rep_fname = report_path + 'subj0' + subj+'ica report.h5'
    rep = mne.open_report(rep_fname)
    
    for session in session_list:
        autoR_figs=[]
        autoR_captions = []
        epochs_ASR_clean = get_epochs_ASR_clean(subj,session)[0]
        if epochs_ASR_clean == False:
            continue
        else:
            ############### step02 ICA components exclusion ##########################
            preProc_ica_fname = preProc_ica_path+'subj0'+subj+'session'+session+'preProc_ica.fif'
            ica = mne.preprocessing.read_ica(preProc_ica_fname)
            ica.exclude=compnts2exclude_ica[subj][session]
            ica.apply(epochs_ASR_clean)

            ################### step03 AutoRejection ##################################
            picks = mne.pick_types(epochs_ASR_clean.info, meg=False, eeg=True, stim=False,
                   eog=False)
            ar = AutoReject(picks=picks,random_state= 11,verbose='tqdm',n_jobs = 8)
            ar=ar.fit(epochs_ASR_clean)
            epochs_autorejected, reject_log = ar.transform(epochs_ASR_clean, return_log=True)
            reject_global_threshold = get_rejection_threshold(epochs_ASR_clean,ch_types=['eeg'],decim=2)
            while autorej_rate(epochs_autorejected) > 0.1:
                epochs_autorejected = epochs_ASR_clean.copy()
                reject_global_threshold['eeg']+=7e-06
                epochs_autorejected.drop_bad(reject=reject_global_threshold)
                print('reset global threshold to {}, rejecting rate turns to {}'.format(reject_global_threshold,
                                                                                        autorej_rate(epochs_autorejected)))

            
            autoR_figs.append(epochs_autorejected.plot_drop_log())
            autoR_captions.append('autorejecting rate')
#             autoR_figs.append(reject_log.plot_epochs(epochs_ASR_clean,scalings=150e-6))
#             autoR_captions.append('a glimpse of autorejecting epochs')

            threshes = ar.threshes_
            unit = r'uV'
            scaling = 1e6

            thres_hist=plt.figure(figsize=(6, 5))
            plt.tick_params(axis='x', which='both', bottom='off', top='off')
            plt.tick_params(axis='y', which='both', left='off', right='off')

            plt.hist(scaling * np.array(list(threshes.values())), 30,
                     color='g', alpha=0.4)
            plt.xlabel('Threshold (%s)' % unit)
            plt.ylabel('Number of sensors')
            plt.tight_layout()
            autoR_figs.append(thres_hist)
            autoR_captions.append('threshold histogram')

#             rep.add_figs_to_section(figs=psd_figs, captions=psd_captions, section = 'session'+session+'spectral plot', 
#                                     image_format = 'svg')
#             rep.add_figs_to_section(figs=ASR_figs, captions=ASR_captions, section = 'session'+session+'ASR plot', 
#                                     image_format = 'svg')
#             rep.add_figs_to_section(figs=ica_figs, captions=ica_captions, section = 'session'+session+'ica plot', 
#                                     image_format = 'svg')
            rep.add_figs_to_section(figs=autoR_figs, captions=autoR_captions, section = 'session'+session+
                                    'autoRejection plot', image_format = 'svg')

            session2conctn_list.append(epochs_autorejected)
            print('one session autoR is done')
    full_epochs_fname = full_epochs_path + 'subj0'+subj+'full_epo.fif'
    full_epochs_autorejected = mne.concatenate_epochs(session2conctn_list)
    psd_full_figs.append(full_epochs_autorejected.plot_psd())
    psd_full_caption.append('psd of epochs after preprocessing')
    rep.add_figs_to_section(figs=psd_full_figs, captions=psd_full_caption, section = 'preprocessed full epochs plot',
                            image_format = 'svg')
    full_epochs_autorejected.save(full_epochs_fname,overwrite = True)
        
#     f_report = report_path + 'subj0' + subj+'.html'
    rep.save(rep_fname, open_browser = False, overwrite = True)

        
