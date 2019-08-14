# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 17:49:16 2019

@author: gansheng.tan
"""

import pickle
import mne
import matplotlib.pyplot as plt
from mne.time_frequency import psd_array_multitaper
from scipy.integrate import simps
import numpy as np
import seaborn as sns
import pandas as pd
precleaned_epochs_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/'
fmin = 1
fmax = 100
def getBpAbsAndRelative4allChannels(epochs,rhythm):
    wavebands = {'alpha':[8,12],'theta':[3,7],'beta':[13,24],'lowG':[25,40],'highG':[60,90]}
    if rhythm in wavebands.keys():
        low,high =  wavebands[rhythm]
    else:
        print('not such rhythm')
    bpAbs_4Epochs=[]
    bpRelative_4Epochs=[]
    data = epochs.get_data(picks=['eeg'])
    for num_epochs in range(data.shape[0]):
        sf = epochs.info['sfreq']
        bpAbs_4allchannels = []
        bpRelative_4allchannels = []
        psd, freqs = psd_array_multitaper(data[num_epochs], sf, fmin = 1, fmax =100,
                          adaptive=True,normalization='full',verbose=0)
        psd= np.log10(psd*10e12)
        freq_res = freqs[1] - freqs[0]
        bp_total = simps(psd, dx=freq_res)
        idx_band = np.logical_and(freqs >= low, freqs <= high)
        bp_abs = simps(psd[:,idx_band], dx=freq_res)
        bp_relative = bp_abs/bp_total
        bpAbs_4Epochs.append(bp_abs)
        bpRelative_4Epochs.append(bp_relative)
    bpAbs_mean4Epochs = np.append([bpAbs_4Epochs[0]],bpAbs_4Epochs[1:],axis = 0).mean(axis=0)
    bpRelative_mean4Epochs = np.append([bpRelative_4Epochs[0]],bpRelative_4Epochs[1:],axis = 0).mean(axis=0)
    return bpAbs_mean4Epochs,bpRelative_mean4Epochs

# alpha bp clustering test - VD OP

subjs=['02', '04','07', '11', '12', '14', '16', '18', '19', '21', '22', '26', '28', '30',
       '32', '34', '36', '37', '38', '40', '42', '50', '51', '52', '53', '54', '55', '56',
       '58', '59','60', '63', '65', '67', '68', '70', '72', '73', '78', '83', '87', '88', 
       '90', '91', '93', '94', '95', '96','10','25','29','39','57','64','69','80','81','82',
       '35','71','79','76','77']

# states_codes={'VD':['111.0','112.0','121.0','122.0','131.0','132.0'],
#               'FA':['211.0','212.0','221.0','222.0','231.0','232.0'],
#               'OP':['311.0','312.0','321.0','322.0','331.0','332.0']}

states_codes={'VD':['111.0','112.0'],
              'FA':['211.0','212.0'],
              'OP':['311.0','312.0']}

bpAbs_mean4Epochs_VD4allsubjs = np.array([])
bpAbs_mean4Epochs_OP4allsubjs = np.array([])
for subj in subjs:
    precleaned_epochs_fname = precleaned_epochs_path + 'subj0'+subj+'full_epo.fif'
    precleaned_epochs = mne.read_epochs(precleaned_epochs_fname, preload=True)
    precleaned_epochs_VD = precleaned_epochs[states_codes['VD']]
    precleaned_epochs_OP = precleaned_epochs[states_codes['OP']]
    bpAbs_mean4Epochs_VD,bpRelative_mean4Epochs_VD = getBpAbsAndRelative4allChannels(precleaned_epochs_VD,'alpha')
    bpAbs_mean4Epochs_OP,bpRelative_mean4Epochs_OP = getBpAbsAndRelative4allChannels(precleaned_epochs_OP,'alpha')
    if len(bpAbs_mean4Epochs_VD4allsubjs)==0:
        bpAbs_mean4Epochs_VD4allsubjs = bpAbs_mean4Epochs_VD
        bpRelative_mean4Epochs_VD4allsubjs = bpRelative_mean4Epochs_VD
    else:
        bpAbs_mean4Epochs_VD4allsubjs = np.vstack((bpAbs_mean4Epochs_VD4allsubjs,bpAbs_mean4Epochs_VD))
        bpRelative_mean4Epochs_VD4allsubjs = np.vstack((bpRelative_mean4Epochs_VD4allsubjs,
                                                        bpRelative_mean4Epochs_VD))
        
    if len(bpAbs_mean4Epochs_OP4allsubjs)==0:
        bpAbs_mean4Epochs_OP4allsubjs = bpAbs_mean4Epochs_OP
        bpRelative_mean4Epochs_OP4allsubjs = bpRelative_mean4Epochs_OP
    else:
        bpAbs_mean4Epochs_OP4allsubjs = np.vstack((bpAbs_mean4Epochs_OP4allsubjs,bpAbs_mean4Epochs_OP))
        bpRelative_mean4Epochs_OP4allsubjs = np.vstack((bpRelative_mean4Epochs_OP4allsubjs,
                                                        bpRelative_mean4Epochs_OP))
        
bpAbs_mean4Epochs2test = [np.expand_dims(bpAbs_mean4Epochs_VD4allsubjs,axis=1),
                          np.expand_dims(bpAbs_mean4Epochs_OP4allsubjs,axis=1)]
bpRelative_mean4Epochs2test = [np.expand_dims(bpRelative_mean4Epochs_VD4allsubjs,axis=1),
                          np.expand_dims(bpRelative_mean4Epochs_OP4allsubjs,axis=1)]

with open("/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/VDOP_alpha_Abs.txt", "wb") as fp:   #Pickling
    pickle.dump(bpAbs_mean4Epochs2test, fp)
with open("/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/VDOP_alpha_Relative.txt", "wb") as fp:   #Pickling
    pickle.dump(bpAbs_mean4Epochs2test, fp)