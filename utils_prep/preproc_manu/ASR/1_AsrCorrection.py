# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:39:43 2019

@author: Manu
"""


import mne
from mne import io
import sys
sys.path.append('C:/_MANU/_U821/Python_Dev/')
import scipy
from util import tools,asr,raw_asrcalibration
import numpy as np


import matplotlib.pyplot as plt
from mne.viz import plot_evoked_topo


fname = 'C:/_MANU/_U821/_wip/ContextOdd/raw/ANDNI_0001.vhdr'

raw = io.read_raw_brainvision(fname, preload = False)
picks_eeg = mne.pick_types(raw.info, meg=False, eeg=True, eog=False,stim=False, exclude='bads')
ListChannels  = np.array(raw.info['ch_names'])
montage = mne.channels.read_montage(kind='standard_1020',ch_names=ListChannels[picks_eeg])
raw = io.read_raw_brainvision(fname,  montage=montage, preload = True)
picks_eeg = mne.pick_types(raw.info, meg=False, eeg=True, eog=False,stim=False, exclude='bads')
raw =raw.pick_types( meg=False, eeg=True, eog=False,stim=True, exclude='bads')



# ASR Calibration    


raworig_Data= raw._data
l_freq = 2
h_freq = 20
Wn = [l_freq/(raw.info['sfreq']/2.), h_freq/(raw.info['sfreq']/2.) ]
b, a = scipy.signal.iirfilter(N=2, Wn=Wn, btype = 'bandpass', analog = False, ftype = 'butter', output = 'ba')
raw._data[picks_eeg,:]=scipy.signal.lfilter(b, a, raworig_Data[picks_eeg,:], axis = 1, zi = None)




rawCalibAsr=raw.copy()
tmin = 30
tmax = 60 #s
rawCalibAsr = rawCalibAsr.crop(tmin=tmin,tmax=tmax)


ChanName4VEOG = ['Fp1','Fp2'] # 2 VEOG
cutoff = 5 # Makoto preprocessing says best between 10 and 20 https://sccn.ucsd.edu/wiki/Makoto%27s_preprocessing_pipeline#Alternatively.2C_cleaning_continuous_data_using_ASR_.2803.2F26.2F2019_updated.29

Yule_Walker_filtering = True
state = raw_asrcalibration.raw_asrcalibration(rawCalibAsr,ChanName4VEOG, cutoff,Yule_Walker_filtering)













# ASR process on epoch
event_id = {'Std': 1, 'Dev': 2}
events_orig,_ =  mne.events_from_annotations(raw)
ixdev = np.array(np.where(events_orig[:,2]==2))
ixstd= ixdev-1

events = events_orig[np.sort(np.array(np.hstack((ixstd , ixdev)))),:]
events = np.squeeze(events, axis=0)
tmin, tmax = -0.2, 0.5
raw4detect = raw.copy()
raw4detect._data,iirstate = asr.YW_filter(raw._data,raw.info['sfreq'],None) ## HERE
epochs4Detect = mne.Epochs(raw4detect, events=events, event_id=event_id, tmin=tmin,tmax=tmax, proj=True,baseline=None, reject=None, picks=picks_eeg)
epochs_filt = mne.Epochs(raw, events=events, event_id=event_id, tmin=tmin,tmax=tmax, proj=None,baseline=None, reject=None, picks=picks_eeg)


Data4detect = epochs4Detect.get_data()
Data2Correct = epochs_filt.get_data()
DataClean = np.zeros((Data2Correct.shape))
for i_epoch in range(Data4detect.shape[0]):
    EpochYR = Data4detect[i_epoch,:,:]
    Epoch2Corr = Data2Correct[i_epoch,:,:]    
    DataClean[i_epoch,:,:] = asr.asr_process_on_epoch(EpochYR,Epoch2Corr,state)
    
    



epochs_clean  =   mne.EpochsArray(DataClean,info=epochs_filt.info,events=events,event_id=event_id)
srate = raw.info['sfreq']
evoked_std = epochs_filt['Std'].average(picks=picks_eeg)
evoked_dev = epochs_filt['Dev'].average(picks=picks_eeg)
    
evoked_clean_std = epochs_clean['Std'].average(picks=picks_eeg)
evoked_clean_dev = epochs_clean['Dev'].average(picks=picks_eeg)

evoked_clean_std.first=-200
evoked_clean_std.last= tmax*srate

evoked_clean_dev.first=-200
evoked_clean_dev.last= tmax*srate

evoked_clean_std.times= np.around(np.linspace(-0.2, tmax, num=DataClean.shape[2]),decimals=3)
evoked_clean_dev.times= np.around(np.linspace(-0.2, tmax, num=DataClean.shape[2]),decimals=3)

evokeds = [evoked_std, evoked_dev, evoked_clean_std, evoked_clean_dev]

colors = 'blue', 'red','steelblue','magenta'


plot_evoked_topo(evokeds, color=colors, title='Std Dev', background_color='w')
plt.show()


evoked_clean_MMN=evoked_clean_std.copy()
evoked_clean_MMN.data = (evoked_clean_dev.data - evoked_clean_std.data)

evoked_MMN =evoked_clean_MMN.copy()
evoked_MMN.data = (evoked_dev.data-evoked_std.data)
evokeds_MMN= [evoked_clean_MMN,evoked_MMN]
colors = 'red', 'black'
plot_evoked_topo(evokeds_MMN, color=colors, title='MMN', background_color='w')
plt.show()




kwargs = dict(times=np.arange(-0.1, 0.40, 0.025), vmin=-1.5, vmax=1.5, layout='auto',
              head_pos=dict(center=(0., 0.), scale=(1., 1.)))

evoked_MMN.plot_topomap(**kwargs)
evoked_clean_MMN.plot_topomap(**kwargs)


















