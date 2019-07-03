# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 18:23:58 2018

@author: Manu
"""

import mne
from mne import io
from mne.preprocessing import ICA
from util import tools

def FitIcaEpoch(src_directory, src_fname,  dest_directory, picks):
    fname = src_directory + src_fname
    
    # Read raw data
    raw =   io.read_raw_fif(fname,preload=True)
    
    # Filter raw data to remove mains and DC offset   
    iir_Butter_params = dict(order=2, ftype='butter', output='sos') 
    raw.filter(l_freq = 0.1, h_freq = 20.,method = 'iir', iir_params=iir_Butter_params)
    
    # Epoch filtered data on stimuli of interest
    events = mne.find_events(raw)
    tmin, tmax = -0.1, 1
    
    
    
    epochs = mne.Epochs(raw, events=events, event_id=None, tmin=tmin,tmax=tmax, preload=True,proj=True,baseline=None, reject=None, picks=picks)
    
    PercentageOfEpochsRejected = 2.0
    ThresholdPeak2peak = tools.RejectThresh(epochs,PercentageOfEpochsRejected)
    
    reject = {'eeg': ThresholdPeak2peak}    
    epochs = mne.Epochs(raw, events=events, event_id=None, tmin=tmin,tmax=tmax, preload=True,proj=True,baseline=None, reject=reject, picks=picks)
   
    # Fit ica on epochs without big artifact
    ica = ICA(n_components=raw.info['nchan']-1, method='fastica',random_state = 30).fit(epochs, picks=picks, decim = 10)
    
    # Store ica solution
    save_fname = dest_directory + src_fname[:-8] +"-ica.fif"
    ica.save(save_fname)
    
    
    
def FitIcaRaw(src_directory, src_fname,  dest_directory, picks):
    fname = src_directory + src_fname
    
    # Read raw data
    raw =   io.read_raw_fif(fname,preload=True)
    
    # Filter raw data to remove mains and DC offset   
    iir_Butter_params = dict(order=2, ftype='butter', output='sos') 
    raw.filter(l_freq = 0.1, h_freq = 20.,method = 'iir', iir_params=iir_Butter_params)       
    
    # Fit ica on epochs without big artifact
    ica = ICA(n_components=raw.info['nchan']-1, method='fastica',random_state = 30).fit(raw, picks=picks, decim = 10)
    
    # Store ica solution
    save_fname = dest_directory + src_fname[:-8] +"-ica.fif"
    ica.save(save_fname)