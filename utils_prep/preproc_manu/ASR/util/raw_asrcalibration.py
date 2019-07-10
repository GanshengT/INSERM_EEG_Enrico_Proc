# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 10:22:01 2019

@author: Manu
"""

import sys
sys.path.append('C:/_MANU/_U821/Python_Dev/')
from util import tools,asr
import numpy as np

from mne.preprocessing import peak_finder




def raw_asrcalibration(raw,ChanName4VEOG, cutoff,Yule_Walker_filtering):
# Compute ASR calibration from raw object
    #
    # raw : The raw data used for ASR calibration
    # ChanName4VEOG : Nome of channels to estimate a vertical EOG (blink detection) ex : ChanName4VEOG = ['Fp1','Fp2']
    # cutoff : parameter that determines the rejection threshold
    # Yule_Walker_filtering : Apply Yule_Walker filter or not 
      
    # Apply Yule-Walker filter
    rawCalibAsr=raw.copy()
    if Yule_Walker_filtering:
        rawCalibAsr._data,iirstate = asr.YW_filter(rawCalibAsr._data,rawCalibAsr.info['sfreq'],None)
    
    
    ## Blink detection
    # compute virtual vertical EOG to detect blinks
    ChanName4HEOG_l = None
    ChanName4HEOG_r = None
    rawWithVirtEog = tools.AddVirtualEogChannels(rawCalibAsr,ChanName4VEOG,ChanName4HEOG_l,ChanName4HEOG_r)
    
    # Find time location of blinks
    rawVEOG =rawWithVirtEog.pick_channels(['VEOG'])
    VEOG_data = np.squeeze(rawVEOG.get_data())
    peak_locs, peak_eeg = peak_finder(VEOG_data,thresh=75e-6)
    
    lengthblink = 0.5*raw.info['sfreq']; #500ms (in ms because the findpeaks return values in ms)
    startremoveblink = peak_locs-(lengthblink/2)
    stopremoveblink = peak_locs+(lengthblink/2)
    NbsampCalibAsrWindow = len(VEOG_data)
    startremoveblink = np.abs((startremoveblink>0)*startremoveblink)
    stopremoveblink =  (stopremoveblink>NbsampCalibAsrWindow-1)*NbsampCalibAsrWindow + (stopremoveblink<NbsampCalibAsrWindow-1)*stopremoveblink
    Mask=np.zeros(NbsampCalibAsrWindow)
    for ix2remove in range(len(startremoveblink)):
        Mask[int(startremoveblink[ix2remove]):int(stopremoveblink[ix2remove])]=1
    
    rawdata_noblink = np.delete(raw.get_data(),np.where(Mask),axis=1)
    SignalCalib=np.delete(rawdata_noblink,np.where(np.abs(rawdata_noblink)>50e-6)[1],axis=1)

  
    ref_maxbadchannels = 0.2;
    ref_tolerances = [-3.5,5.5];
    ref_wndlen = 1;
    
    SignalClean,sample_mask = asr.clean_windows(SignalCalib,raw.info['sfreq'],ref_maxbadchannels,ref_tolerances,ref_wndlen);
    srate = raw.info['sfreq']
    
    
    state = asr.asr_calibrate(SignalClean,srate,cutoff)
    return state