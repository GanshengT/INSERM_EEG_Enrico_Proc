# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 19:17:39 2018

@author: Manu
"""

import os
import mne
from mne import io
from mne.preprocessing import ICA
from mne.preprocessing import create_eog_epochs
from mne.preprocessing.ica import corrmap  # noqa

import numpy as np
import matplotlib.pyplot as plt
from util import tools

def VirtualEog(raw_fname, ica_fname, fig_directory, pick, ChanName4VEOG, ChanName4HEOG_l,ChanName4HEOG_r):
    
    # Read raw data
    head, raw_f = os.path.split(raw_fname)
    raw =   io.read_raw_fif(raw_fname,preload=True)
    
#    FlagVEOG = False
#    FlagHEOG_L = False
#    FlagHEOG_R = False
#    FlagHEOG = True 
#    FlagICAscore = True
#    
#    # Create virtual Vertical EOG 
#    if ChanName4VEOG is not None:               
#        rawSelecChan4Veog  = raw.copy().pick_channels(ChanName4VEOG)
#        rawVEogData = np.zeros((1, rawSelecChan4Veog.n_times), dtype='float')
#        rawVEogData[0,:] = (rawSelecChan4Veog.get_data(picks=range(len(ChanName4VEOG))).sum(axis=0))
#        FlagVEOG = True
#
#     # Create virtual Horizontal EOG
#    if (ChanName4HEOG_l is not None) : 
#        rawSelecChan4Heog_l = raw.copy().pick_channels(ChanName4HEOG_l)
#        rawHEogL_Data = np.zeros((1, rawSelecChan4Heog_l.n_times), dtype='float')
#        rawHEogL_Data[0,:] = (rawSelecChan4Heog_l.get_data(picks=range(len(ChanName4HEOG_l))).sum(axis=0))
#        FlagHEOG_L = True
#        
#        
#    if (ChanName4HEOG_r is not None):
#        rawSelecChan4Heog_r = raw.copy().pick_channels(ChanName4HEOG_r)
#        rawHEogR_Data = np.zeros((1, rawSelecChan4Heog_r.n_times), dtype='float')
#        rawHEogR_Data[0,:] = (rawSelecChan4Heog_r.get_data(picks=range(len(ChanName4HEOG_r))).sum(axis=0))
#        FlagHEOG_R = True
#        
#        
#    if FlagHEOG_L:
#        if FlagHEOG_R:
#            rawHEogData = rawHEogL_Data - rawHEogR_Data
#        else:
#            rawHEogData = rawHEogL_Data
#    else:
#        if FlagHEOG_R:
#            rawHEogData = rawHEogR_Data
#        else:
#            FlagHEOG = False
#    
#    rawWithVirtEOG = raw.copy()
#    if FlagVEOG:
#        infoVEog = mne.create_info(['VEOG'], rawWithVirtEOG.info['sfreq'], ['eeg'])
#        VEogRawArray  = mne.io.RawArray(rawVEogData, infoVEog)
#        rawWithVirtEOG.add_channels([VEogRawArray], force_update_info=True)
#            
#    if FlagHEOG:
#        infoHEog = mne.create_info(['HEOG'], rawWithVirtEOG.info['sfreq'], ['eeg'])
#        HEogRawArray  = mne.io.RawArray(rawHEogData, infoHEog)
#        rawWithVirtEOG.add_channels([HEogRawArray], force_update_info=True)


    rawWithVirtEOG = tools.AddVirtualEogChannels(raw,ChanName4VEOG,ChanName4HEOG_l,ChanName4HEOG_r)

    picks_eeg = mne.pick_types(rawWithVirtEOG.info, meg=False, eeg=True, eog=False,stim=False, exclude='bads')
    
    
    iir_Butter_params = dict(order=2, ftype='butter', output='sos')  
    rawWithVirtEOG.filter(l_freq = 0.1, h_freq = 20.,method = 'iir', iir_params=iir_Butter_params,picks=picks_eeg)
    dict_scaling=dict(eeg=100e-6)
    rawWithVirtEOG.plot(duration=20,n_channels=rawWithVirtEOG.info['nchan']-1,scalings=dict_scaling)
    
            
    # Epoch filtered data on stimuli of interest
    events = mne.find_events(rawWithVirtEOG)
    tmin, tmax = -0.1, 1
    epochs = mne.Epochs(rawWithVirtEOG, events=events, event_id=None, tmin=tmin,tmax=tmax, preload=True,proj=True,baseline=None, reject=None, picks=picks_eeg)
   
    PercentageOfEpochsRejected = 2.0
    ThresholdPeak2peak = tools.RejectThresh(epochs,PercentageOfEpochsRejected)
    reject = {'eeg': ThresholdPeak2peak}   
    
    
    Veog_inds=[]
    Heog_inds=[]
    
    ica_directory, ica_f = os.path.split(ica_fname)
    ica = mne.preprocessing.read_ica(ica_fname)
    
    if FlagVEOG:
        Veog_average = create_eog_epochs(rawWithVirtEOG,ch_name ='VEOG', reject=reject, picks=picks_eeg,verbose=False).average()
        Veog_inds, Veog_scores = ica.find_bads_eog(epochs,ch_name ='VEOG')
        
    if FlagHEOG:
        Heog_average = create_eog_epochs(rawWithVirtEOG,ch_name ='HEOG', reject=reject, picks=picks_eeg,verbose=False).average()
        Heog_inds, Heog_scores = ica.find_bads_eog(epochs,ch_name ='HEOG')



    ica.plot_sources(rawWithVirtEOG,picks= Veog_inds+Heog_inds,stop  =30 )
    fnamesavefigICAcomponents = fig_directory + raw_f[:-8] + "-icacomp.jpg"
    plt.savefig(fnamesavefigICAcomponents)


    ica.plot_scores(Veog_scores, exclude=Veog_inds)  # look at r scores of components
    fnamesavefigIcaScoreVEOG = fig_directory + raw_f[:-8] + "-icascoreVeog.jpg"
    plt.savefig(fnamesavefigIcaScoreVEOG)
    
    
    ica.plot_scores(Heog_scores, exclude=Heog_inds)  # look at r scores of components
    fnamesavefigIcaScoreHEOG = fig_directory + raw_f[:-8] + "-icascoreHeog.jpg"
    plt.savefig(fnamesavefigIcaScoreHEOG)
    
    
    eog_inds = Veog_inds+Heog_inds
    ica.plot_components(eog_inds,ch_type='eeg')
    fnamesavefigICATopocomponents = fig_directory + raw_f[:-8] + "-icatopo.jpg"
    plt.savefig(fnamesavefigICATopocomponents)
    
    if FlagVEOG:
        if FlagHEOG:
            IcaScore2save = {'Veog_scores': Veog_scores[Veog_inds], 'Veog_inds' :Veog_inds,  'Heog_scores': Heog_scores[Heog_inds], 'Heog_inds' :Heog_inds}
        else:
            IcaScore2save = {'Veog_scores': Veog_scores[Veog_inds], 'Veog_inds' :Veog_inds}
    else:
        if FlagHEOG:
            IcaScore2save = {'Heog_scores': Heog_scores[Heog_inds], 'Heog_inds' :Heog_inds}
        else:
            FlagICAscore = False
            
    if FlagICAscore:
        SaveFileNameicascore = os.path.join(ica_directory,raw_f[:-8] +  "-scoreica.npy")
        np.save(SaveFileNameicascore, IcaScore2save)    
        
    SaveFileNameicaWeightVar = os.path.join(ica_directory,raw_f[:-8] +  "-icaweightsvar.npy")
    IcaWeightsVar2save = {'LogVarIcWeights' : np.log(tools.ComputeVarICAWeights(ica))}
    np.save(SaveFileNameicaWeightVar, IcaWeightsVar2save)    
    
    
    ica.exclude.extend(eog_inds)
    # Store ica solution
    save_fname = os.path.join(ica_directory,  raw_f[:-8] +"-icaexclude.fif")
    ica.save(save_fname)
    
    
    
    
    
    
def CorrmapWithTemplate(List_IcaFifFiles,templateVeog,templateHeog):
    
    icas = [mne.preprocessing.read_ica( fname) for fname in List_IcaFifFiles]
    
    if (templateVeog is not None):
        corrmap(icas, template=templateVeog, label='blinks',show=True, threshold=0.95, ch_type='eeg',verbose=False)
        for ica in icas:
            ica.exclude.extend(ica.labels_['blinks'])
            
            
    
    if (templateHeog is not None):
        corrmap(icas, template=templateHeog, label='saccades',show=True, threshold=0.95, ch_type='eeg',verbose=False)
        for ica in icas:
                ica.exclude.extend(ica.labels_['saccades'])
                
                
    if ((templateVeog is not None) or (templateHeog is not None)):
        isess = 0
        for ica in icas:
                fname = List_IcaFifFiles[isess]
                save_fname = fname[:-8] +"-icaexcludecorrmap.fif"
                ica.save(save_fname)
                isess = isess + 1 
