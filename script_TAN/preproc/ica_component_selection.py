#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================
ICA componont visualisation on Enrico data using MNE
===============================================

This is a parallel process of preprocessing, the goal is to visualize the ICA components and select rejecting ones from them:
        ====> input : ica file for one subject
        1)detect eog components and visualize the correlation
        2)detect other types of artifacts
        3)give out excluding components for each subject
        ====> output = fif file that save ica object and a rejecting component dict


Note: exploration version

Suggestions:
1) can be used as a function to be included in automated process
2) 

Updated on June 2019

@author: Gansheng TAN aegean0045@outlook.com, Françoise Lecaignard, francoise.lecaignard@inserm.fr
"""
##############################################################  Set-up ######################################################################

import mne
import importlib
import numpy as np
import os
from autoreject import AutoReject
from autoreject import get_rejection_threshold 
mne.set_log_level('WARNING')


##########################
# ICA visualisation
##########################

def ica_component_selection(subj,state,HEOG_corr=False, VEOG_corr=True,second_check=False):
    
    
    ########################### 
    # Variable definition
    ############################

    ################ inserm computer #####################################
    raw_data_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/raw_data/'
    filt_data_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/filt_data/'
    ica_file_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/ica_file/'
    montage_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/raw_data/Biosemi64_MAS_EOG.locs'
    epochs_autorejed_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/autorejed_epochs_data/'
    epochs_icaexclded_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/icaexclded_epochs_data/'
    ##################################################################################################

    ########### local laptop###################################################
    # raw_data_path = 'E:/important file/INSERM/data_eeglab/raw_data/'
    # filt_data_path = 'E:/important file/INSERM/data_eeglab/filt_data/'
    # ica_file_path = 'E:/important file/INSERM/data_eeglab/ica_file/'
    # montage_path = "C:/Users/aegea/python/INSERM/raw_data/Biosemi64_MAS_EOG.locs"
    ############################################################################
    
    reject_raw_data={'07':['OP1','OP2'], '10':['FA1','VD1','VD2'], '21':['VD1','FA2','VD2'],
                '22':['OP2'], '36':['OP1'], '57':['OP2','FA2'], '82':['FA2','OP2','VD2']}
    #example: reject_raw_data={'94':['FA1']}
    bad_channel={'94':{'FA1':['Pz']}}
    
    ########################### visualization parameters #########################
    
    ###########################
    
    
    
    ###################################### read epochs ##################################################
    
    
    if subj in reject_raw_data.keys() & state in reject_raw_data[subj]:
        print("the recording file is of low quality, we do not do ica on this raw data")
    else:
        if second_check ==True:
            filename_epochs = epochs_icaexclded_path+ 'subj0'+subj+'_'+state+'_icaexclded_epo.fif'
            montage_biosemi=mne.channels.read_montage(montage_path)
            epochs = mne.read_epochs(filename_epochs,preload=True)
        else:
            filename_epochs = epochs_autorejed_path + 'subj0'+subj+'_'+state+'_autorejed_epo.fif'
            montage_biosemi=mne.channels.read_montage(montage_path)
            epochs = mne.read_epochs(filename_epochs,preload=True)
    
        
        f_ica = ica_file_path + 'subj0'+subj+'_'+state+'_component_ica.fif'
        ica = mne.preprocessing.read_ica(f_ica)
    

    
    
        ############### EOG channel verification
        title = 'Sources related to %s artifacts (red)'
        # detect EOG by correlation
        n_max_eog = 1

        if VEOG_corr ==True:
            eog_inds, scores = ica.find_bads_eog(epochs,ch_name='VEOG')
            ica.plot_scores(scores, exclude=eog_inds, title=title % 'Veog', labels='Veog')

            show_picks = np.abs(scores).argsort()[::-1][:5]
            

            #ica.plot_sources(epochs.average(), exclude=eog_inds, title=title % 'Veog')
            if eog_inds ==[]:
                print("no Veog component is found, data is cleaned from eog")
            else:
                ica.plot_components(eog_inds, title=title % 'Veog', colorbar=True)
                eog_inds = eog_inds[:n_max_eog]


        if HEOG_corr ==True:
            eog_inds, scores = ica.find_bads_eog(epochs,ch_name='HEOG')
            ica.plot_scores(scores, exclude=eog_inds, title=title % 'Heog', labels='Heog')

            show_picks = np.abs(scores).argsort()[::-1][:5]

#             ica.plot_sources(epochs, exclude=eog_inds, title=title % 'Heog')
            if eog_inds ==[]:
                print("no Heog component is found, data is cleaned from eog")
            else:
                ica.plot_components(eog_inds, title=title % 'Heog', colorbar=True)
                eog_inds = eog_inds[:n_max_eog]
            #ica.exclude += eog_inds   add eog artifact in the exclude list.

#             eog_evoked = mne.preprocessing.create_eog_epochs(epochs, tmin=-.5, tmax=.5, picks=picks).average()
#             ica.plot_sources(eog_evoked, exclude=eog_inds)  # plot EOG sources + selection   this plot is to plot ica components
#             ica.plot_overlay(eog_evoked, exclude=eog_inds)  # plot EOG cleaning  this plot is compare channels' signal
#                                                                             #after and before removing eog components

        # check the amplitudes do not change
        #ica.plot_overlay(raw_filt)  # EOG artifacts remain # not very useful, we do not gain a lot of info from it
        
        ica.get_sources(epochs).plot(picks='all')
# fail to plot only left components
#         if view_compn_exclded == True:
#             f_ica = ica_file_path + 'subj0'+subj+'_'+state+'_component_exclded_ica.fif'
#             ica = mne.preprocessing.read_ica(f_ica)
#         else:
#             print("plotting all ica components")
        ica.plot_components(inst=epochs)
        #ica.plot_properties()
        print("select the removing components and fill the ica_selection dict")
    return True

 
