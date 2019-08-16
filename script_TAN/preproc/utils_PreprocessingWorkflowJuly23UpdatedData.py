# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 14:29:23 2019

@author: gansheng.tan
"""



# # Below are personalised functions -> to py script and import


import warnings
import os
import mne
from utils_ASR import *
import os
import mne
from utils_ASR import *

<<<<<<< HEAD

def get_epochs_ASR_clean(subj,session):

    ##################### OS path in INSERM computer #####################################################################
    raw_data_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/raw_data/'
    montage_fname = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/raw_data/Biosemi64_MAS_EOG.locs'
=======
def get_epochs_ASR_clean(subj,session):

    ##################### OS path in INSERM computer #####################################################################
#     raw_data_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/raw_data/'
#     montage_fname = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/raw_data/Biosemi64_MAS_EOG.locs'
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
    # # report_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/report/'
    # full_epochs_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/'
    #
    ##################### OS path in cluster ######################################################################
<<<<<<< HEAD
#    raw_data_path = '/mnt/data/gansheng/raw_data/'
#    montage_fname = '/mnt/data/gansheng/raw_data/Biosemi64_MAS_EOG.locs'
=======
    raw_data_path = '/mnt/data/gansheng/raw_data/'
    montage_fname = '/mnt/data/gansheng/raw_data/Biosemi64_MAS_EOG.locs'
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
#     preProc_ica_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/preProc_ica/'
#     report_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/report/'
    # full_epochs_path = '/mnt/data/gansheng/preClean_data/'


    ########################################## Algorithme parameter ############################################
    cutoff = 10
    ########################################## Initialization parameter##########################################
    state_list = ['VD','FA','OP']
    power_freq_array = [50]
<<<<<<< HEAD
    reject_raw_data_session1 = {'74':['FA','OP','VD'],'62':['FA','OP','VD']}
    reject_raw_data_session2 = {'74':['FA','OP','VD'],'62':['FA','OP','VD']}

    # bad channel rejection is not apllied in the preproc, bad channels will be defined by eyes later
    #62,74 IS SO SPECIAL, MISS 75 76 77
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
=======
    reject_raw_data_session1 = {'29':['VD','FA','OP'],'30':['FA'],'36':['OP'],'74':['FA','OP','VD']}
    reject_raw_data_session2 = {'74':['FA','OP','VD'],'55':['VD']}

    # bad channel rejection is not apllied in the preproc, bad channels will be defined by eyes later
    #62,74 IS SO SPECIAL, MISS 75 76 77
    bad_channels={'02':{'1':['P2']},
                  '04':{'2':['FC6']},
                  '07':{'1':['Iz'],'2':['F8','T7']},
                  '10':{'2':['F2']},
                  '12':{'1':['P2']},
                  '14':{'1':['P2'],'2':['P2']},
                  '19':{'1':['P2','T7'],'2':['P2']},
                  '21':{'2':['Iz']},
                  '22':{'1':['T7','TP7','Iz','P1','POz','Pz'],'2':['T7','TP7','Iz']},
                  '25':{'1':['T8','FP1'],'2':['FC1','C3','PO4','F2','Pz']},
                  '26':{'1':['T8'],'2':['CPz']},
                  '28':{'2':['CP2','PO7','Oz','POz']},
                  '29':{'2':['F6']},
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                  '32':{'1':['P2'],'2':['P2','TP7']},
                  '34':{'1':['P2'],'2':['P2']},
                  '35':{'1':['T7','T8'],'2':['T8','PO8']},
                  '36':{'1':['P2','PO4'],'2':['P2','PO4']},
                  '37':{'1':['Iz']},
                  '38':{'2':['C5','FC5','TP8']},
<<<<<<< HEAD
                  '39':{'1':['P2','F8','AF8','Fp1','AF7'],'2':['P2','FT8','AF8','T8','P10']},
=======
                  '39':{'1':['P2'],'2':['P2','FT8']},
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
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
<<<<<<< HEAD
=======
                  '61':{'1':['P2']},
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                  '63':{'2':['PO8']},
                  '64':{'1':['C1']},
                  '65':{'1':['P2'],'2':['P2']},
                  '67':{'1':['FT8']},
                  '68':{'1':['P2'],'2':['P2']},
                  '70':{'1':['PO4','O2','FC3','FC5','F4','F6'],'2':['PO4','O2','FC5','FC3']},
                  '71':{'1':['P2','Iz'],'2':['P2','Iz','C1','Cz']},
                  '73':{'2':['FCz','FT8']},
<<<<<<< HEAD
                #'75':{'1':['C6','P2','FT8','AF8','CP1','P9','PO4','O2']},
                  '76':{'2':['T7']},
                  '77':{'1':['F6'],'2':['O1','Oz','F6','O2']},
                  '78':{'1':['P2'],'2':['P2']},
                  '79':{'1':['P2','POz'],'2':['P2','POz','T7','Fp1','AF7']},
=======
                  '78':{'1':['P2'],'2':['P2']},
                  '79':{'1':['P2','POz'],'2':['P2','POz','T7']},
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                  '81':{'1':['Iz','Oz','Pz','CPz','PO4','P2','POz'],'2':['Iz','Oz','POz','CPz','P2','PO4','FC1','C1','Pz']},
                  '82':{'1':['P2'],'2':['AFz']},
                  '83':{'1':['T7'],'2':['T7']},
                  '87':{'2':['P2']},
                  '88':{'1':['FC2','T8'],'1':['F4','P8','CP4']},
                  '90':{'1':['T7','P2'],'2':['P2']},
                  '91':{'1':['P2'],'2':['P2']},
<<<<<<< HEAD
                  '93':{'1':['FC5','Fp1','F3','PO4'],'2':['Fp1','F3','FC3','PO4']},
                  '94':{'1':['Fp1','F6','AF8','Fp2','T7','T8'],'2':['Fp1','F6','AF8','Fp2','T7']},
=======
                  '93':{'1':['PO4'],'2':['PO4']},
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                  '95':{'1':['P2'],'2':['P2']}
                 }
    reject_state = []
    conctn_list = []
    conctn_anno_list=[]
    montage_biosemi=mne.channels.read_montage(montage_fname)
    ASR_figs = []
    ASR_captions = []
    psd_figs=[]
    psd_captions=[]
    if subj in eval('reject_raw_data_session'+session).keys():
        reject_state = eval('reject_raw_data_session'+session)[subj]
        print("the rejected states of subject {} in session {} are {}".format(subj,session,reject_state))
    for state in state_list:
        raw_fname = raw_data_path + 'subj0'+subj+'_'+state+session+'_mast.set'
        if state in reject_state:
            continue
        elif os.path.isfile(raw_fname)==False:
            print('raw data file is missing: subject{}, session{}, state{}'.format(subj,session,state))

        else:
            raw = mne.io.read_raw_eeglab(raw_fname,montage_biosemi,verbose='INFO',preload=True,eog='auto')
            if raw.info['sfreq'] != 512:
                raw.resample(sfreq=512, npad = 'auto')
            events = mne.events_from_annotations(raw)
#             if subj == '36' and session == '2' and state == 'VD':
#                 raw.crop(tmin=220)
#                 events = mne.events_from_annotations(raw)
#                 events = events_time_forward(events,220)

            events_coding=events[1]
            events=np.asarray(events[0])  
<<<<<<< HEAD
            if '254' not in events_coding.keys():
                raw.annotations.append(onset=0,duration=0,description='254')
                events = mne.events_from_annotations(raw)
                events_coding=events[1]
                events=np.asarray(events[0])  
            if '255' not in events_coding.keys():
                raw.annotations.append(onset=events[-1][0]/512,duration=0,description='255')
=======
            if '254.0' not in events_coding.keys():
                raw.annotations.append(onset=0,duration=0,description='254.0')
                events = mne.events_from_annotations(raw)
                events_coding=events[1]
                events=np.asarray(events[0])  
            if '255.0' not in events_coding.keys():
                raw.annotations.append(onset=events[-1][0]/512,duration=0,description='255.0')
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                events = mne.events_from_annotations(raw)
                events_coding=events[1]
                events=np.asarray(events[0])  

<<<<<<< HEAD
            events_code_start = events_coding['254']
            start = events[events[:,2]==events_code_start][0][0]
            events_code_end = events_coding['255']
=======
            events_code_start = events_coding['254.0']
            start = events[events[:,2]==events_code_start][0][0]
            events_code_end = events_coding['255.0']
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
            stop = events[events[:,2]==events_code_end][0][0]

            raw_cut_filt = raw.copy()
            raw_cut_filt.crop(tmin = start/raw.info['sfreq'], tmax = stop/raw.info['sfreq'])
            raw_cut_filt.notch_filter(freqs=power_freq_array)
            raw_cut_filt.filter(l_freq=1,h_freq=100)
            psd_figs.append(raw_cut_filt.plot_psd(show = False))
            psd_captions.append('subject '+subj+"'s "+'psd plot after cut and filtering in session' 
                                +session+ ' during '+state+' state')
            ############ annotation engineering ################
            index_dlt=0
            for i in range(raw_cut_filt.annotations.__len__()):
<<<<<<< HEAD
                if (raw_cut_filt.annotations.__getitem__(i-index_dlt)['description']) not in ['131','132','255']:
=======
                if (raw_cut_filt.annotations.__getitem__(i-index_dlt)['description']) not in ['131.0','132.0','255.0']:
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                    raw_cut_filt.annotations.delete(i-index_dlt)
                    index_dlt+=1                       
                else: 
                    continue
            mne_annotation_recode_by_adding(session=session,state=state,annotations=raw_cut_filt.annotations)

            conctn_anno_list.append(raw_cut_filt.annotations)
            conctn_list.append(raw_cut_filt)

################### Concatenation process #################################
    if len(conctn_list)== 0:
        print('nothing 2 concatenate')
        return False,False,False,False,False
    else:        
        full_array = conctn_list[0]._data
        full_info = conctn_list[0].info
        del conctn_list[0]
        for raw2conctn in conctn_list:
            full_array = np.concatenate((full_array,raw2conctn._data),axis=1)
        raw_full = mne.io.RawArray(full_array,info = full_info)
        full_annotation = conctn_anno_list[0]
        del conctn_anno_list[0]
        for annos2conctn in conctn_anno_list:
            mne_annotation_postpone (pptime=full_annotation.__getitem__(full_annotation.__len__()-1)['onset'], 
                                         annotations=annos2conctn)
            full_annotation = full_annotation.__add__(annos2conctn)

        raw_full.set_annotations(full_annotation)
        if subj in bad_channels.keys() and session in bad_channels[subj].keys() :
            raw_full.info['bads']=bad_channels[subj][session]


###########raw_full now is for one session 

        ############### step01: epochs engineering - calibration-epochs-ASR #################################
        rawCalibAsr = raw_full.copy()
        rawCalibAsr = rawCalibAsr.crop(tmin=10,tmax=150)
        rawCalibAsr_noYW = rawCalibAsr.copy()
        rawCalibAsr._data,iirstate = YW_filter(rawCalibAsr._data, rawCalibAsr.info['sfreq'],None)
        rawVEOG= rawCalibAsr.copy()
        rawVEOG = rawVEOG.pick_channels(['VEOG'])
        VEOG_data = np.squeeze(rawVEOG.get_data())
<<<<<<< HEAD
        peak_locs, peak_eeg = mne.preprocessing.peak_finder(VEOG_data, thresh = 160e-6)
=======
        peak_locs, peak_eeg = mne.preprocessing.peak_finder(VEOG_data, thresh = 100e-6)
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
        lengthblink = 0.5*rawCalibAsr.info['sfreq']
        startremoveblink = peak_locs-(lengthblink/2)
        stopremoveblink = peak_locs+(lengthblink/2)
        NbsampCalibAsrWindow = len(VEOG_data)
        startremoveblink = np.abs((startremoveblink>0)*startremoveblink)
        stopremoveblink =  (stopremoveblink>NbsampCalibAsrWindow-1)*NbsampCalibAsrWindow + (stopremoveblink<NbsampCalibAsrWindow-1)*stopremoveblink
        Mask=np.zeros(NbsampCalibAsrWindow)
        for ix2remove in range(len(startremoveblink)):
            Mask[int(startremoveblink[ix2remove]):int(stopremoveblink[ix2remove])]=1
        rawCalibAsr_noYW.pick_types(eeg=True)
        rawdata_noblink = np.delete(rawCalibAsr_noYW.get_data(),np.where(Mask),axis=1)

        SignalCalib=np.delete(rawdata_noblink,np.where(np.abs(rawdata_noblink)>75e-6)[1],axis=1)
        ref_maxbadchannels = 0.2
        ref_tolerances = [-3.5,5.5]
        ref_wndlen = 1
        SignalClean,sample_mask = clean_windows(SignalCalib,rawCalibAsr.info['sfreq'],ref_maxbadchannels,ref_tolerances,ref_wndlen)
        SignalClean_raw = mne.io.RawArray(SignalClean,rawCalibAsr_noYW.info)
        ASR_figs.append(SignalClean_raw.plot(n_channels = 64, scalings = 150e-6))
        srate = rawCalibAsr.info['sfreq']
        cutoff = cutoff
        asr_state = asr_calibrate(SignalClean,srate,cutoff)
        if subj in bad_channels.keys() and session in bad_channels[subj].keys() :
            raw_bad = raw_full.copy()
            raw_bad.pick(raw_bad.info['bads'])
        raw4detect = raw_full.copy()
        raw4detect.pick_types(eeg=True)
        raw_full_eeg=raw_full.copy()
        raw_full_eeg.pick_types(eeg=True)
        raw_full_eog=raw_full.copy()
        raw_full_eog.pick_types(eog=True)
        raw4detect._data,iirstate = YW_filter(raw4detect._data,raw4detect.info['sfreq'],None)


        events = mne.events_from_annotations(raw_full)
        for i in range(len(events[0][:,2])):
            events[0][i][2]=int(float(dict_getValue(events[1],events[0][i][2])))

        for key in events[1].keys():
            events[1][key] = int(float(key))

        events_time_description = events[0]   
        i=0
        while i <len(events_time_description[:,2]):
            if events_time_description[i][2]==255:
                events_time_description=np.delete(events_time_description,i,0)
            else:
                i+=1
        events_dict=events[1]
<<<<<<< HEAD
        events_dict=removeItem_from_dict(events_dict,'255')
=======
        events_dict=removeItem_from_dict(events_dict,'255.0')
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0

        sfreq=raw_full.info['sfreq']   
        i=0
        overflow=False

        while overflow==False:
            if events_time_description[i,0]+sfreq*2>=events_time_description[-1][0]:
                overflow=True
            elif events_time_description[i+1,0]-events_time_description[i,0]>sfreq*2:
                events_time_description=np.insert(events_time_description,i+1,
                                                  [sfreq*2+events_time_description[i,0],0,events_time_description[i,2]],
                                                  axis=0)
                i+=1
            else:
                i+=1


        events=(events_time_description,events_dict)

        epochs4detect=mne.Epochs(raw4detect,events=events[0],event_id = events[1],tmin=0, tmax=2,preload=True)
        epochs_full=mne.Epochs(raw_full_eeg,events=events[0], event_id = events[1],tmin=0, tmax=2,preload=True)
        epochs_eog_raw_filt = mne.Epochs(raw_full_eog,events=events[0], event_id = events[1],tmin=0, tmax=2,preload=True)
        if subj in bad_channels.keys() and session in bad_channels[subj].keys() :
            epochs_bad_channels = mne.Epochs(raw_bad,events=events[0], event_id = events[1],tmin=0, tmax=2,preload=True) 
        Data4detect = epochs4detect.get_data()
        Data2correct = epochs_full.get_data()
        DataClean = np.zeros((Data2correct.shape))
        num_epochs2correct = 0
        num_epochscorrected = 0
        for i_epoch in range(Data2correct.shape[0]):
            num_epochs2correct+=1
            Epoch4detect = Data4detect[i_epoch,:,:] 
            Epoch2corr = Data2correct[i_epoch,:,:]    
            DataClean[i_epoch,:,:],reconstruct = asr_process_on_epoch(Epoch2corr,Epoch4detect,asr_state)

            if reconstruct ==True:
                num_epochscorrected +=1
        print('ASR correcting rate is {}'.format(num_epochscorrected/num_epochs2correct))
        ASR_captions.append('ASR signalClean of'+subj+' in session '+session+
                           ' with correcting rate '+'%.5f' % (num_epochscorrected/num_epochs2correct))

        epochs_ASR_clean =   mne.EpochsArray(DataClean,info=epochs_full.info,events=events[0],event_id = events[1])
        epochs_ASR_clean.add_channels([epochs_eog_raw_filt])
        if subj in bad_channels.keys() and session in bad_channels[subj].keys() :
            epochs_ASR_clean.add_channels([epochs_bad_channels])
            epochs_ASR_clean.interpolate_bads()
        return epochs_ASR_clean, psd_figs, psd_captions, ASR_figs, ASR_captions

<<<<<<< HEAD
    

=======
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0

def autorej_rate(epochs_autorejected):
    epoch_num = 0
    drop_epoch = 0
    for epoch in epochs_autorejected.drop_log:
        epoch_num+=1
        if not epoch:
            continue
        else:
            drop_epoch+=1
    rejecting_rate = drop_epoch/epoch_num
    return rejecting_rate

def events_time_forward(events,forward_time):
    ### forward_time in second
    for event in events[0]:
        event[0]=event[0]-forward_time*512
    return events

def removeItem_from_dict(d,key):
    r = dict(d)
    del r[key]
    return r

def dict_getValue(dictionary, search_value):
    for key,value in dictionary.items():
        if value==search_value:
            return key

def mne_annotation_empty(annotations):
    index_dlt = 0
    for i in range(annotations.__len__()):
        annotations.delete(i-index_dlt)
        index_dlt+=1
    return annotations

def mne_annotation_postpone (pptime, annotations):
    onset = []
    duration = []
    description = []
    index_dlt = 0
    for i in range(annotations.__len__()):
        onset.append(annotations.__getitem__(i)['onset']+pptime)
        duration.append(0.0)
        description.append(annotations.__getitem__(i)['description'])
    for i in range(annotations.__len__()):
        annotations.delete(i-index_dlt)
        index_dlt+=1
    annotations.append(onset=onset,duration=duration,description=description)
    print ('annotation time shift succeed')
    
        

def mne_annotation_recode_by_adding(session,state,annotations):
    onset = []
    duration = []
    description = []
    for i in range(annotations.__len__()):
<<<<<<< HEAD
        if annotations.__getitem__(i)['description'] in ['131','132']:
=======
        if annotations.__getitem__(i)['description'] in ['131.0','132.0']:
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
            onset,duration,description = mne_annotation_recode_info_extract(session=session,state=state,
                                                                        original_annotation = 
                                                                        annotations.__getitem__(i),
                                                                       onset=onset,duration=duration,
                                                                        description=description)
        else:
            continue
    index_dlt = 0
    for i in range(annotations.__len__()):
<<<<<<< HEAD
        if annotations.__getitem__(i-index_dlt)['description'] in ['131','132']:
=======
        if annotations.__getitem__(i-index_dlt)['description'] in ['131.0','132.0']:
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
            annotations.delete(i-index_dlt)
            index_dlt+=1
        else:
            continue
    onset.append(0.0)
    duration.append(0.0)
    description.append(mne_annotation_add_baseline(session=session,state=state))
    annotations.append(onset=onset,duration=duration,description=description)
    print ('annotation engineering succeed')
    return True

def mne_annotation_add_baseline(session,state):
    if session == '1':
        if state == 'VD':
            return '111.0'
        elif state == 'FA':
            return '211.0'
        elif state == 'OP':
            return '311.0'
        else:
            warnings.warn("unknown state detected", DeprecationWarning)
    elif session == '2':
        if state == 'VD':
            return '112.0'
        elif state == 'FA':
            return '212.0'
        elif state == 'OP':
            return '312.0'
        else:
            warnings.warn("unknown state detected", DeprecationWarning)
    else:
        warnings.warn("add baseline function only apply on rawfile having 2 sessions", DeprecationWarning)
    return '999.0'
        

def mne_annotation_recode_info_extract(session,state,original_annotation,onset,duration,description):
    if session =='1':
        if state == 'VD':
<<<<<<< HEAD
            if original_annotation['description']=='131':
=======
            if original_annotation['description']=='131.0':
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('121.0')
                
<<<<<<< HEAD
            elif original_annotation['description']=='132':
=======
            elif original_annotation['description']=='132.0':
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('131.0')
            else:
                print('this function only detect safe and threat period, please check original annotations')
        elif state == 'FA':
<<<<<<< HEAD
            if original_annotation['description']=='131':
=======
            if original_annotation['description']=='131.0':
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('221.0')

<<<<<<< HEAD
            elif original_annotation['description']=='132':
=======
            elif original_annotation['description']=='132.0':
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('231.0')
            else:
                print('this function only detect safe and threat period, please check original annotations')
        elif state == 'OP':
<<<<<<< HEAD
            if original_annotation['description']=='131':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('321.0')
            elif original_annotation['description']=='132':
=======
            if original_annotation['description']=='131.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('321.0')
            elif original_annotation['description']=='132.0':
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('331.0')
            else:
                print('this function only detect VD, FA, OP states, please check original annotations')
    elif session =='2':
        if state == 'VD':
<<<<<<< HEAD
            if original_annotation['description']=='131':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('122.0')
            elif original_annotation['description']=='132':
=======
            if original_annotation['description']=='131.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('122.0')
            elif original_annotation['description']=='132.0':
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('132.0')
            else:
                print('this function only detect safe and threat period, please check original annotations')
        elif state == 'FA':
<<<<<<< HEAD
            if original_annotation['description']=='131':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('222.0')
            elif original_annotation['description']=='132':
=======
            if original_annotation['description']=='131.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('222.0')
            elif original_annotation['description']=='132.0':
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('232.0')
            else:
                print('this function only detect safe and threat period, please check original annotations')
        elif state == 'OP':
<<<<<<< HEAD
            if original_annotation['description']=='131':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('322.0')
            elif original_annotation['description']=='132':
=======
            if original_annotation['description']=='131.0':
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('322.0')
            elif original_annotation['description']=='132.0':
>>>>>>> 7d5e7154242b95c4072aefc73cf20394ce2e22a0
                onset.append(original_annotation['onset'])
                duration.append(original_annotation['duration'])
                description.append('332.0')
            else:
                print('this function only detect VD, FA, OP states, please check original annotations')
    else:
        print('3rd session dected, please check annotations')
    return(onset,duration,description)
        
