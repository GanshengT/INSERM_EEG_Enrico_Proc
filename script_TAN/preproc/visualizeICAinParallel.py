# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 17:32:41 2019

@author: gansheng.tan
"""

##### Visualize ica components and exclude them subjectively
from mne.preprocessing import read_ica
#### GET EPOCHS TO VISUALIZE ICA
import os
import mne
from utils_ASR import *
from utils_PreprocessingWorkflowJuly23UpdatedData import *
# %matplotlib qt

def h5ToHtml(rep_fname):
    import mne
    rep = mne.open_report(rep_fname)
    rep_fname = rep_fname[:-1]
    rep_fname = rep_fname[:-1]
    rep.save(rep_name+'html',open_browser=True, overwrite=False)

def get_epochs_ASR_clean2seeICA(subj,session):

    ##################### OS path in INSERM computer #####################################################################
    raw_data_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/raw_data/'
    montage_fname = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/raw_data/Biosemi64_MAS_EOG.locs'
#     report_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/report/'
    # full_epochs_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/'
    #
    ##################### OS path in cluster ######################################################################
#     raw_data_path = '/mnt/data/gansheng/raw_data/'
#     montage_fname = '/mnt/data/gansheng/raw_data/Biosemi64_MAS_EOG.locs'
#     preProc_ica_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/preProc_ica/'
#     report_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/report/'
    # full_epochs_path = '/mnt/data/gansheng/preClean_data/'


    ########################################## Algorithme parameter ############################################
    cutoff = 10
    ########################################## Initialization parameter##########################################
    state_list = ['VD','FA','OP']
    power_freq_array = [50]
    reject_raw_data_session1 = {'74':['FA','OP','VD'],'62':['FA','OP','VD']}
    reject_raw_data_session2 = {'74':['FA','OP','VD'],'62':['FA','OP','VD']}

    # bad channel rejection is not apllied in the preproc, bad channels will be defined by eyes later
    #62,74 IS SO SPECIAL, MISS 75 76 77
    bad_channels={'02':{'1':['P2']},
                  '04':{'2':['FC6']},
                  '07':{'1':['Iz'],'2':['F8','T7']},
                  '10':{'2':['F2']},
                  '12':{'1':['P2']},
                  '14':{'1':['P2'],'2':['P2']},
                  '16':{'1':['Fp1','AF7','Fp2','AF8'],'2':['Fp1','AF7','Fp2','AF8','F6']},
                  '19':{'1':['P2','T7'],'2':['P2']},
                  '21':{'2':['Iz']},
                  '22':{'1':['T7','TP7','Iz','P1','POz','Pz'],'2':['T7','TP7','Iz']},
                  '25':{'1':['T8','Fp1'],'2':['FC1','C3','PO4','F2','Pz']},
                  '26':{'1':['T8'],'2':['CPz']},
                  '28':{'2':['CP2','PO7','Oz','POz']},
                  '29':{'2':['F6']},
                  '32':{'1':['P2'],'2':['P2','TP7']},
                  '34':{'1':['P2'],'2':['P2']},
                  '35':{'1':['T7','T8'],'2':['T8','PO8']},
                  '36':{'1':['P2','PO4'],'2':['P2','PO4']},
                  '37':{'1':['Iz']},
                  '38':{'2':['C5','FC5','TP8']},
                  '39':{'1':['P2'],'2':['P2','FT8']},
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
                  '61':{'1':['P2']},
                  '63':{'2':['PO8']},
                  '64':{'1':['C1']},
                  '65':{'1':['P2'],'2':['P2']},
                  '67':{'1':['FT8']},
                  '68':{'1':['P2'],'2':['P2']},
                  '70':{'1':['PO4','O2','FC3','FC5','F4','F6'],'2':['PO4','O2','FC5','FC3']},
                  '71':{'1':['P2','Iz'],'2':['P2','Iz','C1','Cz']},
                  '73':{'2':['FCz','FT8']},
                  '78':{'1':['P2'],'2':['P2']},
                  '79':{'1':['P2','POz'],'2':['P2','POz','T7','Fp1','AF7']},
                  '81':{'1':['Iz','Oz','Pz','CPz','PO4','P2','POz'],'2':['Iz','Oz','POz','CPz','P2','PO4','FC1','C1','Pz']},
                  '82':{'1':['P2'],'2':['AFz']},
                  '83':{'1':['T7'],'2':['T7']},
                  '87':{'2':['P2']},
                  '88':{'1':['FC2','T8'],'1':['F4','P8','CP4']},
                  '90':{'1':['T7','P2'],'2':['P2']},
                  '91':{'1':['P2'],'2':['P2']},
                  '93':{'1':['PO4'],'2':['PO4']},
                  '95':{'1':['P2'],'2':['P2']}
                 }
    reject_state = []
    conctn_list = []
    conctn_anno_list=[]
    montage_biosemi=mne.channels.read_montage(montage_fname)

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
            if '254' not in events_coding.keys():
                raw.annotations.append(onset=0,duration=0,description='254')
                events = mne.events_from_annotations(raw)
                events_coding=events[1]
                events=np.asarray(events[0])  
            if '255' not in events_coding.keys():
                raw.annotations.append(onset=events[-1][0]/512,duration=0,description='255')
                events = mne.events_from_annotations(raw)
                events_coding=events[1]
                events=np.asarray(events[0])  

            events_code_start = events_coding['254']
            start = events[events[:,2]==events_code_start][0][0]
            events_code_end = events_coding['255']
            stop = events[events[:,2]==events_code_end][0][0]

            raw_cut_filt = raw.copy()
            raw_cut_filt.crop(tmin = start/raw.info['sfreq'], tmax = stop/raw.info['sfreq'])
            raw_cut_filt.notch_filter(freqs=power_freq_array)
            raw_cut_filt.filter(l_freq=1,h_freq=100)

            ############ annotation engineering ################
            index_dlt=0
            for i in range(raw_cut_filt.annotations.__len__()):
                if (raw_cut_filt.annotations.__getitem__(i-index_dlt)['description']) not in ['131','132','255']:
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
        peak_locs, peak_eeg = mne.preprocessing.peak_finder(VEOG_data, thresh = 110e-6)
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
        events_dict=removeItem_from_dict(events_dict,'255')

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
        epochs_ASR_clean =   mne.EpochsArray(DataClean,info=epochs_full.info,events=events[0],event_id = events[1])
        epochs_ASR_clean.add_channels([epochs_eog_raw_filt])
        if subj in bad_channels.keys() and session in bad_channels[subj].keys() :
            epochs_ASR_clean.add_channels([epochs_bad_channels])
            epochs_ASR_clean.interpolate_bads()
        return epochs_ASR_clean

    
# normal ones 110 75 (48subjs)
# subjs=['02', '04', '07', '11', '12', '14', '16', '18', '19', '21', '22', '26', '28', '30',
#        '32', '34', '36', '37', '38', '40', '42', '50', '51', '52', '53', '54', '55', '56', '58', '59',
#        '60', '63', '65', '67', '68', '70', '72', '73', '78', '83', '87', '88', '90', '91', '93', '94', '95', '96']

#subject acquire specific parameters 130 75(10subjs)
#SUBJECTS=(10 25 29 39 57 64 69 80 81 82)

#subject acquire specific parameters 140 752(2subjs)
#SUBJECTS=(35 71)

#subject acquire specific parameters 160 75(1subj)
#SUBJECTS=(79)

# subjs=['02']
# sessions=['1']
# for subj in subjs:
#     for session in sessions:

############################################### change here for new recordings ################################
subj = '34'
session = '1'
############################################### change here for new recordings ################################

preProc_ica_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/preProc_ica/'

ica_captions=[]
ica_figs=[]
preProc_ica_fname = preProc_ica_path+'subj0'+subj+'session'+session+'preProc_ica.fif'
if os.path.isfile(preProc_ica_fname)==False:
    print('no such ica file')

ica = read_ica(preProc_ica_fname)
epochs_ASR_clean = get_epochs_ASR_clean2seeICA(subj,session)

title = 'Sources related to %s artifacts (red)'
eog_inds, scores = ica.find_bads_eog(epochs_ASR_clean,ch_name='VEOG')
ica.plot_scores(scores, exclude=eog_inds, title=title % 'Veog', labels='Veog',show = False)
ica_figs.append(ica.plot_scores(scores, exclude=eog_inds, title=title % 'Veog', labels='Veog',show = False))
ica_captions.append('VEOG component correlation plot')

if eog_inds ==[]:
    print("no Veog component is found")

else:
    ica.plot_components(eog_inds, title=title % 'Veog', colorbar=True, show = False)
    ica_figs.append(ica.plot_components(eog_inds, title=title % 'Veog', colorbar=True, show = False))
    ica_captions.append('VEOG component plot')


eog_inds, scores = ica.find_bads_eog(epochs_ASR_clean,ch_name='HEOG')
ica.plot_scores(scores, exclude=eog_inds, title=title % 'Heog', labels='Heog', show= False)
ica_figs.append(ica.plot_scores(scores, exclude=eog_inds, title=title % 'Heog', labels='Heog', show= False))
ica_captions.append('HEOG component correlation plot')
#ica.plot_sources(epochs.average(), exclude=eog_inds, title=title % 'Veog')
if eog_inds ==[]:
    print("no Heog component is found")
else:
    ica.plot_components(eog_inds, title=title % 'Heog', colorbar=True, show = False)
    ica_figs.append(ica.plot_components(eog_inds, title=title % 'Heog', colorbar=True, show = False))
    ica_captions.append('HEOG component correlation plot')


### if we want to have access to component properties, we need to redo ASR process to get inst = epochs_ASR_cleaned

ica_figs=ica_figs+ ica.plot_components(inst = epochs_ASR_clean)
for i in range(len(ica.plot_components(inst=epochs_ASR_clean))):
#     ica_figs.append(ica.plot_components(inst=epochs_ASR_clean,show=False)[i])
    ica_captions.append('ica_components figure'+str(i))
ica.get_sources(inst = epochs_ASR_clean).plot(picks='all', n_channels = 10, n_epochs = 10)
epochs_ASR_clean.plot(n_channels = 22, n_epochs = 10, scalings = 100e-6)

## if we want a quick check decomment following scripts 
# from mne.preprocessing import read_ica
# from utils_preProcessingWorkflowJuly05 import *

# %matplotlib qt
# plt.close('all')
# preProc_ica_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/preProc_ica/'
# subj = '02'
# session = '1'
# preProc_ica_fname = preProc_ica_path+'subj0'+subj+'session'+session+'preProc_ica.fif'

# ica = read_ica(preProc_ica_fname)
# ica.plot_components()
plt.show()
report_path = '/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/report/'
f_report = report_path + 'subj0' + subj+'ica report.h5'
with mne.open_report(f_report) as rep:
# rep = mne.Report(image_format = 'png', subject = 'subj0'+subj)
    rep.add_figs_to_section(figs=ica_figs, captions=ica_captions, section = 'ica plot session'+session,
                                image_format = 'svg')
    # plt.pause(180)
    # https://labeling.ucsd.edu/tutorial
    rep.save(f_report, open_browser = False, overwrite = True)
# To visualize the report
# h5ToHtml(f_report)

