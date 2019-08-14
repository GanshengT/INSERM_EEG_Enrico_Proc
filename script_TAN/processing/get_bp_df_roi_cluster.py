# create bp for 2 ROI into csv - using panda dataframe
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

# all subject 63
subjs=['02', '04','07', '11', '12', '14', '16', '18', '19', '21', '22', '26', '28', '30',
       '32', '34', '36', '37', '38', '40', '42', '50', '51', '52', '53', '54', '55', '56',
       '58', '59','60', '63', '65', '67', '68', '70', '72', '73', '78', '83', '87', '88', 
       '90', '91', '93', '94', '95', '96','10','25','29','39','57','64','69','80','81','82',
       '35','71','79','76','77']


listNovices = ['02', '04', '07', '10', '11', '12', '14', '16', '18', '19', '21', '22', '26', 
 '28', '29', '30', '32', '34', '35', '36', '37', '38','39', '40', '42', '81', '82', 
 '83', '87', '88', '90', '91', '93', '94', '95', '96']
listExperts = ['25', '50', '51', '52',
 '53', '54', '55', '56', '57', '58', '59', '60', '63', '64', '65', '67','68', '69' ,'70' ,
               '71', '72', '73', '76', '77', '78' ,'79', '80']

ROIs = {'frontal':['Fz','F1','F3','AF3','AFz','AF4','F4','F2'],
       'occipital':['POz','Oz','O1','O2']}

CondStates={'baseline':{'VD':['111.0','112.0'],'FA':['211.0','212.0'],'OP':['311.0','312.0']},
           'safe':{'VD':['121.0','122.0'],'FA':['221.0','222.0'],'OP':['321.0','322.0']},
           'threat':{'VD':['131.0','132.0'],'FA':['231.0','232.0'],'OP':['331.0','332.0']}}

wavebands = {'alpha':[8,12],'theta':[3,7],'beta':[13,24],'lowG':[25,40],'highG':[60,90]}


frames_subject = []
for subj in subjs:
    precleaned_epochs_fname = precleaned_epochs_path + 'subj0'+subj+'full_epo.fif'
    precleaned_epochs = mne.read_epochs(precleaned_epochs_fname, preload=True)
    frames_condition = []
    for condition in CondStates.keys():
        frames_state=[]
        for state in CondStates[condition].keys():
#             precleaned_epochs_engineering = precleaned_epochs.copy()
            precleaned_epochs_selected = precleaned_epochs[CondStates[condition][state]]
            frames_ROI=[]
            for roi in ROIs.keys():
                precleaned_epochs_engineering = precleaned_epochs_selected.copy()
                data = precleaned_epochs_engineering.pick_channels(ROIs[roi])._data
                frames_waveband = []
                for waveband in wavebands.keys():
                    low,high = wavebands[waveband]
                    bpAbs_4Epochs=[]
                    bpRelative_4Epochs=[]
                    for num_epochs in range(data.shape[0]):
                        sf = precleaned_epochs_engineering.info['sfreq']
                        bpAbs_4allchannels = []
                        bpRelative_4allchannels = []
                        psd, freqs = psd_array_multitaper(data[num_epochs], sf, fmin = fmin, fmax =fmax,
                                          adaptive=True,normalization='full',verbose=0)
                        psd = np.log10(psd*10e12)
                        # average over channels
                        psd = psd.mean(axis=0)
                        freq_res = freqs[1] - freqs[0]
                        bp_total = simps(psd, dx=freq_res)
                        idx_band = np.logical_and(freqs >= low, freqs <= high)
                        bp_abs = simps(psd[idx_band], dx=freq_res)
                        bp_relative = bp_abs/bp_total
                        bpAbs_4Epochs.append(bp_abs)
                        bpRelative_4Epochs.append(bp_relative)
                    bpAbs_mean4Epochs = np.array(bpAbs_4Epochs).mean()
                    bpRelative_mean4Epochs = np.array(bpRelative_4Epochs).mean()
                        
                    if subj in listNovices:
                        group = ['novice']
                    elif subj in listExperts:
                        group = ['expert']
                    data_df_psd = {'subj':[subj],'ROI':[roi],'state':[state],'condition':[condition],
                                   'group':group,'rhythm':[waveband],'AbsBP':[bpAbs_mean4Epochs],
                                  'RelativeBP':[bpRelative_mean4Epochs]}
                    df_psd_perWaveband = pd.DataFrame(data_df_psd)
                    frames_waveband.append(df_psd_perWaveband)
                df_psd_perROI = pd.concat(frames_waveband,ignore_index=True)
                frames_ROI.append(df_psd_perROI)
            df_psd_perState = pd.concat(frames_ROI,ignore_index=True)
            frames_state.append(df_psd_perState)
        df_psd_perCondition = pd.concat(frames_state,ignore_index=True)
        frames_condition.append(df_psd_perCondition)
    df_psd_perSubject = pd.concat(frames_condition,ignore_index=True)
    frames_subject.append(df_psd_perSubject)
df_psd = pd.concat(frames_subject,ignore_index=True)
df_psd.to_csv(path_or_buf='/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/df_bp_log_cluster.csv',index = False)
print('finished')