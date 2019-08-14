### get bp matrice for 6 baseline  expert/novice  VD FA OP -baseline

import pickle
import mne
import matplotlib.pyplot as plt
from mne.time_frequency import psd_array_multitaper
from scipy.integrate import simps
import numpy as np
import seaborn as sns
import pandas as pd
import sys
from mne.channels import find_ch_connectivity
from scipy.stats.distributions import f,t
from mpl_toolkits.axes_grid1 import make_axes_locatable
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

def getBpAbs4allChannels(epochs,rhythm):
    wavebands = {'alpha':[8,12],'theta':[3,7],'beta':[13,24],'lowG':[25,40],'highG':[60,90]}
    if rhythm in wavebands.keys():
        low,high =  wavebands[rhythm]
    else:
        print('not such rhythm')
    bpAbs_4Epochs=[]
    data = epochs.get_data(picks=['eeg'])
    for num_epochs in range(data.shape[0]):
        sf = epochs.info['sfreq']
        bpAbs_4allchannels = []
        psd, freqs = psd_array_multitaper(data[num_epochs], sf, fmin = 1, fmax =100,
                          adaptive=True,normalization='full',verbose=0)
        psd= np.log10(psd*10e12)
        freq_res = freqs[1] - freqs[0]
        idx_band = np.logical_and(freqs >= low, freqs <= high)
        bp_abs = simps(psd[:,idx_band], dx=freq_res)
        bpAbs_4Epochs.append(bp_abs)
    bpAbs_mean4Epochs = np.append([bpAbs_4Epochs[0]],bpAbs_4Epochs[1:],axis = 0).mean(axis=0)
    return bpAbs_mean4Epochs

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

listNovices = ['02', '04', '07', '10', '11', '12', '14', '16', '18', '19', '21', '22', '26', 
 '28', '29', '30', '32', '34', '35', '36', '37', '38','39', '40', '42', '81', '82', 
 '83', '87', '88', '90', '91', '93', '94', '95', '96']
listExperts = ['25', '50','51', '52',
 '53', '54', '55', '56', '57', '58', '59', '60', '63', '64', '65', '67','68', '69' ,'70' ,
               '71', '72', '73', '76', '77', '78' ,'79', '80']

group = sys.argv[1]
# group = 'expert'
bpAbs_mean4Epochs_VD4allsubjs = np.array([])
bpAbs_mean4Epochs_FA4allsubjs = np.array([])
bpAbs_mean4Epochs_OP4allsubjs = np.array([])
if group == 'expert':
    subjs = listExperts
elif group == 'novice':
    subjs = listNovices
else:
    print('no such group')
for subj in subjs:
    precleaned_epochs_fname = precleaned_epochs_path + 'subj0'+subj+'full_epo.fif'
    precleaned_epochs = mne.read_epochs(precleaned_epochs_fname, preload=True)
    
#     precleaned_epochs_VD = precleaned_epochs[states_codes['VD']]
#     precleaned_epochs_OP = precleaned_epochs[states_codes['OP']]
#     bpAbs_mean4Epochs_VD = getBpAbs4allChannels(precleaned_epochs_VD,'alpha')
#     bpAbs_mean4Epochs_OP= getBpAbs4allChannels(precleaned_epochs_OP,'alpha')

    precleaned_epochs_VD = precleaned_epochs[states_codes['VD']]
    precleaned_epochs_FA = precleaned_epochs[states_codes['FA']]
    bpAbs_mean4Epochs_VD = getBpAbs4allChannels(precleaned_epochs_VD,'alpha')
    bpAbs_mean4Epochs_FA= getBpAbs4allChannels(precleaned_epochs_FA,'alpha')

#     precleaned_epochs_FA = precleaned_epochs[states_codes['FA']]
#     precleaned_epochs_OP = precleaned_epochs[states_codes['OP']]
#     bpAbs_mean4Epochs_FA = getBpAbs4allChannels(precleaned_epochs_FA,'alpha')
#     bpAbs_mean4Epochs_OP= getBpAbs4allChannels(precleaned_epochs_OP,'alpha')
    
    if len(bpAbs_mean4Epochs_VD4allsubjs)==0:
        bpAbs_mean4Epochs_VD4allsubjs = bpAbs_mean4Epochs_VD
#         bpRelative_mean4Epochs_VD4allsubjs = bpRelative_mean4Epochs_VD
    else:
        bpAbs_mean4Epochs_VD4allsubjs = np.vstack((bpAbs_mean4Epochs_VD4allsubjs,bpAbs_mean4Epochs_VD))
#         bpRelative_mean4Epochs_VD4allsubjs = np.vstack((bpRelative_mean4Epochs_VD4allsubjs,
#                                                         bpRelative_mean4Epochs_VD))
        
#     if len(bpAbs_mean4Epochs_OP4allsubjs)==0:
#         bpAbs_mean4Epochs_OP4allsubjs = bpAbs_mean4Epochs_OP
# #         bpRelative_mean4Epochs_OP4allsubjs = bpRelative_mean4Epochs_OP
#     else:
#         bpAbs_mean4Epochs_OP4allsubjs = np.vstack((bpAbs_mean4Epochs_OP4allsubjs,bpAbs_mean4Epochs_OP))
# #         bpRelative_mean4Epochs_OP4allsubjs = np.vstack((bpRelative_mean4Epochs_OP4allsubjs,
# #                                                         bpRelative_mean4Epochs_OP))

    if len(bpAbs_mean4Epochs_FA4allsubjs)==0:
        bpAbs_mean4Epochs_FA4allsubjs = bpAbs_mean4Epochs_FA
#         bpRelative_mean4Epochs_OP4allsubjs = bpRelative_mean4Epochs_OP
    else:
        bpAbs_mean4Epochs_FA4allsubjs = np.vstack((bpAbs_mean4Epochs_FA4allsubjs,bpAbs_mean4Epochs_FA))
#         bpRelative_mean4Epochs_OP4allsubjs = np.vstack((bpRelative_mean4Epochs_OP4allsubjs,
#                                                         bpRelative_mean4Epochs_OP))
        
# bpAbs_mean4Epochs2test = [np.expand_dims(bpAbs_mean4Epochs_VD4allsubjs,axis=1),
#                           np.expand_dims(bpAbs_mean4Epochs_OP4allsubjs,axis=1)]

bpAbs_mean4Epochs2test = [np.expand_dims(bpAbs_mean4Epochs_VD4allsubjs,axis=1),
                          np.expand_dims(bpAbs_mean4Epochs_FA4allsubjs,axis=1)]

# bpAbs_mean4Epochs2test = [np.expand_dims(bpAbs_mean4Epochs_FA4allsubjs,axis=1),
#                           np.expand_dims(bpAbs_mean4Epochs_OP4allsubjs,axis=1)]


# with open('/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/statistic/VDOP_alpha_baseline_Abs'+group+'.txt', "wb") as fp:   #Pickling
#     pickle.dump(bpAbs_mean4Epochs2test, fp)
    
with open('/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/statistic/VDFA_alpha_baseline_Abs'+group+'.txt', "wb") as fp:   #Pickling
    pickle.dump(bpAbs_mean4Epochs2test, fp)
# with open('/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/statistic/FAOP_alpha_baseline_Abs'+group+'.txt', "wb") as fp:   #Pickling
#     pickle.dump(bpAbs_mean4Epochs2test, fp)
    
precleaned_epochs_fname = precleaned_epochs_path + 'subj004full_epo.fif'
precleaned_epochs = mne.read_epochs(precleaned_epochs_fname, preload=True)
connectivity, ch_names = find_ch_connectivity(precleaned_epochs.info, ch_type='eeg')

p_threshold = 0.05
# threshold = -t.ppf(p_threshold/2,np.expand_dims(bpAbs_mean4Epochs_VD4allsubjs,axis=1).shape[0]-1)
# cluster_stats = mne.stats.spatio_temporal_cluster_1samp_test(np.expand_dims(bpAbs_mean4Epochs_OP4allsubjs,axis=1)-np.expand_dims(bpAbs_mean4Epochs_VD4allsubjs,axis=1), 
#                                                              n_permutations=10000,tail=0,threshold=threshold,
#                                              n_jobs=2, buffer_size=None,verbose=True,
#                                              connectivity=connectivity)

threshold = -t.ppf(p_threshold/2,np.expand_dims(bpAbs_mean4Epochs_VD4allsubjs,axis=1).shape[0]-1)
cluster_stats = mne.stats.spatio_temporal_cluster_1samp_test(np.expand_dims(bpAbs_mean4Epochs_FA4allsubjs,axis=1)-np.expand_dims(bpAbs_mean4Epochs_VD4allsubjs,axis=1), 
                                                             n_permutations=10000,tail=0,threshold=threshold,
                                             n_jobs=2, buffer_size=None,verbose=True,
                                             connectivity=connectivity)

# threshold = -t.ppf(p_threshold/2,np.expand_dims(bpAbs_mean4Epochs_FA4allsubjs,axis=1).shape[0]-1)
# cluster_stats = mne.stats.spatio_temporal_cluster_1samp_test(np.expand_dims(bpAbs_mean4Epochs_OP4allsubjs,axis=1)-np.expand_dims(bpAbs_mean4Epochs_FA4allsubjs,axis=1), 
#                                                              n_permutations=10000,tail=0,threshold=threshold,
#                                              n_jobs=2, buffer_size=None,verbose=True,
#                                              connectivity=connectivity)

T_obs, clusters, p_values, _ = cluster_stats
print(clusters)
# good_cluster_inds = np.array(range(len(clusters)))
# precleaned_epochs_fname = precleaned_epochs_path + 'subj004full_epo.fif'
# precleaned_epochs = mne.read_epochs(precleaned_epochs_fname, preload=True)
# pos = mne.find_layout(precleaned_epochs.info).pos
# for i_clu, clu_idx in enumerate(good_cluster_inds):
#     # unpack cluster information, get unique indices
#     time_inds, space_inds = np.squeeze(clusters[clu_idx])
#     ch_inds = np.unique(space_inds)
#     time_inds = np.unique(time_inds)

#     # get topography for bp-mean
# #     bp_map = np.squeeze((np.expand_dims(bpAbs_mean4Epochs_OP4allsubjs,axis=1)-np.expand_dims(bpAbs_mean4Epochs_VD4allsubjs,axis=1)).mean(axis=0))
# #     bp_map = np.squeeze((np.expand_dims(bpAbs_mean4Epochs_FA4allsubjs,axis=1)-np.expand_dims(bpAbs_mean4Epochs_VD4allsubjs,axis=1)).mean(axis=0))
#     bp_map = np.squeeze((np.expand_dims(bpAbs_mean4Epochs_OP4allsubjs,axis=1)-np.expand_dims(bpAbs_mean4Epochs_FA4allsubjs,axis=1)).mean(axis=0))

    
#     # create spatial mask
#     mask = np.zeros((bp_map.shape[0], 1), dtype=bool)
#     mask[ch_inds, :] = True

#     # initialize figure
#     fig, ax_topo = plt.subplots(1, 1, figsize=(10, 3))

#     # plot average test statistic and mark significant sensors
#     image, _ = mne.viz.plot_topomap(bp_map, pos, mask=mask, axes=ax_topo, cmap='Reds',
#                             vmin=np.min, vmax=np.max, show=False)
#     divider = make_axes_locatable(ax_topo)

#     # add axes for colorbar
#     ax_colorbar = divider.append_axes('right', size='5%', pad=0.05)
#     plt.colorbar(image, cax=ax_colorbar)
# #     ax_topo.set_xlabel('Averaged baseline alpha bandpower OP-VD for {}'.format(group))
    
# #     ax_topo.set_xlabel('Averaged baseline alpha bandpower FA-VD for {}'.format(group))

#     ax_topo.set_xlabel('Averaged baseline alpha bandpower OP-FA for {}'.format(group))
    
#     mne.viz.tight_layout(fig=fig)
#     fig.subplots_adjust(bottom=.05)
# #     plt.show()
# #     fig.savefig('/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/statistic/OP-VD_alpha_baseline topoplot.png')
# #     fig.savefig('/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/statistic/FA-VD_alpha_baseline topoplot'+'i_clu'+'.png')
#     fig.savefig('/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/statistic/OP-FA_alpha_baseline topoplot'+'i_clu'+'.png')

