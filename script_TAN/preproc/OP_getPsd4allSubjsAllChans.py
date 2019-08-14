# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 18:04:25 2019

@author: gansheng.tan
"""

### get psd matrix ###

## to save at the end 
# >>> from tempfile import TemporaryFile
# >>> outfile = TemporaryFile()
# >>> x = np.arange(10)
# >>> np.save(outfile, x)
# >>> _ = outfile.seek(0) # Only needed here to simulate closing & reopening file
# >>> np.load(outfile)
# array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
import numpy as np
subjs=['02', '04','07', '11', '12', '14', '16', '18', '19', '21', '22', '26', '28', '30',
       '32', '34', '36', '37', '38', '40', '42', '50', '51', '52', '53', '54', '55', '56',
       '58', '59','60', '63', '65', '67', '68', '70', '72', '73', '78', '83', '87', '88', 
       '90', '91', '93', '94', '95', '96','10','25','29','39','57','64','69','80','81','82',
       '35','71','79','76','77']

states_codes={'VD':['111.0','112.0'],
              'FA':['211.0','212.0'],
              'OP':['311.0','312.0']}

# create VD baseline
psdList_mean4epochs4allSubjs = []
for subj in subjs:
    precleaned_epochs_fname = precleaned_epochs_path + 'subj0'+subj+'full_epo.fif'
    precleaned_epochs = mne.read_epochs(precleaned_epochs_fname, preload=True)
    precleaned_epochs_FA = precleaned_epochs[states_codes['FA']]
    data = precleaned_epochs_FA.get_data(picks=['eeg'])
    psd4epochs = []
    for num_epochs in range(data.shape[0]):
        sf = precleaned_epochs_FA.info['sfreq']
        psd, freqs = psd_array_multitaper(data[num_epochs], sf, fmin = 1, fmax =100,
                          adaptive=True,normalization='full',verbose=0)
        psd= np.log10(psd*10e12)
        psd4epochs.append(psd)
    psd_mean4epochs = np.append([psd4epochs[0]],psd4epochs[1:],axis = 0).mean(axis=0)
    psdList_mean4epochs4allSubjs.append(psd_mean4epochs)
psd_final = np.append([psdList_mean4epochs4allSubjs[0]],psdList_mean4epochs4allSubjs[1:],axis = 0)

np.save('/home/gansheng.tan/process_mne/INSERM_EEG_Enrico_Proc/data_eeglab/full_epochs_data/freq_channels4allSubjsFA.npy', psd_final)