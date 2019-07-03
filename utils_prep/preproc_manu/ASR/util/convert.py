# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 17:47:16 2018

@author: Manu
"""
from mne import io
#class Convert():
################################################################################
def EEGLab2Mne(src_directory, src_fname, montage_fname, dest_directory):
    fname = src_directory + src_fname
    raw = io.eeglab.read_raw_eeglab(fname, montage=montage_fname, preload = True)
#    raw = io.eeglab.read_raw_eeglab(fname, eog=(), montage=montage_fname, preload = True)
    save_fname = dest_directory + src_fname[:-4] +"-raw.fif"
    raw.save(save_fname,overwrite=True)

################################################################################
def Vhdr2Mne(src_directory, src_fname, montage_fname, dest_directory):
    fname = src_directory + src_fname
    raw = io.read_raw_brainvision(fname,  montage=montage_fname, preload = True)
    save_fname = dest_directory + src_fname[:-5] +"-raw.fif"
    raw.save(save_fname,overwrite=True)
