#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 13:10:04 2018

@author: benjamin.ador
"""
import time
t0 = time.perf_counter()

from header import *
from preproc import load_preproc

from fnmatch import filter
from mne.report import Report
from mne.time_frequency import psd_welch, psd_multitaper
from typing import List, Union


#%% PSD TO DATAARRAY FUNCTION

def psd_batch(task: str, subjects: List[str]=None, states: List[str]=['RS', 'FA', 'OM'], n_blk={'RS': 1, 'FA': 2, 'OM': 2}, verbose=False):
    """
    
    """
    verb = mne.set_log_level(verbose, return_old_level=True)
    if not subjects:
        subjects = sorted(get_subjlist(task))
    
    state_blk = []
    for state in states:
        state_blk.extend([state + str(b+1) for b in range(n_blk[state])])
    
    for st,state in enumerate(states):
        logger.info(state)
        for su,sub in enumerate(tqdm(subjects)):
            blocks = get_blocks(sub, task=task, state=state)
            for b,blk in enumerate(blocks):
                raw = load_preproc(task, sub, state, blk, exclude_eog=True, exclude_ecg=True, ICA_kwargs={'method': 'fastica'})
                psds, freqs = psd_welch(raw, n_fft=2400, fmax=int(raw.info['sfreq']/4), n_jobs=4) #2**13
                if not st and not su and not b:
                    channels = [ch.split('-')[0] for c,ch in enumerate(raw.ch_names) if c in mne.pick_types(raw.info, ref_meg=False)]
                    PSD = xr.DataArray(np.zeros((len(state_blk), len(subjects), len(channels), freqs.size)), 
                                       dims=['state', 'subject', 'chan', 'freq'], 
                                       coords={'state':state_blk, 'subject':subjects, 'chan':channels, 'freq':freqs})
                PSD.loc[state+str(b+1), sub] = psds
    
    PSD.to_netcdf(path=op.join(Analysis_path, task, 'meg', 'Alpha', 'PSD.nc'))
    mne.set_log_level(verb)
    return PSD


#%% COMPUTE OR LOAD PSD

task = 'SMEG'
subjects = get_subjlist(task)
subjects.sort()

#warnings.filterwarnings("ignore",category=DeprecationWarning)
#PSD = psd_average(task, subjects)

PSD = xr.open_dataarray(op.join(Analysis_path, task, 'meg', 'Alpha', 'PSD.nc'))
PSD.load()

PSD_norm = PSD/PSD.mean(['chan', 'freq'])


#%%

def plot_psd(data, ave_dim, new_fig=True):
    """
    
    """
    average = data.mean(ave_dim)
    sem = data.std(ave_dim)/np.sqrt(data[ave_dim].size)
    
    if new_fig:
        fig = plt.figure()
    plt.semilogx(data.freq, average)
    plt.fill_between(data.freq, average+sem, average-sem, alpha=.4)
    
    if new_fig:
        return fig


#%%

subjects = {'all': sorted(get_subjlist(task))}
no_blk2 = ['002', '004', '007', '016'] #FA2 and OM2 don't exist for these subjects (PSD=0)

subjects['expert'] = list(); subjects['novice'] = list()
for sub in subjects['all']:
    if expertise(sub) is 'N':
        subjects['novice'].append(sub)
    elif expertise(sub) is 'E':
        subjects['expert'].append(sub)
    else:
        warnings.warn('Expertise of subject {} is not specified, check masterfile.'.format(sub))

channels = {'all': PSD.chan.values}
channels['left'] = filter(channels['all'], 'ML*')
channels['left frontal'] = filter(channels['all'], 'MLF*')
channels['right'] = filter(channels['all'], 'MR*')
channels['central'] = filter(channels['all'], 'MZ[!O]*') #midline: all but occipital
channels['frontal'] = filter(channels['all'], 'M*F*')
channels['occipital'] = filter(channels['all'], 'M*O*')


#%% SELECTION PARAMETERS

state = 'FA1'

sub_key = 'all'
sub = '004'

chan_key = 'all' #select a subset of channels

fmin = 1 #PSD.freq.values[0]
fmax = PSD.freq.values[-1]


#%% SINGLE SUBJECT

#data = PSD.loc[state, sub, channels[chan_key], fmin:fmax]
#plot_psd(data, 'chan')
#plt.title('Average PSD of subject {sub} for state {s}\nAverage of {ch} channels'.format(s=state, sub=sub, ch=chan_key))


#%% GRAND AVERAGE

#data = PSD.loc[state, subjects[sub_key], channels[chan_key], fmin:fmax].mean('chan')
#plot_psd(data, 'subject')
#plt.title('Average PSD over {sub} subjects for state {s}\nAverage of {ch} channels'.format(s=state, sub=sub_key, ch=chan_key))
#plt.savefig(op.join(Analysis_path, task, 'meg', 'Plots', 'PSD', '{}-{}_subs.png'.format(state, sub_key)))


#%% PLOT BY GROUP

#plt.figure()
#keys = ['novice', 'expert']
#for sub_key in keys:
#    data = PSD.loc[state, subjects[sub_key], channels[chan_key], fmin:fmax].mean('chan')
#    plot_psd(data, 'subject', new_fig=False)
#plt.title('Average PSD for state {s}\nAverage of {ch} channels'.format(s=state, sub=sub_key, ch=chan_key))
#plt.legend(keys)
#plt.savefig(op.join(Analysis_path, task, 'meg', 'Plots', 'PSD', '{}-{}.png'.format(state, '+'.join(keys))))


#%% PLOT BY STATE

#plt.figure()
#keys = ['novice', 'expert']average = data.mean('subject')
#sub_key = 'all'
#for stt in PSD.state.values:
##for sub_key in keys:
#    data = PSD.loc[stt, subjects[sub_key], channels[chan_key], fmin:fmax].mean('chan')
#    average = data.mean('subject')
#    sem = data.std('subject')/sqrt(data.subject.size)
#    plt.semilogx(data.freq, average)
#    plt.fill_between(data.freq, average+sem, average-sem, alpha=.4)
#plt.title('Average PSD over {sub} subjects\nAverage of {ch} channels'.format(s=state, sub=sub_key, ch=chan_key))
#plt.legend(PSD.state.values)
#plt.savefig(op.join(Analysis_path, task, 'meg', 'Plots', 'PSD', '{}-{}.png'.format(state, '+'.join(keys))))


#%% ALPHA PEAK

fmin = 8
fmax = 13

lateral = ['*', 'L', 'R', 'Z']
ante_post = ['*', 'F', 'T', 'C', 'P', 'O']

alpha_tsv = 'group\tsubject\tstate\tlateral\tante_post\tpeak_freq\tpeak_val\tpeak_norm\n'
alpha_fname = op.join(Analysis_path, task, 'meg', 'Alpha', 'alpha_peak_{}_{}.tsv'.format(fmin, fmax))

for sub in tqdm(PSD.subject.values):
    report = Report(subject=sub)
    
    for stt in PSD.state.values:
        if sub in no_blk2 and fnmatch.fnmatch(stt, '*2'):
            continue #Old subjects only did one session of meditation states
        
        fig, axes = plt.subplots(len(lateral)-1, len(ante_post), sharex=True, sharey='col', figsize=(25,10))
        
        for row,lat in enumerate(lateral):
            for col,a_p in enumerate(ante_post):
                chs = filter(PSD.chan.values.tolist(), 'M'+lat+a_p+'*')
                if lat is 'Z':
                    chs = filter(PSD.chan.values.tolist(), 'M'+lat+'[!O]*') #Exclude Occipital channels from midline
                
                psd_sel = PSD.loc[stt, sub, chs, fmin:fmax].mean('chan')
                f_peak = psd_sel.freq[psd_sel.argmax()].values
                
                plot_data = PSD_norm.loc[stt, sub, chs, fmin:fmax]
                
                peak_val = psd_sel[psd_sel.argmax()].values
                peak_norm = plot_data.mean('chan')[psd_sel.argmax()].values
                alpha_tsv += '{}\t{}\t{}\t{}\t{}\t{:.1f}\t{}\t{}\n'.format(expertise(sub), sub, stt, lat, a_p, f_peak, peak_val, peak_norm)
                
                if lat is 'Z':
                    figZ = plot_psd(plot_data, 'chan')
                    plt.axvline(f_peak, color='r')
                    plt.axvline(ref_peak, linestyle='--', color='k')
                    plt.title('Midline: alpha peak at {:.1f} Hz'.format(f_peak))
                    figZ.set_size_inches((5,3))
                    break #Only one set of channels for midline
                
                plt.sca(axes[row,col])
                plot_psd(plot_data, 'chan', new_fig=False)
                plt.axvline(f_peak, color='r')
                if lat is '*' and a_p is '*':
                    ref_peak = f_peak
                else:
                    plt.axvline(ref_peak, linestyle='--', color='k')
                plt.title('{}{} ROI: {:.1f} Hz'.format(lat, a_p, f_peak))
        
        fig.suptitle('Subject {}, {}'.format(sub, stt))
        fig.set_tight_layout(True)
#        fig.set_size_inches((25,10))
        
        report.add_figs_to_section([fig, figZ], ['{}: {}-{} Hz PSD per ROI'.format(stt,fmin,fmax), '{}: Midline'.format(stt)], section=stt)
    report.save(op.join(Analysis_path, task, 'meg', 'Alpha', '{}-alpha_report.html'.format(sub)), open_browser=False, overwrite=True)

with open(alpha_fname, 'w') as fid:
    fid.write(alpha_tsv)
    fid.close()


#%% LOAD STATS
task = 'SMEG'
tests = ['FA_vs_RS', 'OM_vs_RS', 'FA_vs_OM', 'FA_vs_RS+E', 'OM_vs_RS+E', 'FA_vs_OM+E', 'FA_vs_RS+N', 'OM_vs_RS+N', 'FA_vs_OM+N', 'N_vs_E+FA-RS', 'N_vs_E+OM-RS', 'N_vs_E+FA-OM']
f1, f2 = .5, 100
stats = dict()
NS = []
for t in tests:
    file = op.join(Analysis_path, task, 'meg', 'Stats', 'PSD', '{}-{}Hz_{}_stat-ave.fif'.format(f1, f2, t))
    if op.isfile(file):
        stats[t] = mne.read_evokeds(file)
    else:
        NS.append(t)


#%% PLOT EVOKED
        
k = 'FA_vs_RS+E'
fmin = 3
fmax = 8
evoked = stats[k][-1]
#freqs = np.linspace(fmin,fmax,20)
freqs = evoked.times[np.logical_and(fmin <= evoked.times, evoked.times <= fmax)]
freqs = freqs[::len(freqs)//6 + 1 if len(freqs)%20 else 0]
#f = evoked.time_as_index(fmin)[0]
#freqs = evoked.times[f:f+20]
os.makedirs(op.join(Analysis_path, task, 'meg', 'Plots', 'PSD'), exist_ok=True)

for k in stats.keys():
    logger.info(k)
    evoked = stats[k][-1]
    evoked.plot_joint(times='peaks', title=k+' - {} significant clusters'.format(len(stats[k])-1),
                      topomap_args={'time_format': '%0.1f Hz', 'scalings': {'mag':1}, 'units': 'T value'},
                      ts_args={'scalings': {'mag':1}, 'units': {'mag':'T value'}})
    plt.savefig(op.join(Analysis_path, task, 'meg', 'Plots', 'PSD', 'joint_{}.svg'.format(k)), dpi='figure', transparent=True)
    plt.close()
    
    fig, axes = plt.subplots(5,7)
    bins = []
    locs = []
    for n, (fmin, fmax) in enumerate(zip([.5, 4, 10, 17, 30], [3, 9, 15, 27, 100])):
#    for n, (fmin, fmax) in enumerate(zip([.5, 3, 8, 13, 30], [2.5, 7, 12, 25, 100])):
        freqs = evoked.times[np.logical_and(np.logical_and(fmin <= evoked.times, evoked.times <= fmax), evoked.data.mean(0) != 0)]
        if freqs.size:
            freqs = freqs[::len(freqs)//6 if len(freqs)//6 else 1][:6]
        bins.extend(freqs)
        locs.extend(axes[n,:freqs.size])
    evoked.plot_topomap(times=bins, time_unit='s', time_format='%0.1f Hz', scalings={'mag':1}, units='T value', title=k, axes=locs)#, colorbar=False)#not n)
    fig.set_size_inches(15,15)
    fig.savefig(op.join(Analysis_path, task, 'meg', 'Plots', 'PSD', 'topo_{}.svg'.format(k)), dpi='figure', transparent=True)
    plt.close(fig)

clu_range = dict()
clu_range_i = dict()
for c,clu in enumerate(stats[k][:-1]):
    start = np.min(np.where(clu.data)[1])
    stop = np.max(np.where(clu.data)[1])
    clu_range[c] = (np.round(clu.times[start], 1), np.round(clu.times[stop], 1))
    clu_range_i[c] = (start, stop)

c = 0
cluster = stats[k][c].copy()
cluster.crop(*clu_range[c])
cluster.plot_topomap(time_unit = 's', title = cluster.comment)


#%% RUNNING TIME

t1 = time.perf_counter()
T = t1 - t0
print(colored(time.strftime('Finished %c',time.localtime()),'blue'))
print(colored('Elapsed time: {d}d {h}h {m}min {s}s'.format(s=round(T%60), m=round((T - T%60)%(60*60)/60), h=round((T - T%(60*60))%(24*60*60)/(60*60)), d=round((T - T%(24*60*60))/(24*60*60))), 'green'))