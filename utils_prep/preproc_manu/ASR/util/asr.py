# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 08:28:09 2019

@author: Manu
"""
import numpy as np
import scipy
from util import tools
from scipy import signal
from scipy import linalg
from numpy import matlib

def clean_windows(Signal,srate,max_bad_channels,zthresholds,window_len):
# Remove periods with abnormally high-power content from continuous data.
# [Signal,Mask] = clean_windows(Signal,MaxBadChannels,PowerTolerances,WindowLength,WindowOverlap,MaxDropoutFraction,Min)
#
# This function cuts segments from the data which contain high-power artifacts. Specifically,
# only windows are retained which have less than a certain fraction of "bad" channels, where a channel
# is bad in a window if its power is above or below a given upper/lower threshold (in standard 
# deviations from a robust estimate of the EEG power distribution in the channel).
#
# In:
#   Signal         : Continuous data set, assumed to be appropriately high-passed (e.g. >1Hz or
#                    0.5Hz - 2.0Hz transition band)
#
#   MaxBadChannels : The maximum number or fraction of bad channels that a retained window may still
#                    contain (more than this and it is removed). Reasonable range is 0.05 (very clean
#                    output) to 0.3 (very lax cleaning of only coarse artifacts). Default: 0.2.
#
#   PowerTolerances: The minimum and maximum standard deviations within which the power of a channel
#                    must lie (relative to a robust estimate of the clean EEG power distribution in 
#                    the channel) for it to be considered "not bad". Default: [-3.5 5].
#
#
#   The following are detail parameters that usually do not have to be tuned. If you can't get
#   the function to do what you want, you might consider adapting these to your data.
#
#   WindowLength    : Window length that is used to check the data for artifact content. This is 
#                     ideally as long as the expected time scale of the artifacts but not shorter 
#                     than half a cycle of the high-pass filter that was used. Default: 1.
#
#   WindowOverlap : Window overlap fraction. The fraction of two successive windows that overlaps.
#                   Higher overlap ensures that fewer artifact portions are going to be missed (but
#                   is slower). (default: 0.66)
# 
#   MaxDropoutFraction : Maximum fraction that can have dropouts. This is the maximum fraction of
#                        time windows that may have arbitrarily low amplitude (e.g., due to the
#                        sensors being unplugged). (default: 0.1)
#
#   MinCleanFraction : Minimum fraction that needs to be clean. This is the minimum fraction of time
#                      windows that need to contain essentially uncontaminated EEG. (default: 0.25)
#
#   
#   The following are expert-level parameters that you should not tune unless you fully understand
#   how the method works.
#
#   TruncateQuantile : Truncated Gaussian quantile. Quantile range [upper,lower] of the truncated
#                      Gaussian distribution that shall be fit to the EEG contents. (default: [0.022 0.6])
#
#   StepSizes : Grid search stepping. Step size of the grid search, in quantiles; separately for
#               [lower,upper] edge of the truncated Gaussian. The lower edge has finer stepping
#               because the clean data density is assumed to be lower there, so small changes in
#               quantile amount to large changes in data space. (default: [0.01 0.01])
#
#   ShapeRange : Shape parameter range. Search range for the shape parameter of the generalized
#                Gaussian distribution used to fit clean EEG. (default: 1.7:0.15:3.5)
#
# Out:
#   SignalClean : data set with bad time periods removed.
#
#   Mask   : mask of retained samples (logical array)


    window_overlap = 0.66
    max_dropout_fraction = 0.1
    min_clean_fraction = 0.25
    truncate_quant = [0.0220,0.6000]
    step_sizes = [0.01,0.01]
    shape_range = np.linspace(1.7,3.5,13)
    max_bad_channels = np.round(Signal.shape[0]*max_bad_channels);
    
#    Signal = Signal *1e6
    [C,S] = Signal.shape;
    N = int(window_len*srate);
    wnd = np.arange(0,N);
    offsets = np.int_(np.arange(0,S-N,np.round(N*(1-window_overlap))))

    print('Determining time window rejection thresholds...')
    print('for each channel...')
    
    
    wz=np.array([])
    for ichan in range(C):
        X = Signal[ichan,:]**2
        Y=[]
        for joffset in offsets:
            Y.append(np.sqrt(np.sum(X[joffset:joffset+N])/N))
            
        Y=np.transpose(Y)
        mu,sig,alpha,beta = tools.fit_eeg_distribution(Y, min_clean_fraction, max_dropout_fraction,truncate_quant, step_sizes,shape_range)
        if (ichan==0):
           wz = (Y-mu)/sig
        else:
            wz=np.vstack((wz,np.array((Y-mu)/sig)))

    # sort z scores into quantiles
    swz = np.sort(wz,axis=0)
    
    #    determine which windows to remove
    if (np.max(zthresholds)>0):
        remove_mask1 = swz[-(np.int(max_bad_channels)+1),:] > np.max(zthresholds)
              
    if (np.min(zthresholds)<0):
        remove_mask2 = swz[1+np.int(max_bad_channels-1),:] < np.min(zthresholds)
        
    remove_mask=np.logical_or(remove_mask1, remove_mask2)
    removed_windows = np.where(remove_mask)
    
    sample_maskidx = []
    for iremoved in range(len(removed_windows[0])):
        if (iremoved==0):
            sample_maskidx=np.arange(offsets[removed_windows[0][iremoved]],offsets[removed_windows[0][iremoved]]+N)
        else:
            sample_maskidx=np.vstack((sample_maskidx,(np.arange(offsets[removed_windows[0][iremoved]],offsets[removed_windows[0][iremoved]]+N))))

    sample_mask2remove = np.unique(sample_maskidx)
    
    SignalClean = np.delete(Signal,sample_mask2remove,1)
    sample_mask = np.ones((1, S), dtype=bool)
    sample_mask[0,sample_mask2remove]=False
    
    return SignalClean,sample_mask




def YW_filter(Data,srate,iirstate_in):
#     FilterB, FilterA : Coefficients of an IIR filter that is used to shape the spectrum of the signal
#                      when calculating artifact statistics. The output signal does not go through
#                      this filter. This is an optional way to tune the sensitivity of the algorithm
#                      to each frequency component of the signal. The default filter is less
#                      sensitive at alpha and beta frequencies and more sensitive at delta (blinks)
#                      and gamma (muscle) frequencies. Default: 
#                      [b,a] = yulewalk(8,[[0 2 3 13 16 40 min(80,srate/2-1)]*2/srate 1],[3 0.75 0.33 0.33 1 1 3 3]);
    [C,S] = Data.shape
    F=np.array([0,2,3,13,16,40,np.minimum(80.0,(srate/2.0)-1.0),srate/2.0])*2.0/srate
    M = np.array([3,0.75,0.33,0.33,1,1,3,3])
    B,A = tools.yulewalk(8,F,M)
    
    # apply the signal shaping filter and initialize the IIR filter state
    DataFilt = np.zeros((C,S))
    iirstate = np.zeros((C,len(A)-1))
    zi = signal.lfilter_zi(B, A)
    for ichan in range(C):
        if (iirstate_in is None):
#            DataFilt[ichan,:], iirstate[ichan,:] = signal.lfilter(B,A,Data[ichan,:],zi=zi*0)#zi*Data[ichan,0])
            DataFilt[ichan,:], iirstate[ichan,:] = signal.lfilter(B,A,Data[ichan,:],zi=zi*Data[ichan,0])
        else:
            DataFilt[ichan,:], iirstate[ichan,:] = signal.lfilter(B,A,Data[ichan,:],zi=iirstate_in[ichan,:])

            


    return DataFilt, iirstate









def asr_calibrate(Data,srate,cutoff):
    # Calibration function for the Artifact Subspace Reconstruction (ASR) method.
    # State = asr_calibrate(Data,SamplingRate,Cutoff,BlockSize,FilterB,FilterA,WindowLength,WindowOverlap,MaxDropoutFraction,MinCleanFraction)
    #
    # The input to this data is a multi-channel time series of calibration data. In typical uses the
    # calibration data is clean resting EEG data of ca. 1 minute duration (can also be longer). One can
    # also use on-task data if the fraction of artifact content is below the breakdown point of the
    # robust statistics used for estimation (50# theoretical, ~30# practical). If the data has a
    # proportion of more than 30-50# artifacts then bad time windows should be removed beforehand. This
    # data is used to estimate the thresholds that are used by the ASR processing function to identify
    # and remove artifact components.
    #
    # The calibration data must have been recorded for the same cap design from which data for cleanup
    # will be recorded, and ideally should be from the same session and same subject, but it is possible
    # to reuse the calibration data from a previous session and montage to the extent that the cap is
    # placed in the same location (where loss in accuracy is more or less proportional to the mismatch
    # in cap placement).
    #
    # The calibration data should have been high-pass filtered (for example at 0.5Hz or 1Hz using a
    # Butterworth IIR filter).
    #
    # In:
    #   Data : Calibration data [#channels x #samples]; *zero-mean* (e.g., high-pass filtered) and
    #          reasonably clean EEG of not much less than 30 seconds length (this method is typically
    #          used with 1 minute or more).
    #
    #   SamplingRate : Sampling rate of the data, in Hz.
    #
    #
    #   The following are optional parameters (the key parameter of the method is the RejectionCutoff):
    #
    #   RejectionCutoff: Standard deviation cutoff for rejection. Data portions whose variance is larger
    #                    than this threshold relative to the calibration data are considered missing
    #                    data and will be removed. The most aggressive value that can be used without
    #                    losing too much EEG is 2.5. A quite conservative value would be 5. Default: 5.
    #
    #   Blocksize : Block size for calculating the robust data covariance and thresholds, in samples;
    #               allows to reduce the memory and time requirements of the robust estimators by this 
    #               factor (down to Channels x Channels x Samples x 16 / Blocksize bytes). Default: 10
    #
    #   FilterB, FilterA : Coefficients of an IIR filter that is used to shape the spectrum of the signal
    #                      when calculating artifact statistics. The output signal does not go through
    #                      this filter. This is an optional way to tune the sensitivity of the algorithm
    #                      to each frequency component of the signal. The default filter is less
    #                      sensitive at alpha and beta frequencies and more sensitive at delta (blinks)
    #                      and gamma (muscle) frequencies. Default: 
    #                      [b,a] = yulewalk(8,[[0 2 3 13 16 40 min(80,srate/2-1)]*2/srate 1],[3 0.75 0.33 0.33 1 1 3 3]);
    #
    #   WindowLength : Window length that is used to check the data for artifact content. This is 
    #                  ideally as long as the expected time scale of the artifacts but short enough to 
    #				   allow for several 1000 windows to compute statistics over. Default: 0.5.
    #
    #   WindowOverlap : Window overlap fraction. The fraction of two successive windows that overlaps.
    #                   Higher overlap ensures that fewer artifact portions are going to be missed (but
    #                   is slower). Default: 0.66
    #
    #   MaxDropoutFraction : Maximum fraction of windows that can be subject to signal dropouts 
    #                        (e.g., sensor unplugged), used for threshold estimation. Default: 0.1
    #
    #   MinCleanFraction : Minimum fraction of windows that need to be clean, used for threshold
    #                      estimation. Default: 0.25
    #
    #
    # Out:
    #   State : initial state struct for asr_process
    
    [C,S] = Data.shape
    blocksize = 10
    window_len = 0.5
    window_overlap = 0.66
    
    max_dropout_fraction = 0.1
    min_clean_fraction = 0.25
#    F=np.array([0,2,3,13,16,40,np.minimum(80.0,(srate/2.0)-1.0),srate/2.0])*2.0/srate
#    M = np.array([3,0.75,0.33,0.33,1,1,3,3])
#    B,A = tools.yulewalk(8,F,M)
#    
#    # apply the signal shaping filter and initialize the IIR filter state
#    SigFilt = np.zeros((C,S))
#    iirstate = np.zeros((C,len(A)-1))
#    zi = signal.lfilter_zi(B, A)
#    for ichan in range(C):
#        SigFilt[ichan,:], iirstate[ichan,:] = signal.lfilter(B,A,Data[ichan,:],zi=zi*0)#zi*Data[ichan,0])
    Data = Data.T
    U = np.zeros((len(np.arange(0,S,blocksize)),C*C))
    for k in range(blocksize):
        rangevect = np.minimum(S-1,np.arange(k,S+k,blocksize))
        Xrange = Data[rangevect,:]
        for ic in range(C):
            islice = np.arange((ic*C),((ic+1)*C),1,dtype=int)
            U[:,islice] = U[:,islice] + (Xrange*np.transpose(np.matlib.repmat(Xrange[:,ic],C,1)))
        

    # get the mixing matrix M
    M = scipy.linalg.sqrtm(np.real(np.reshape(tools.block_geometric_median(U/blocksize,1),(C,C))));
    
    # window length for calculating thresholds
    N = int(np.round(window_len*srate))
    
    
    
    # get the threshold matrix T
    print('Determining per-component thresholds...');
    
    D,Vtmp = scipy.linalg.eig(M)
    V=Vtmp[:,np.argsort(D)]
    
    X = np.abs(np.dot(Data,V));
    
    offsets = np.int_(np.arange(0,S-N,np.round(N*(1-window_overlap))))
    truncate_quant = [0.0220,0.6000]
    step_sizes = [0.01,0.01]
    shape_range = np.linspace(1.7,3.5,13)
    mu=np.zeros(C)
    sig=np.zeros(C)
    for ichan in range(C):
        rms = X[:,ichan]**2
        Y=[]
        for joffset in offsets:
            Y.append(np.sqrt(np.sum(rms[joffset:joffset+N])/N))
            
        Y=np.transpose(Y)
        mu[ichan],sig[ichan],alpha,beta = tools.fit_eeg_distribution(Y, min_clean_fraction, max_dropout_fraction,truncate_quant, step_sizes,shape_range)
    
    T = np.dot(np.diag(mu + cutoff*sig),V.T)
#    print('mu',mu)
#    print('sig',sig)
#    
    print('done.');
    calibASRparam= {'M':M,'T':T}
    return calibASRparam
            #'cov',[],'carry',[],'iir',iirstate,'last_R',[],'last_trivial',true}
    # initialize the remaining filter state
    #state = struct('M',M,'T',T,'B',B,'A',A,'cov',[],'carry',[],'iir',iirstate,'last_R',[],'last_trivial',true);
    
    
    
    
    
def asr_process_on_epoch(epoch2correct, epochYWfiltered,state):
# Processing function for the Artifact Subspace Reconstruction (ASR) method.
# EpochClean = asr_process_on_epoch(epoch2correct, epochYWfiltered,state)
#
# This function is used to clean multi-channel signal using the ASR method. The required inputs are 
# the data matrix, the sampling rate of the data, and the filter state (as initialized by
# asr_calibrate).         
    
    
    [C,S] = epochYWfiltered.shape
    epochYWfiltered = scipy.signal.detrend(epochYWfiltered, axis=1, type='constant')
    Xcov = np.cov(epochYWfiltered,bias=True)
    D,Vtmp = np.linalg.eig(Xcov)
    V=np.real(Vtmp[:,np.argsort(D)])
    D=np.real(D[np.argsort(D)])
    
    
    maxdims = int(np.fix(0.66*C))
    
    #determine which components to keep (variance below directional threshold or not admissible for rejection)
    keep=(D<np.sum(np.dot(state['T'],V)**2,axis=0)) +  ((np.arange(C))<(C-maxdims))
    
    
    trivial = keep.all()
    # update the reconstruction matrix R (reconstruct artifact components using the mixing matrix)
    if trivial:
        R = np.eye(C)
    else:
        VT = (np.dot(V.T,state['M']))
        demux = np.zeros((C,C))
        for icov in range(C):
            demux[icov,:] = VT[:,icov]*keep
        demux = np.transpose(demux)
        R = np.dot(np.dot(state['M'],np.linalg.pinv(demux)),V.T)
    
    EpochClean = np.dot(R,epoch2correct)

    return EpochClean