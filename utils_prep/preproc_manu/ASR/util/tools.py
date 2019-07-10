# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 09:11:04 2018

@author: Manu
"""

import numpy as np
import mne
import scipy.special
from scipy.linalg import toeplitz
from util import tools
from scipy import signal
from numpy import linalg
import h5py



def RejectThresh(epochs,PercentageOfEpochsRejected):
    NbEpoch2Keep = np.fix(epochs.__len__() * (1.0-(PercentageOfEpochsRejected/100)))
    Mat_epoch =  epochs.get_data()
    MinWOI = Mat_epoch.min(axis=2)
    MaxWOI = Mat_epoch.max(axis=2)
    Peak2Peak = MaxWOI-MinWOI
    MaxPeak2PeakROI = Peak2Peak.max(axis=1)    
    MaxPeak2PeakROI.sort()
    Treshold = MaxPeak2PeakROI[int(NbEpoch2Keep)]
    
    return Treshold

def ComputeVarICAWeights(ica):
    fast_dot = np.dot
    VarIcaWeigths = np.zeros(ica.n_components_, dtype='float')
    for icomp in range(ica.n_components_):
        maps = fast_dot(ica.mixing_matrix_[:, icomp].T,ica.pca_components_[:ica.n_components_])
        mapsNorm=maps/(np.max(np.abs([maps.max(),maps.min()])))
        VarIcaWeigths[icomp] = np.var(mapsNorm)
        
    return VarIcaWeigths



def AddVirtualEogChannels(raw,ChanName4VEOG,ChanName4HEOG_l,ChanName4HEOG_r):
    FlagVEOG = False
    FlagHEOG_L = False
    FlagHEOG_R = False
    FlagHEOG = True 
    
    if ChanName4VEOG is not None:               
        rawSelecChan4Veog  = raw.copy().pick_channels(ChanName4VEOG)
        rawVEogData = np.zeros((1, rawSelecChan4Veog.n_times), dtype='float')
        rawVEogData[0,:] = (rawSelecChan4Veog.get_data(picks=range(len(ChanName4VEOG))).sum(axis=0))
        FlagVEOG = True

     # Create virtual Horizontal EOG
    if (ChanName4HEOG_l is not None) : 
        rawSelecChan4Heog_l = raw.copy().pick_channels(ChanName4HEOG_l)
        rawHEogL_Data = np.zeros((1, rawSelecChan4Heog_l.n_times), dtype='float')
        rawHEogL_Data[0,:] = (rawSelecChan4Heog_l.get_data(picks=range(len(ChanName4HEOG_l))).sum(axis=0))
        FlagHEOG_L = True
        
        
    if (ChanName4HEOG_r is not None):
        rawSelecChan4Heog_r = raw.copy().pick_channels(ChanName4HEOG_r)
        rawHEogR_Data = np.zeros((1, rawSelecChan4Heog_r.n_times), dtype='float')
        rawHEogR_Data[0,:] = (rawSelecChan4Heog_r.get_data(picks=range(len(ChanName4HEOG_r))).sum(axis=0))
        FlagHEOG_R = True
        
    if FlagHEOG_L:
        if FlagHEOG_R:
            rawHEogData = rawHEogL_Data - rawHEogR_Data
        else:
            rawHEogData = rawHEogL_Data
    else:
        if FlagHEOG_R:
            rawHEogData = rawHEogR_Data
        else:
            FlagHEOG = False
    
    rawWithVirtEOG = raw.copy()
    if FlagVEOG:
        infoVEog = mne.create_info(['VEOG'], rawWithVirtEOG.info['sfreq'], ['eeg'])
        VEogRawArray  = mne.io.RawArray(rawVEogData, infoVEog)
        rawWithVirtEOG.add_channels([VEogRawArray], force_update_info=True)
            
    if FlagHEOG:
        infoHEog = mne.create_info(['HEOG'], rawWithVirtEOG.info['sfreq'], ['eeg'])
        HEogRawArray  = mne.io.RawArray(rawHEogData, infoHEog)
        rawWithVirtEOG.add_channels([HEogRawArray], force_update_info=True)
    
    return rawWithVirtEOG





#Estimate the mean and standard deviation of clean EEG from contaminated data.
def fit_eeg_distribution(X,min_clean_fraction,max_dropout_fraction,quants,step_sizes,beta):
    
    # sort data so we can access quantiles directly  
    X = np.sort(X)
    n = len(X);
    
    # calc z bounds for the truncated standard generalized Gaussian pdf and pdf rescaler
    quants=np.array(quants)
    zbounds=[]
    rescale=[]
    for b in range(len(beta)):
        zbounds.append( np.sign(quants-1/2) * scipy.special.gammaincinv(1/beta[b],np.sign(quants-1/2)*(2*quants-1))**(1/beta[b]))
        rescale.append(beta[b]/(2*scipy.special.gamma(1/beta[b])))
        
    # determine the quantile-dependent limits for the grid search
    lower_min = np.min(quants)                   # we can generally skip the tail below the lower quantile
    max_width = np.diff(quants)                   # maximum width is the fit interval if all data is clean
    min_width = min_clean_fraction*max_width   # minimum width of the fit interval, as fraction of data


    rowval = np.array(np.round(n*(np.arange(lower_min,lower_min+max_dropout_fraction+ (step_sizes[0]*1e-9),step_sizes[0]))))
    colval = np.array(np.arange(0,int(np.round(n*max_width))))
    newX=[]
    for iX in range(len(colval)):
        newX.append(X[np.int_(iX+rowval)])
    
    X1=newX[0]
    newX=newX-np.matlib.repmat(X1,len(colval),1)
    
    opt_val = np.inf;
    
    for m in (np.round(n*np.arange(max_width,min_width,-step_sizes[1]))):
        mcurr = int(m-1)
        nbins = int(np.round(3*np.log2(1+m/2)))
        
        rowval = np.array(nbins/newX[mcurr])
        H=newX[0:int(m)]*np.matlib.repmat(rowval,int(m),1)
        
        HistAll=[]        
        for ih in range(len(rowval)):
            histcurr = np.histogram(H[:,ih], bins=np.arange(0,nbins+1))
            HistAll.append(histcurr[0])
        
        HistAll=np.int_(np.transpose(HistAll))
        HistAll = np.vstack((HistAll,np.zeros(len(rowval), dtype=int)))
        logq = np.log(HistAll + 0.01)
        
        # for each shape value...
        for b in range(len(beta)):
            bounds = zbounds[b]    
            x= bounds[0]+(np.arange(0.5,nbins+0.5)/nbins*np.diff(bounds))
            p = np.exp(-np.abs(x)**beta[b])*rescale[b]
            p=p/np.sum(p)
            
            # calc KL divergences        
            kl = np.sum((np.transpose(np.matlib.repmat(p,logq.shape[1],1)))*((np.transpose(np.matlib.repmat(np.log(p),logq.shape[1],1))) - logq[0:-1,:]),axis=0) + np.log(m)
            
            # update optimal parameters
            min_val = np.min(kl)
            idx =np.argmin(kl)
            
            if (min_val < opt_val):
                opt_val = min_val;
                opt_beta = beta[b];
                opt_bounds = bounds;
                opt_lu = [X1[idx],(X1[idx] + newX[int(m-1),idx])]
                
    
    # recover distribution parameters at optimum
    alpha = (opt_lu[1]-opt_lu[0])/np.diff(opt_bounds);
    mu = opt_lu[0]-opt_bounds[0]*alpha;
    beta = opt_beta;
    
    # calculate the distribution's standard deviation from alpha and beta
    sig = np.sqrt((alpha**2)*scipy.special.gamma(3/beta)/scipy.special.gamma(1/beta))
    
    return mu,sig,alpha,beta



def polystab(a):
    #POLYSTAB Polynomial stabilization.
    #   POLYSTAB(A), where A is a vector of polynomial coefficients,
    #   stabilizes the polynomial with respect to the unit circle;
    #   roots whose magnitudes are greater than one are reflected
    #   inside the unit circle.
    #
    #   # Example:
    #   #   Convert a linear-phase filter into a minimum-phase filter with the 
    #   #   same magnitude response.
    #
    #   h = fir1(25,0.4);               # Window-based FIR filter design
    #   flag_linphase = islinphase(h)   # Determines if filter is linear phase
    #   hmin = polystab(h) * norm(h)/norm(polystab(h)); 
    #   flag_minphase = isminphase(hmin)# Determines if filter is minimum phase
    
    v = np.roots(a);
    i =  np.where(v!=0)
    vs = 0.5*(np.sign(np.abs(v[i])-1)+1);
    v[i] = (1-vs)*v[i] + vs/np.conj(v[i]);
    ind = np.where(a!=0)
    b = a[ind[0][0]]*np.poly(v);
    
#     Return only real coefficients if input was real:
    if not(np.sum(np.imag(a))):
    	b = np.real(b)
    
    return b


def numf(h,a,nb):
    #NUMF	Find numerator B given impulse-response h of B/A and denominator A
    #   NB is the numerator order.  This function is used by YULEWALK.
  

    nh = np.max(h.size); 
    xn=np.concatenate((1,np.zeros((1,nh-1))),axis=None)    
    impr = signal.lfilter(np.array([1.0]),a,xn)
    toeplitz(impr,np.concatenate((1,np.zeros((1,nb))),axis=None))
    
    b = np.linalg.lstsq(toeplitz(impr,np.concatenate((1,np.zeros((1,nb))),axis=None)), h.T,rcond=None)[0].T

    return b





def denf(R,na):
    #DENF	Compute denominator from covariances.
    #   A = DENF(R,NA) computes order NA denominator A from covariances 
    #   R(0)...R(nr) using the Modified Yule-Walker method.  
    #   This function is used by YULEWALK.    
    nr = np.max(np.size(R));
    Rm = toeplitz(R[na:nr-1],R[na:0:-1])
    Rhs = - R[na+1:nr];
    A = np.concatenate((1,np.linalg.lstsq(Rm, Rhs.T,rcond=None)[0].T),axis=None)
    return A





def yulewalk(Order, F, M):
#YULEWALK Recursive filter design using a least-squares method.
#   [B,A] = YULEWALK(N,F,M) finds the N-th order recursive filter
#   coefficients B and A such that the filter:
#   	                      -1             -(n-1) 
#   	   B(z)   b(1) + b(2)z + .... + b(n)z
#   	   ---- = ---------------------------
#   	                      -1             -(n-1)
#   	   A(z)    1   + a(1)z + .... + a(n)z
#
#   matches the magnitude frequency response given by vectors F and M.
#   Vectors F and M specify the frequency and magnitude breakpoints for
#   the filter such that PLOT(F,M) would show a plot of the desired
#   frequency response. The frequencies in F must be between 0.0 and 1.0,
#   with 1.0 corresponding to half the sample rate. They must be in
#   increasing order and start with 0.0 and end with 1.0. 
#
#   # Example:
#   #   Design an 8th-order lowpass filter and overplot the desired  
#   #   frequency response with the actual frequency response.
#
#   f = [0 0.6 0.6 1];      # Frequency breakpoints 
#   m = [1 1 0 0];          # Magnitude breakpoints
#   [b,a] = yulewalk(8,f,m);# Filter design using a least-squares method
#   [h,w] = freqz(b,a,128); # Frequency response of filter
#   plot(f,m,w/pi,abs(h),'--')
#   legend('Ideal','yulewalk Designed')
#   title('Comparison of Frequency Response Magnitudes')
#
#   See also FIR1, BUTTER, CHEBY1, CHEBY2, ELLIP, FREQZ and FILTER.

#   The YULEWALK function performs a least squares fit in the time
#   domain. The denominator coefficients {a(1),...,a(NA)} are computed
#   by the so called "modified Yule Walker" equations, using NR
#   correlation coefficients computed by inverse Fourier transformation
#   of the specified frequency response H.
#   The numerator is computed by a four step procedure. First, a numerator
#   polynomial corresponding to an additive decomposition of the power 
#   frequency response is computed. Next, the complete frequency response
#   corresponding to the numerator and denominator polynomials is
#   evaluated. Then a spectral factorization technique is used to
#   obtain the impulse response of the filter. Finally, the numerator
#   polynomial is obtained by a least squares fit to this impulse
#   response. For a more detailed explanation of the algorithm see 
#   B. Friedlander and B. Porat, "The Modified Yule-Walker Method
#   of ARMA Spectral Estimation," IEEE Transactions on Aerospace
#   Electronic Systems, Vol. AES-20, No. 2, pp. 158-173, March 1984.

    npt = 512
    lap = np.fix(npt/25)       
    mf = F.size
    mm = M.size
    npt = npt + 1  # For [dc 1 2 ... nyquist].
    Ht = np.array(np.zeros((1,npt)))
    nint=mf-1;
    df = np.diff(F)
    
    nb = 0
    Ht[0][0]=M[0]
    for i in range(nint):
        if (df[i] == 0):
            nb = nb - lap/2
            ne = nb + lap
        else:
            ne = int(np.fix(F[i+1]*npt))-1
    
        j=np.arange(nb,ne+1)
        if (ne == nb):
            inc = 0
        else:
            inc = (j-nb)/(ne-nb)
        
        Ht[0][nb:ne+1] = np.array(inc*M[i+1] + (1 - inc)*M[i])
        nb = ne + 1;


    Ht = np.concatenate((Ht,Ht[0][-2:0:-1]), axis=None)
    n = Ht.size
    n2 = np.fix((n+1)/2)
    nb = Order
    nr = 4*Order;
    nt = np.arange(0,nr)
    
    #    compute correlation function of magnitude squared response    
    R = np.real(np.fft.ifft(Ht*Ht))
    R  = R[0:nr]*(0.54+0.46*np.cos(np.pi*nt/(nr-1)))     # pick NR correlations 
    
    #     Form window to be used in extracting the right "wing" of two-sided covariance sequence
    Rwindow = np.concatenate((1/2,np.ones((1,int(n2-1))),np.zeros((1,int(n-n2)))),axis=None) 
    A = tools.polystab(tools.denf(R,Order));            	# compute denominator
    
    Qh = tools.numf(np.concatenate((R[0]/2,R[1:nr]),axis=None),A,Order);	# compute additive decomposition
    
    _,Ss = 2*np.real(scipy.signal.freqz(Qh,A, worN=n,whole=True)) # compute impulse response
    
    hh = np.real(np.fft.ifft(np.exp(np.fft.fft( Rwindow*np.real(np.fft.ifft(np.log(Ss))) ))))
    B  = np.real(tools.numf(hh[0:nr],A,nb));
    
    return B,A
        
        
    
    






def block_geometric_median(X,blocksize):
    # Calculate a blockwise geometric median (faster and less memory-intensive 
    # than the regular geom_median function).
    #
    # This statistic is not robust to artifacts that persist over a duration that
    # is significantly shorter than the blocksize.
    #
    # In:
    #   X : the data (#observations x #variables)
    #   blocksize : the number of successive samples over which a regular mean 
    #               should be taken
    #   tol : tolerance (default: 1.e-5)
    #   y : initial value (default: median(X))
    #   max_iter : max number of iterations (default: 500)
    #
    # Out:
    #   g : geometric median over X
    #
    # Notes:
    #   This function is noticably faster if the length of the data is divisible by the block size.
    #   Uses the GPU if available.
    # 
    
    if (blocksize > 1):
        o,v=X.shape       # #observations & #variables
        r = np.mod(o,blocksize)  # #rest in last block
        b = int((o-r)/blocksize)   # #blocks
        Xreshape = np.zeros((b+1,v))
        if (r > 0):
           Xreshape[0:b,:] = np.reshape( np.sum(np.reshape(X[0:(o-r),:],(blocksize,b*v)),axis=0), (b,v))
           Xreshape[b,:] = np.sum(X[(o-r+1):o,:],axis=0)*(blocksize/r)
        else:
           Xreshape =   np.reshape(np.sum(np.reshape(X,(blocksize,b*v)),axis=0),(b,v))
        X=Xreshape
    tol = 1.e-5
    y = np.median(X,axis=0)
    max_iter = 500
    y = geometric_median(X,tol,y,max_iter)/blocksize
    
    return y

def geometric_median(X,tol,y,max_iter):
# Calculate the geometric median for a set of observations (mean under a Laplacian noise distribution)
# This is using Weiszfeld's algorithm.
#
# In:
#   X : the data, as in mean
#   tol : tolerance (default: 1.e-5)
#   y : initial value (default: median(X))
#   max_iter : max number of iterations (default: 500)
#
# Out:
#   g : geometric median over X
    
      
    for i in range(max_iter):        
        invnorms =1/np.sqrt(np.sum((X - np.matlib.repmat(y,X.shape[0],1))**2,axis=1))
        oldy = y
        y = np.sum(X * np.transpose(np.matlib.repmat(invnorms,X.shape[1],1)),axis=0) / np.sum(invnorms)
        
        if ((np.linalg.norm(y-oldy)/np.linalg.norm(y)) < tol):
            break
    
    return y 
    




def moving_average(N,X,Zi):
    # Run a moving-average filter along the second dimension of the data.
    # X,Zf = moving_average(N,X,Zi)
    #
    # In:
    #   N : filter length in samples
    #   X : data matrix [#Channels x #Samples]
    #   Zi : initial filter conditions (default: [])
    #
    # Out:
    #   X : the filtered data
    #   Zf : final filter conditions
    #
    #                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
    #                           2012-01-10
    #
    #if nargin <= 2 || isempty(Zi)
    #    Zi = zeros(size(X,1),N); end
    #
    # pre-pend initial state & get dimensions
    #Y = [Zi X]; M = size(Y,2);
    # get alternating index vector (for additions & subtractions)
    #I = [1:M-N; 1+N:M];
    # get sign vector (also alternating, and includes the scaling)
    #S = [-ones(1,M-N); ones(1,M-N)]/N;
    # run moving average
    #X = cumsum(bsxfun(@times,Y(:,I(:)),S(:)'),2);
    # read out result
    #X = X(:,2:2:end);
    #
    #if nargout > 1
    #    Zf = [-(X(:,end)*N-Y(:,end-N+1)) Y(:,end-N+2:end)]; end
    #
    #
    
    [C,S] = X.shape

    if (Zi==None):
        Zi = np.zeros((C,N))
    
    # pre-pend initial state & get dimensions    
    Y = np.concatenate((Zi, X), axis=1)
    [CC,M]=Y.shape
    
    # get alternating index vector (for additions & subtractions)    
    I = np.vstack((np.arange(0,M-N),np.arange(N,M)))
    
    # get sign vector (also alternating, and includes the scaling)
    S = np.vstack((- np.ones((1,M-N)),np.ones((1,M-N))))/N
    
    # run moving average
    YS = np.zeros((C,S.shape[1]*2))
    for i in range(C):
        YS[i,:] = Y[i,I.flatten(order='F')]*S.flatten(order='F')
    
    X = np.cumsum(YS,axis=1)
    # read out result
    X = X[:,1::2]
    
    Zf =np.transpose(np.vstack(((-((X[:,-1]*N)-Y[:,-N])),np.transpose(Y[:,-N+1:]))))
    
    return X,Zf
    
    



def ReadHDF5(Template_H5Filename):
    f = h5py.File(Template_H5Filename, 'r')
    
    print('## Lecture du fichier {}'.format(Template_H5Filename))
    
    dictio = {}
    for element in f:
        groupe = f[element]
                    
        for element in groupe:
            dictio[element] = groupe[element]
            


    return dictio      