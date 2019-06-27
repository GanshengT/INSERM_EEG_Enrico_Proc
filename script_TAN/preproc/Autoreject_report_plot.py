import mne
import importlib
import numpy as np
import os
from autoreject import AutoReject
from autoreject import compute_thresholds
from autoreject import get_rejection_threshold 
import matplotlib.pyplot as plt  # noqa
import matplotlib.patches as patches  # noqa
from autoreject import set_matplotlib_defaults  # noqa


def Autoreject_report_plot(ar,epochs,epochs_autorejed,reject_log,plot_err_matrix = False,plot_epochs_rej_intpl = True,
                           plot_epochs_rejed =False,plot_drop_rate = True,plot_epochs_bfaft_compare = False,
                           plot_chan_hist_thresh = False):
                               
    mne.set_log_level('WARNING')

    ########################### second Autorejection plot ############################################################################
    if plot_err_matrix == True:
        set_matplotlib_defaults(plt, style='seaborn-white')
        loss = ar.loss_['eeg'].mean(axis=-1)  # losses are stored by channel type.

        plt.matshow(loss.T * 1e6, cmap=plt.get_cmap('viridis'))
        plt.xticks(range(len(ar.consensus)), ['%.1f' % c for c in ar.consensus])
        plt.yticks(range(len(ar.n_interpolate)), ar.n_interpolate)

        # Draw rectangle at location of best parameters
        ax = plt.gca()
        idx, jdx = np.unravel_index(loss.argmin(), loss.shape)
        rect = patches.Rectangle((idx - 0.5, jdx - 0.5), 1, 1, linewidth=2,
                                 edgecolor='r', facecolor='none')
        ax.add_patch(rect)
        ax.xaxis.set_ticks_position('bottom')
        plt.xlabel(r'Consensus percentage $\kappa$')
        plt.ylabel(r'Max sensors interpolated $\rho$')
        plt.title('Mean cross validation error (x 1e6)')
        plt.colorbar()
        plt.show()

    if plot_epochs_rej_intpl == True:
        reject_log.plot_epochs(epochs)
    if plot_epochs_rejed == True:
        epochs_autorejed.plot()
    if plot_drop_rate ==True:
        epochs_autorejed.plot_drop_log()
    if plot_epochs_bfaft_compare == True:
        evoked_clean = epochs_autorejed.average()
        evoked = epochs.average()
        set_matplotlib_defaults(plt)

        fig, axes = plt.subplots(2, 1, figsize=(6, 6))

        for ax in axes:
            ax.tick_params(axis='x', which='both', bottom='off', top='off')
            ax.tick_params(axis='y', which='both', left='off', right='off')

        evoked.plot(exclude=[], axes=axes[0], show=False)
        axes[0].set_title('Before autoreject')
        evoked_clean.plot(exclude=[], axes=axes[1])
        axes[1].set_title('After autoreject')
        plt.tight_layout()

    if plot_chan_hist_thresh == True:
        threshes = ar.threshes_
        set_matplotlib_defaults(plt)
        unit = r'uV'
        scaling = 1e6

        plt.figure(figsize=(6, 5))
        plt.tick_params(axis='x', which='both', bottom='off', top='off')
        plt.tick_params(axis='y', which='both', left='off', right='off')

        plt.hist(scaling * np.array(list(threshes.values())), 30,
                 color='g', alpha=0.4)
        plt.xlabel('Threshold (%s)' % unit)
        plt.ylabel('Number of sensors')
        plt.tight_layout()
        plt.show() 
