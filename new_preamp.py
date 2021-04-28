import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd
from scipy.optimize import curve_fit

def func(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y

class SiPMAna:
    """
    SiPMAna is a analyzer to do different analysis of the
    SiPM test data take from Tektronix scope
    """
    def __init__(self, file_lowV, file_highV, nrows):
        """
        look for a data file assuming the files are stored with the same
        pattern
        """
        ds_low = pd.read_csv(file_lowV, skiprows=22, header=None, delimiter='\t', nrows=nrows)
        ds_high = pd.read_csv(file_highV, skiprows=22, header=None, delimiter='\t', nrows=nrows)
        self.ds_low = ds_low
        self.ds_high = ds_high


    def PlotEvent(self):
        """
        correct the baseline by the mean value
        and plot the baseline
        """
        from scipy import signal
        b, a = signal.butter(5, 0.05, btype='low', analog=False)
        fig, ax = plt.subplots(3, 2)
        x = np.arange(len(self.ds_low[1]))*0.4
        x_short = np.arange(10000)*0.4
        filtered_wf1 = signal.filtfilt(b, a, np.array(self.ds_low[1]))
        filtered_wf2 = signal.filtfilt(b, a, np.array(self.ds_low[2]))
        ax[0, 0].plot(x, np.array(self.ds_low[1]), color='blue', label='Raw')
        ax[0, 0].plot(x, filtered_wf1, color='red', label='Filtered')
        ax[0, 0].set_xlabel('Time (ns)')
        ax[0, 0].set_ylabel('Magnitude (V)')
        ax[0, 0].set_title('Channel 1 30 V Bias')
        ax[0, 0].set_ylim(np.min(self.ds_low[1])-0.01, np.max(self.ds_low[1]) + 0.02)
        ax[0, 0].legend(ncol=2)
        ax[0, 1].plot(x, np.array(self.ds_low[2]), color='blue')
        ax[0, 1].plot(x, filtered_wf2, color='red')
        ax[0, 1].set_xlabel('Time (ns)')
        ax[0, 1].set_ylabel('Magnitude (V)')
        ax[0, 1].set_title('Channel 2 30 V Bias')
        ax[0, 1].set_ylim(np.min(self.ds_low[2])-0.01, np.max(self.ds_low[2]) + 0.02)
        ax[1, 0].plot(x_short, np.array(self.ds_high[1])[:10000], color='blue', label='Raw Waveform')
        ax[1, 0].plot(x_short, signal.filtfilt(b, a, np.array(self.ds_high[1])[:10000]), color='red', label='Filtered Waveform')
        ax[1, 0].set_xlabel('Time (ns)')
        ax[1, 0].set_ylabel('Magnitude (V)')
        ax[1, 0].set_title('Channel 1 54 V Bias')
        #ax[1, 0].legend()
        ax[1, 1].plot(x_short, np.array(self.ds_high[2])[:10000], color='blue', label='Raw Waveform')
        ax[1, 1].plot(x_short, signal.filtfilt(b, a, np.array(self.ds_high[2])[:10000]), color='red', label='Filtered Waveform')
        ax[1, 1].set_xlabel('Time (ns)')
        ax[1, 1].set_ylabel('Magnitude (V)')
        ax[1, 1].set_title('Channel 1 54 V Bias')
        #ax[1, 1].legend()
        n, bins, patch = ax[2, 0].hist(np.array(self.ds_low[1]), bins=np.linspace(np.min(self.ds_low[1]), np.max(self.ds_low[1]), 16), histtype='step', color='blue', label='Raw')
        ax[2, 0].hist(filtered_wf1, bins=np.linspace(np.min(self.ds_low[1]), np.max(self.ds_low[1]), 16), histtype='step', color='red', label='Filtered')
        ax[2, 0].text(np.min(self.ds_low[1]), np.mean(n), r'$\sigma$=%.4f V' % np.std(filtered_wf1))
        ax[2, 0].set_xlabel('Baseline (V)')
        ax[2, 0].set_title('Channel 1 Baseline')
        #ax[2, 0].legend()
        n, bins, pathc = ax[2, 1].hist(np.array(self.ds_low[2]), bins=np.linspace(np.min(self.ds_low[2]), np.max(self.ds_low[2]), 16), histtype='step', color='blue', label='Raw')
        ax[2, 1].hist(filtered_wf2, bins=np.linspace(np.min(self.ds_low[2]), np.max(self.ds_low[2]), 16), histtype='step', color='red', label='Filtered')
        ax[2, 1].text(np.min(self.ds_low[2]), np.mean(n), r'$\sigma$=%.4f V' % np.std(filtered_wf2))
        ax[2, 1].set_xlabel('Baseline (V)')
        ax[2, 1].set_title('Channel 2 Baseline')
        #ax[2, 1].legend()
        plt.tight_layout()
        fig.savefig('wfs.png')

    def LowFilter(self):
        from scipy import signal
        fwfs = []
        b, a = signal.butter(5, 0.05, btype='low', analog=False)
        for wf in self.waveforms:
            fwf= signal.filtfilt(b, a, wf)
            fwfs.append(fwf)
        self.filteredwfs = np.array(fwfs)


if __name__ == '__main__':
    #ana = SiPMAna('/workfs2/juno/weiwl/new_board/sipm3_30.txt', '/workfs2/juno/weiwl/new_board/sipm3.txt', nrows=1500000 )
    ana = SiPMAna('/workfs2/juno/weiwl/new_board/9sipm_30v_3.5v.txt', '/workfs2/juno/weiwl/new_board/9sipm_54v_3.5v.txt', nrows=100000 )
    ana.PlotEvent()
