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
    def __init__(self, datafile, ch, starters, nrows):
        """
        look for a data file assuming the files are stored with the same
        pattern
        """
        wfs = np.zeros((len(starters), nrows))
        filtered_wfs = np.zeros((len(starters), nrows))
        from scipy import signal
        b, a = signal.butter(5, 0.05, btype='low', analog=False)
        for starter in starters:
            ds = pd.read_csv(datafile, skiprows=starter, header=None, delimiter='\t', nrows=nrows)
            print(np.max(ds[ch][:1000]))
            wfs[starters.index(starter), :] = ds[ch]
            filtered_wfs[starters.index(starter), :] = signal.filtfilt(b, a, np.array(ds[ch]))
        self.wfs = wfs
        self.filtered_wfs = filtered_wfs
        self.ch = ch


    def PlotEvent(self):
        """
        correct the baseline by the mean value
        and plot the baseline
        """
        fig, ax = plt.subplots(3, 1)
        x = np.arange(self.wfs.shape[1])*0.4
        print(self.wfs.shape)
        for i in range(self.wfs.shape[0]):
            if i == 0:
                ax[0].plot(x, self.wfs[i,:], color='blue', label='Raw')
                ax[0].plot(x, self.filtered_wfs[i,:], color='red', label='Filtered')
            ax[0].plot(x, self.wfs[i,:], color='blue')
            ax[0].plot(x, self.filtered_wfs[i,:], color='red')
            ax[0].set_xlabel('Time (ns)')
            ax[0].set_ylabel('Magnitude (V)')
            ax[0].set_title('Channel {} 53 V Bias'.format(self.ch))
            ax[0].set_ylim(np.min(self.wfs) - 0.01, np.max(self.wfs) + 0.02)
        ax[0].legend(ncol=2)

        baseline = self.wfs[:,1000]
        filtered_baseline = self.filtered_wfs[:, 1000]
        n, bins, patch = ax[1].hist(baseline.flatten(), bins=np.linspace(np.min(baseline), np.max(baseline), 16), histtype='step', color='blue', label='Raw')
        ax[1].hist(filtered_baseline.flatten(), bins=np.linspace(np.min(baseline), np.max(baseline), 16), histtype='step', color='red', label='Filtered')
        ax[1].text(np.min(baseline), np.mean(n), r'$\sigma$=%.4f V' % np.std(filtered_baseline.flatten()))
        ax[1].set_xlabel('Baseline (V)')
        ax[1].set_title('Channel {} Baseline'.format(self.ch))
        ax[1].legend()
        plt.tight_layout()
        fig.savefig('micro_ch{}.png'.format(self.ch))

    def CalcQE(self):
        QEs = np.zeros(self.filtered_wfs.shape[0])
        for i in range(self.filtered_wfs.shape[0]):
            QEs[i] = np.sum(self.filtered_wfs[i,1800:2500]) - np.mean(self.filtered_wfs[i,:1000])*700
        self.QEs = QEs
        fig, ax = plt.subplots()
        ax.hist(QEs, histtype='step', color='blue', bins=np.linspace(np.min(QEs), np.max(QEs), 100))
        ax.set_xlabel('Integration of waveform')
        ax.set_ylabel('Events')
        fig.savefig('QE.png')



if __name__ == '__main__':
    #ana = SiPMAna('/workfs2/juno/weiwl/new_board/sipm3_30.txt', '/workfs2/juno/weiwl/new_board/sipm3.txt', nrows=1500000 )
    datafile = '/scratchfs/exo/zepengli94/nexo/micro_ch2.txt'
    lines =  []
    with open(datafile) as dfile:
        for num, line in enumerate(dfile, 1):
            if 'X_Value' in line:
                lines.append(num)
    ana = SiPMAna('/scratchfs/exo/zepengli94/nexo/micro_ch2.txt', ch=2, starters = lines, nrows=9999)
    ana.PlotEvent()
    ana.CalcQE()
