import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.0*sigma**2))

class ASIC_Ana:
    def __init__(self, ftrigger, fout):
        self.dict_chan = {}
        points = 500/33.
        for i in range(33):
            self.dict_chan[str(i)] = []
        for i in range(500):
            chan = int(i/points)
            if i - chan*points  > 1 and i - chan*points < 33:
                self.dict_chan[str(chan)].append(i)
        self.trigger = pd.read_csv(ftrigger, skiprows=4)
        self.out = pd.read_csv(fout, skiprows=4)
        trig = np.array(self.trigger['Ampl'])
        out_amp = np.array(self.out['Ampl'])
        self.amp_chan = {}
        for i in range(33):
            self.amp_chan[str(i)] = np.array([])

        for i in range(len(out_amp) - 500):
            if trig[i] < 0.3 and trig[i+1] > 0.3:
                for j in range(33):
                    chans = self.dict_chan[str(j)]
                    val = 0
                    batch = np.zeros(len(chans))
                    for n, chan in enumerate(chans):
                        if i+chan < len(trig):
                            batch[n] = out_amp[i+chan]
                    if i < 1100:
                        self.amp_chan[str(j)] = np.array([batch])
                    else:
                        self.amp_chan[str(j)] = np.vstack([self.amp_chan[str(j)], batch])

    def LowFilter(self, cutoff):
        from scipy import signal
        b, a = signal.butter(5, cutoff, fs=2000000, btype='low', analog=False)
        self.filt_chan = {}
        for key in self.amp_chan:
            self.filt_chan[key] = self.amp_chan[key]
            for n in range(self.filt_chan[key].shape[1]):
                self.filt_chan[key][:, n] = signal.filtfilt(b, a, self.filt_chan[key][:, n])

    def plot_tick(self, savename):
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(self.trigger['Time'][:1000], self.trigger['Ampl'][:1000], c='red')
        ax[1].plot(self.out['Time'][:1000], self.out['Ampl'][:1000], c='blue')
        plt.savefig(savename)
        plt.close()

    def plot_channels(self, chs, savename, Filtered: bool =False):
        fig, ax = plt.subplots(len(chs), 2)
        if len(chs) == 1:
            sch = str(chs[0])
            if Filtered:
                ax[0].plot(np.arange(self.filt_chan[sch].shape[0]), self.filt_chan[sch][:,0])
            else:
                ax[0].plot(np.arange(self.amp_chan[sch].shape[0]), self.amp_chan[sch][:,0])
            ax[0].set_xlabel('Time ($\mu s$)')
            ax[0].set_ylabel('Ch %d magnitude (V)' % chs[0])
            lwf = self.amp_chan[sch].shape[0]
            if Filtered:
                ax[1].plot(np.arange(int(lwf/2)-150, int(lwf/2)+150), self.filt_chan[sch][int(lwf/2)-150:int(lwf/2)+150, 0],color='red', label='Filtered' )
                ax[1].plot(np.arange(int(lwf/2)-150, int(lwf/2)+150), self.amp_chan[sch][int(lwf/2)-150:int(lwf/2)+150, 0], color='blue', label='Raw' )
            else:
                ax[1].plot(np.arange(int(lwf/2)-150, int(lwf/2)+150), self.amp_chan[sch][int(lwf/2)-150:int(lwf/2)+150, 0] )
            ax[1].set_xlabel('Time ($\mu s$)')

        else:
            for i in range(len(chs)):
                sch = str(chs[i])
                if Filtered:
                    ax[i,0].plot(np.arange(self.filt_chan[sch].shape[0]), self.filt_chan[sch][:,0])
                else:
                    ax[i,0].plot(np.arange(self.amp_chan[sch].shape[0]), self.amp_chan[sch][:,0])
                ax[i,0].set_xlabel('Time ($\mu s$)')
                ax[i,0].set_ylabel('Ch %d magnitude (V)' % chs[i])
                lwf = self.amp_chan[sch].shape[0]
                if Filtered:
                    ax[i,1].plot(np.arange(int(lwf/2)-150, int(lwf/2)+150), self.filt_chan[sch][int(lwf/2)-150:int(lwf/2)+150, 0],color='red', label='Filtered' )
                    ax[i,1].plot(np.arange(int(lwf/2)-150, int(lwf/2)+150), self.amp_chan[sch][int(lwf/2)-150:int(lwf/2)+150, 0], color='blue', label='Raw' )
                else:
                    ax[i,1].plot(np.arange(int(lwf/2)-150, int(lwf/2)+150), self.amp_chan[sch][int(lwf/2)-150:int(lwf/2)+150, 0] )
                ax[i,1].set_xlabel('Time ($\mu s$)')

        plt.savefig(savename)
        plt.close()

    def CalcNoise(self, savename, dp:int=5, Filtered:bool=False):
        fig1, ax1 = plt.subplots(6,6, figsize=(15,15))
        for i in range(6):
            for j in range(6):
                if j*6+i > 31:
                    continue
                lwf = self.amp_chan[str(i)].shape[0]
                ch_mean = np.mean(self.amp_chan[str(j*6+i)][:int(lwf/2)-100, dp])
                ch_std = np.std(self.amp_chan[str(j*6+i)][:int(lwf/2)-100, dp])
                v_bins = np.linspace(ch_mean - 3*ch_std, ch_mean + 3*ch_std, 30)
                n, bins, patches = ax1[i, j].hist(self.amp_chan[str(j*6+i)][:int(lwf/2)-100, dp], bins = v_bins, histtype='step', color='blue')
                n_avg, bins_avg, patches = ax1[i, j].hist(np.mean(self.amp_chan[str(j*6+i)][:int(lwf/2)-100, 3:-3], axis=1), bins = v_bins, histtype='step', color='red')
                if Filtered:
                    n, bins, patches = ax1[i, j].hist(self.filt_chan[str(j*6+i)][:int(lwf/2)-100, dp], bins = v_bins, histtype='step', color='blue')
                    n_avg, bins_avg, patches = ax1[i, j].hist(np.mean(self.filt_chan[str(j*6+i)][:int(lwf/2)-100, 3:-3], axis=1), bins = v_bins, histtype='step', color='red')
                bin_centers = (bins_avg[:-1] + bins_avg[1:])/2
                p0=[4000, np.mean(bins), ch_std]
                coeff, var_matric = curve_fit(gauss, bin_centers, n_avg, p0=p0)
                ax1[i, j].plot(bin_centers, gauss(bin_centers, *coeff), color='red')
                ax1[i, j].text(bins[1], 100, r'$\sigma$=%.4f V' % coeff[2])
        plt.tight_layout()
        fig1.savefig(savename)
        plt.close()

    def NoiseFreq(self, ch):
        from scipy.fft import fft, fftfreq
        T = 1/2000000.
        N = self.amp_chan[str(ch)].shape[0]
        batchsize = self.amp_chan[str(ch)].shape[1]
        batch_noise = []
        fig, ax = plt.subplots()
        for n in range(batchsize):
            yf = fft(self.amp_chan[str(ch)][:, n])
            xf = fftfreq(N, T)[:N//2]
            batch_noise.append(yf)
            ax.plot(xf, 2.0/N * np.abs(yf[0:N//2]), color='gray')
        mean_noise = np.mean(np.absolute(np.array(batch_noise)), axis=0)
        avg_wf = np.mean(self.amp_chan[str(ch)][:,3:-3], axis=1)
        avg_noise = fft(avg_wf)
        ax.plot(xf, 2.0/N * mean_noise[0:N//2], color='black', label='Mean of noise')
        ax.plot(xf, 2.0/N * avg_noise[0:N//2], color='red', label='Noise of mean wf')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Arbitrary Unit.')
        ax.set_ylim(0.1**7, 0.1**3)
        ax.legend()
        fig.savefig('plots/noise_fft_%d.pdf' % ch)
        plt.close()

if __name__ == '__main__':
    #ftrigger = '/scratchfs/exo/zepengli94/nexo/pcb_20210418/C3c4trigger_c3s1_c2sh31_c1out_DC50_1ms_CalON10fC_CLKON66MHz_10ch00000.csv'
    #fout = '/scratchfs/exo/zepengli94/nexo/pcb_20210418/C1c4trigger_c3s1_c2sh31_c1out_DC50_1ms_CalON10fC_CLKON66MHz_10ch00000.csv'
    ftrigger = '/scratchfs/exo/zepengli94/nexo/pcb_v20_20210420/C3c4trigger_c3s1_c2sh31_c1out_DC50_1ms_CalON10fC_CLKON66MHz_25ch00000.csv'
    fout = '/scratchfs/exo/zepengli94/nexo/pcb_v20_20210420/C1c4trigger_c3s1_c2sh31_c1out_DC50_1ms_CalON10fC_CLKON66MHz_25ch00000.csv'
    ana = ASIC_Ana(ftrigger, fout)
    ana.plot_tick('plots/OneTick.pdf' )
    for n in range(5, 12):
        ana.CalcNoise('plots/Noise_dp_%d.pdf' %  n, n)
    ana.LowFilter(200000)
    for i in range(14):
        ana.plot_channels([i], 'plots/chans_%d.pdf' % i)
        ana.plot_channels([i], 'plots/filt_chans_%d.pdf' % i, True)
        ana.NoiseFreq(i)
