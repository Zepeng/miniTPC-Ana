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
        points = 1000/33.
        for i in range(33):
            self.dict_chan[str(i)] = []
        for i in range(1000):
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

        for i in range(len(out_amp) - 1000):
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
        ax[0].plot(self.trigger['Time'][:3000], self.trigger['Ampl'][:3000], c='red')
        ax[1].plot(self.out['Time'][:3000], self.out['Ampl'][:3000], c='blue')
        plt.savefig(savename)
        plt.close()

    def plot_channels(self, chs, savename, Filtered: bool =False):
        fig, ax = plt.subplots(len(chs), 2)
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

    def CalcNoise(self, savename):
        fig1, ax1 = plt.subplots(6,6, figsize=(15,15))
        for i in range(6):
            for j in range(6):
                if j*6+i > 31:
                    continue
                lwf = self.amp_chan[str(i)].shape[0]
                n, bins, patches = ax1[i, j].hist(np.array(self.amp_chan[str(j*6+i)][:int(lwf/2)-100, 0]), bins=20, histtype='step', color='blue')
                bin_centers = (bins[:-1] + bins[1:])/2
                p0=[200, np.mean(bins), 0.001]
                coeff, var_matric = curve_fit(gauss, bin_centers, n, p0=p0)
                ax1[i, j].plot(bin_centers, gauss(bin_centers, *coeff), color='red')
                ax1[i, j].text(bins[1], 100, r'$\sigma$=%.4f V' % coeff[2])
        plt.tight_layout()
        fig1.savefig(savename)
        plt.close()

    def CalcFiltNoise(self, savename):
        fig1, ax1 = plt.subplots(6,6, figsize=(15,15))
        for i in range(6):
            for j in range(6):
                if j*6+i > 31:
                    continue
                lwf = self.filt_chan[str(i)].shape[0]
                n, bins, patches = ax1[i, j].hist(np.array(self.filt_chan[str(j*6+i)][:int(lwf/2)-100, 0]), bins=20, histtype='step', color='blue')
                bin_centers = (bins[:-1] + bins[1:])/2
                p0=[200, np.mean(bins), 0.001]
                coeff, var_matric = curve_fit(gauss, bin_centers, n, p0=p0)
                ax1[i, j].plot(bin_centers, gauss(bin_centers, *coeff), color='red')
                ax1[i, j].text(bins[1], 100, r'$\sigma$=%.4f V' % coeff[2])
        plt.tight_layout()
        fig1.savefig(savename)
        plt.close()
    def NoiseFreq(self, ch):
        from scipy.fft import fft, fftfreq
        T = 1/1000000.
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
        ax.plot(xf, 2.0/N * mean_noise[0:N//2], color='black')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Arbitrary Unit.')
        ax.set_ylim(0.1**7, 0.1**3)
        fig.savefig('plots/noise_fft_%d.pdf' % ch)
        plt.close()

if __name__ == '__main__':
    ftrigger = '/workfs2/exo/zepengli94/nexo/pcb_20210416/C3c4trigger_c3s1_c2sh31_c1out_DC50_1ms_CalON10fC_CLKON33MHz_10ch00001.csv'
    fout = '/workfs2/exo/zepengli94/nexo/pcb_20210416/C1c4trigger_c3s1_c2sh31_c1out_DC50_1ms_CalON10fC_CLKON33MHz_10ch00001.csv'
    ana = ASIC_Ana(ftrigger, fout)
    ana.LowFilter(200000)
    for i in range(14):
        ana.plot_tick('plots/OneTick_%d.pdf' % i)
        ana.plot_channels([i, 21], 'plots/chans_%d.pdf' % i)
        ana.plot_channels([i, 21], 'plots/filt_chans_%d.pdf' % i, True)
        ana.CalcNoise('plots/Noise_%d.pdf' % i)
        ana.CalcFiltNoise('plots/FiltNoise_%d.pdf' % i)
        ana.NoiseFreq(i)
