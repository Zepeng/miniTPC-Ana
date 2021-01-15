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
            if abs(chan*points - i) > 5:
                self.dict_chan[str(chan)].append(i)
        self.trigger = pd.read_csv(ftrigger, skiprows=4)
        self.out = pd.read_csv(fout, skiprows=4)
        trig = np.array(self.trigger['Ampl'])
        out_amp = np.array(self.out['Ampl'])
        print(trig.shape, out_amp.shape)
        self.amp_chan = {}
        for i in range(33):
            self.amp_chan[str(i)] = []

        for i in range(len(out_amp) - 1000):
            if trig[i] < 0.3 and trig[i+1] > 0.3:
                for j in range(33):
                    chans = self.dict_chan[str(j)]
                    val = 0
                    for chan in chans:
                        if i+chan < len(trig):
                            val += out_amp[i+chan]
                    self.amp_chan[str(j)].append(val/len(chans))

    def plot_tick(self, savename):
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(self.trigger['Time'][:3000], self.trigger['Ampl'][:3000], c='red')
        ax[1].plot(self.out['Time'][:3000], self.out['Ampl'][:3000], c='blue')
        plt.savefig(savename)

    def plot_channels(self, chs, savename):
        fig, ax = plt.subplots(len(chs), 2)
        for i in range(len(chs)):
            sch = str(chs[i])
            ax[i,0].plot(np.arange(len(self.amp_chan[sch])), self.amp_chan[sch])
            ax[i,0].set_xlabel('Time ($\mu s$)')
            ax[i,0].set_ylabel('Ch %d magnitude (V)' % chs[i])
            lwf = len(self.amp_chan[sch])
            ax[i,1].plot(np.arange(int(lwf/2)-150, int(lwf/2)+150), self.amp_chan[sch][int(lwf/2)-150:int(lwf/2)+150] )
            ax[i,1].set_xlabel('Time ($\mu s$)')

        plt.savefig(savename)

    def CalcNoise(self, savename):
        fig1, ax1 = plt.subplots(6,6, figsize=(15,15))
        for i in range(6):
            for j in range(6):
                if j*6+i > 31:
                    continue
                lwf = len(self.amp_chan[str(i)])
                n, bins, patches = ax1[i, j].hist(np.array(self.amp_chan[str(j*6+i)][:int(lwf/2)-100]), bins=20, histtype='step', color='blue')
                bin_centers = (bins[:-1] + bins[1:])/2
                p0=[200, np.mean(bins), 0.001]
                coeff, var_matric = curve_fit(gauss, bin_centers, n, p0=p0)
                ax1[i, j].plot(bin_centers, gauss(bin_centers, *coeff), color='red')
                ax1[i, j].text(bins[1], 100, r'$\sigma$=%.4f V' % coeff[2])
        plt.tight_layout()
        fig1.savefig(savename)

if __name__ == '__main__':
    ftrigger = '/scratchfs/exo/zepengli94/pcbtile20210115/C3c4trigger_c3s1_c1out_c2s1_AC1M_50HzInCalOFFRealOFF33MHz_-5usdiv100mVdiv00001.csv'
    fout = '/scratchfs/exo/zepengli94/pcbtile20210115/C1c4trigger_c3s1_c1out_c2s1_AC1M_50HzInCalOFFRealOFF33MHz_-5usdiv100mVdiv00005.csv'
    ana = ASIC_Ana(ftrigger, fout)
    ana.plot_tick('OneTick_N0.pdf')
    ana.plot_channels([20, 21, 22], 'chans_N0.pdf')
    ana.CalcNoise('Noise_N0.pdf')
