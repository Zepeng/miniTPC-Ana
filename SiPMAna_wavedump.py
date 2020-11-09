import numpy as np
import matplotlib.pyplot as plt

def func(x, a, b, c, d, e):
    return a*np.exp(-1*(x-b)**2/(2*c**2)) + d*np.exp(-1*(x-e)**2/(2*c**2))

class SiPMAna:
    waveforms = 0
    filteredwfs = 0
    charges = 0
    filepath = ''
    qdcs = 0
    def __init__(self):
        self.waveforms = 0
        self.filteredwfs = 0
        self.charges = 0
        self.filepath = 0
    def LoadData(self, fpath, ch):
        self.filepath = fpath + '/wave_%d.txt' % ch
        qdc = []
        import os
        if os.path.exists('qdc_%d.npy' % ch):
            return
        with open(self.filepath) as f:
            i = 0
            wfs = []
            count = 0
            wf = []
            for line in f:
                item = line.split()[0]
                wf.append(float(item))
                if i%1024 ==0 and i>1:
                    wfs.append(wf)
                    if ch == 0 or ch ==6:
                        qdc.append(self.CalcQ(self.LowFilter(wf), 600))
                    else:
                        qdc.append(self.CalcQ(self.LowFilter(wf), 450))
                    wf = []
                    count += 1
                i+=1
        self.charges = np.array(qdc)
        np.save('qdc_%d.npy' % ch,  np.array(qdc))

    def SumEvent(self):
        channels = []
        for i in range(12):
            channels.append(np.load('qdc_%d.npy' % i))
        npchannels = np.array(channels)
        EventEs = []
        bi1MeV = []
        for i in range(channels[0].shape[0]):
            eventq = -1*np.sum(npchannels[:, i])
            EventEs.append(eventq)
            if eventq > 1600000 and eventq < 2000000:
                bi1MeV.append(eventq)
        fig, ax = plt.subplots()
        from scipy.stats import norm
        from scipy.optimize import curve_fit
        bins=np.linspace(1600000,2000000,100)
        data_entries, bins_1 = np.histogram(bi1MeV, bins=np.linspace(1600000,2000000,100))
        binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
        (mu, sigma) = norm.fit(bi1MeV)
        popt, pcov = curve_fit(func, xdata=binscenters, ydata=data_entries, p0=[1400,175000,0.04*1730000, 1, 1850000])
        ax.hist(bi1MeV, bins)
        #print(mu, sigma)
        #print(binscenters, data_entries)
        print(popt)
        ax.plot(binscenters, func(binscenters, *popt))
        fig.savefig('energy_n40_27V.png')
    def PlotWFs(self, entries = [0]):
        for entry in entries:
            fig, ax = plt.subplots()
            if entry > np.size(self.waveforms, 0):
                print('Entry not found in saved waveforms')
                return
            ax.plot(np.arange(len(self.waveforms[entry])), self.waveforms[entry])
            fig.savefig('wfs_%d.png' % entry)

    def LowFilter(self, wf):
        from scipy import signal
        b, a = signal.butter(5, 0.1, btype='low', analog=False)
        fwf= signal.filtfilt(b, a, wf)
        return fwf

    def CalcQ(self, wf, window):
        return np.sum(wf[400:400 + window]) - window*np.mean(wf[100:300])

if __name__ == '__main__':
    ana = SiPMAna()
    for ch in range(12):
        SiPMAna.LoadData(ana, './FADC_Bi207_N40_27V', ch)
    ana.SumEvent()
