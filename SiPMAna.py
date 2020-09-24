import numpy as np
import matplotlib.pyplot as plt
import scipy
import os, sys

class SiPMAna:
    """
    SiPMAna is a analyzer to do different analysis of the
    SiPM test data take from TD5751
    """
    waveforms = 0
    filteredwfs = 0
    charges = 0
    wfpath = ''
    qdcpath = ''
    wflen = 0
    bSaveFig = False
    def __init__(self, filedir, wflen, ch=0, savefig=False):
        self.waveforms = 0 # waveforms
        self.filteredwfs = 0 # waveform after software filter
        self.charges = 0 # a software QDC
        self.qdc = 0 # QDC saved from TD
        self.wflen = wflen #length of waveform used in TD5751
        self.bSaveFig = savefig # set to True to save figures
        """
        look for a data file assuming the files are stored with the same
        pattern
        """
        filename = filedir + '/ch' + str(ch) + '.txt'
        if os.path.exists(filename):
            self.wfpath = filename
        else:
            print('Waveform file NOT found!')
            sys.exit()
        qdcname = filedir + '/testTitle_QDC.txt'
        if os.path.exists(qdcname):
            self.qdcpath = qdcname
        else:
            print('QDC file NOT found!')
            sys.exit()

    def LoadWFs(self):
        with open(self.wfpath) as f:
            i = 0
            wfs = []
            for line in f:
                wf = []
                for item in line.split():
                    wf.append(float(item))
                wfs.append(wf)
                i+=1
                if i % 1000 ==0:
                    print('Loaded %d of waveforms from file.' % i)
            self.waveforms = np.array(wfs)
            print('Loaded %d waveforms from %s' % (i, self.wfpath))

ana = SiPMAna(filedir='/junofs/users/weiwl/sipmtest/105-46/', ch=0, wflen=5005)
SiPMAna.LoadWFs(ana)
"""
    def PlotWFs(self, entries = [0]):
        for entry in entries:
            fig, ax = plt.subplots()
            if entry > np.size(self.waveforms, 0):
                print('Entry not found in saved waveforms')
                return
            ax.plot(np.arange(len(self.waveforms[entry])), self.waveforms[entry])
            fig.savefig('wfs_%d.png' % entry)
    def PlotBaseline(self):
        baseline = []
        for wf in self.filteredwfs:
            for i in range(100, 250):
                baseline.append(wf[i])
        fig, ax = plt.subplots()
        ax.hist(np.array(baseline), bins = np.arange(945, 960, 0.2))
        from scipy.stats import norm
        (mu, sigma) = norm.fit(baseline)
        print(mu, sigma)
        fig.savefig('baseline.png')

    def LowFilter(self):
        from scipy import signal
        fwfs = []
        b, a = signal.butter(5, 0.05, btype='low', analog=False)
        for wf in self.waveforms:
            fwf= signal.filtfilt(b, a, wf)
            fwfs.append(fwf)
        self.filteredwfs = np.array(fwfs)

    def PlotFilterWFs(self, entries = [0]):
        for entry in entries:
            fig, ax = plt.subplots()
            if entry > np.size(self.filteredwfs, 0):
                print('Entry not found in saved waveforms')
                return
            ax.plot(np.arange(len(self.filteredwfs[entry])), self.filteredwfs[entry])
            #ax.set_xlim(700, 1600)
            fig.savefig('filteredwfs_%d.png' % entry)

    def CalcQ(self):
        Qs = []
        for wf in self.filteredwfs:
            baseline = np.mean(wf[100:300])
            q = np.sum(np.mean(wf[350:700]) - baseline)
            Qs.append(q)
        self.charges = np.array(Qs)
        fig, ax = plt.subplots()
        ax.hist(Qs, bins=np.arange(-4, 0, 0.05))
        #ax.set_yscale('log')
        fig.savefig('qdc.png')


if __name__ == '__main__':
    ana = SiPMAna()
    SiPMAna.LoadData(ana, '/junofs/users/caogf/FBK_NUV_PDE/2019_7_8_230304/ch0.txt')
    SiPMAna.PlotWFs(ana, [0, 1, 2])
    ana.LowFilter()
    SiPMAna.PlotFilterWFs(ana, np.arange(0,10))
    ana.CalcQ()
    ana.PlotBaseline()
"""
