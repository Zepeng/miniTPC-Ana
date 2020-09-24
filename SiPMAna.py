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
    channel = 0
    wfpath = ''
    qdcpath = ''
    wflen = 0
    bSaveFig = False
    bDebug = False
    def __init__(self, filedir, wflen, ch=0, savefig=False, isDebug=False):
        self.waveforms = 0 # waveforms
        self.filteredwfs = 0 # waveform after software filter
        self.charges = 0 # a software QDC
        self.qdc = 0 # QDC saved from TD
        self.wflen = wflen #length of waveform used in TD5751
        self.bSaveFig = savefig # set to True to save figures
        self.channel = ch
        self.bDebug = isDebug
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
        self.LoadWFs()

    def LoadWFs(self):
        #four channels are saved in four different files, \
        #and each line is a single waveform in the file
        with open(self.wfpath) as f:
            i = 0
            wfs = []
            for line in f:
                wf = []
                for item in line.split():
                    wf.append(float(item))
                wfs.append(wf)
                i+=1
                if self.bDebug and i >= 5000:
                    print('Loading 5000 waveforms in debug mode.')
                    break
                if i % 1000 ==0:
                    print('Loaded %d of waveforms from file.' % i)
            self.waveforms = np.array(wfs)
            print('Loaded %d waveforms from %s' % (i, self.wfpath))

    def LoadQDC(self):
        #four channels are saved in a single file but there are some garbage\
        #in the QDC file. Select the useful data by line number
        with open(self.qdcpath) as f:
            i = 0
            qdc = []
            for line in f:
                qdc = []
                for item in line.split():
                    #There are some inf in the data file. I don't know the reason\
                    #just skip that part of data.
                    if float(item) < np.Inf:
                        qdc.append(float(item))
                fig, ax = plt.subplots()
                ax.hist(np.array(qdc), bins=np.linspace(np.mean(qdc)-2*np.std(qdc),\
                        np.mean(qdc) + 2*np.std(qdc), 40),histtype='step', color='blue')
                ax.set_xlabel('QDC')
                ax.set_ylabel('Events')
                fig.savefig('qdc_%d.pdf' % i)
                i += 1

    def PlotWFs(self, entries = [0]):
        for entry in entries:
            fig, ax = plt.subplots()
            if entry > np.size(self.waveforms, 0):
                print('Entry not found in saved waveforms')
                return
            ax.plot(np.arange(len(self.waveforms[entry])), self.waveforms[entry])
            ax.set_xlabel('Time (ns)')
            ax.set_ylabel('ADC (mV)')
            fig.savefig('wfs_%d_%d.png' % (self.channel, entry) )

    def PlotBaseline(self):
        """
        correct the baseline by the mean value
        and plot the baseline
        """
        baseline = []
        for wf in self.filteredwfs:
            baseline.append(np.array(wf[100:2100]) - np.mean(wf[100:2100]))
        fig, ax = plt.subplots()
        ax.hist(np.array(baseline).flatten(), bins=np.linspace(np.mean(baseline) - 4*np.std(baseline)\
                , np.mean(baseline)+4*np.std(baseline), 50), histtype='step', color='blue')
        from scipy.stats import norm
        (mu, sigma) = norm.fit(baseline)
        drawtext = 'Mean: %.3f\n $\sigma$: %.3f' % (mu, sigma)
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        ax.text(xmax - 0.3*(xmax - xmin), ymax - 0.2*(ymax-ymin), drawtext)
        print(mu, sigma)
        ax.set_xlabel('Baseline subtracted mean (mV)')
        ax.set_ylabel('N points')
        fig.savefig('baseline.png')

    def EventBaseline(self):
        """
        correct the baseline by the mean value
        and plot the baseline in each event
        """
        mean = []
        sigma = []
        for wf in self.filteredwfs:
            baseline = np.array(wf[100:2100]) - np.mean(wf[100:2100])
            mean.append(np.mean(baseline))
            sigma.append(np.std(baseline))
        fig, ax = plt.subplots()
        ax.hist(np.array(sigma), bins = 30, histtype='step', color='blue')
        ax.set_xlabel('$\sigma$ of baseline (mV)')
        ax.set_ylabel('Events')
        fig.savefig('event_baseline.pdf')

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
            ax.set_xlim(2500, 3200)
            ax.set_xlabel('Time (ns)')
            ax.set_ylabel('ADC (mV)')
            fig.savefig('filteredwfs_%d.png' % entry)

    def CalcQ(self):
        Qs = []
        for wf in self.filteredwfs:
            baseline = np.mean(wf[100:2000])
            q = np.sum(np.array(wf[2500:2750]) - baseline)
            Qs.append(q)
        self.charges = np.array(Qs)
        fig, ax = plt.subplots()
        ax.hist(Qs, bins=100, histtype='step')
        ax.set_xlabel('Software QDC')
        ax.set_ylabel('Events')
        #ax.set_yscale('log')
        fig.savefig('qdc.png')


if __name__ == '__main__':
    ana = SiPMAna(filedir='/junofs/users/weiwl/sipmtest/105-47.3/', ch=1, wflen=5005, isDebug=True)
    ana.LowFilter()
    SiPMAna.PlotFilterWFs(ana, [0, 1, 2, 3, 4, 5])
    SiPMAna.EventBaseline(ana)
    SiPMAna.CalcQ(ana)
