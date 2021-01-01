import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def readADC(infile):
    df = pd.read_csv(infile)
    ADC = np.reshape(df.to_numpy(), (1, -1))[0]
    fig, ax = plt.subplots()
    ax.plot(np.arange(len(ADC)), ADC)
    fig.savefig('ADC.png')
    return ADC

def readScope(infile):
    skiprows = [i for i in range(22)]
    for skip in range(10022, 110124):
        skiprows.append(skip)
    df = pd.read_fwf(infile, skiprows=skiprows)
    scope = df.iloc[:-1, 0].to_numpy()
    scope = scope.astype(np.float)
    fig, ax = plt.subplots()
    ax.plot(np.arange(len(scope)), scope)
    fig.savefig('scope.png')
    return scope

scope = readScope('/afs/ihep.ac.cn/users/w/weiwl/osc.txt')
from scipy import interpolate
f = interpolate.interp1d(np.arange(len(scope)) - 10, scope)
def func(x, a, b, c, d):
    return c*f(a*x+b)+d

def fit():
    ADC = readADC('/afs/ihep.ac.cn/users/w/weiwl/wave_11.txt')
    from scipy.optimize import curve_fit
    p0 = [2500/750., 0, 0, 0]
    popt, pcov = curve_fit(func, np.arange(len(ADC)), ADC, p0, bounds=((2500./800, -1, -np.inf, -np.inf), (2500./700, np.inf, np.inf, np.inf)) )
    print(popt, pcov)
    fig, axs = plt.subplots(3, figsize=(9,15))
    axs[0].plot(np.arange(len(ADC)), ADC, 'b-', label='ADC channel 11')
    axs[0].plot(np.arange(len(ADC)), func(np.arange(len(ADC)), *popt), 'g--',
            label='fit to scope' )
    axs[0].set_xlabel('Time (1/750 MHz-1)')
    axs[0].legend()
    axs[1].plot(np.arange(len(ADC)), ADC - func(np.arange(len(ADC)), *popt), color='blue')
    axs[1].set_xlabel('Time (1/750 MHz-1)')
    axs[1].set_ylabel('Residual of fit')
    axs[2].hist(ADC - func(np.arange(len(ADC)), *popt), bins=20, histtype='step', color='blue')
    axs[2].set_xlabel('Residual')
    plt.tight_layout()
    fig.savefig('fit_ADC_scope.pdf')
fit()
