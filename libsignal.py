import numpy as np
import pandas as pd

def psd(t, y, averaging=None, nperseg=None):
    """ Perform single-sided PSD of signal y 
    INPUTS:
     - t : time vector
     - y : signal vector
     - averaging: None or 'welch' (for a smoother spectrum)
     - nperseg  : argument for averaging with the Welch method, typically 2^10 or so
    """
    t = np.asarray(t)
    y = np.asarray(y)
    N  = len(y)
    dt = (t[-1]-t[0])/(N-1)
    T  = N*dt
    fs = 1/dt
    if averaging is None:
        Y          = np.fft.rfft(y)
        freqs      = np.fft.rfftfreq(len(y), dt) # [Hz]
        df         = freqs[1]-freqs[0]
        A          = np.abs(Y) / N  # double-sided "amplitudes"
        S_Xf       = A**2 *N*dt     # double-sided PSD
        S_Xf[1:-1] *= 2        # single-sided PSD 
        #A[1:-1]    *= 2       # single-sided "amplitudes"
        #phi         = np.angle(Y)  # Phases
    elif averaging=='welch':
        from scipy.signal import welch
        def fnextpow2(x):
            return 
        if nperseg is None:
            nFFTAll = 2 ** int(np.ceil(np.log2(max(N, 1))))
            nperseg = nFFTAll // 2
            nperseg = max(128, min(nperseg, N // 2))       # reasonable bounds
        freqs, S_Xf = welch(y, fs, nperseg=nperseg)
    else:
        raise NotImplementedError(averaging)
    return freqs, S_Xf


def autoCorrCoeff(x, nMax=None, dt=1, method='corrcoef'):
    """ 
    Compute auto correlation of a signal, rho_xx(tau)
    INPUTS:
      - x     : input signal
      - nMax  : number of values to return (typically nMax << len(x))
      - dt    : optional time step to scale the delay tau
      - method: optional method, between 'direct-sum', or 'corrcoeff'
    OUTPUTS:
      - rho: auto-correlation coefficient vector
      - tau: delay vector
    """
    x    = np.asarray(x)
    x    = x.copy() - np.mean(x)
    var  = np.var(x)
    n    = len(x)
    if nMax is None:
        nMax = n
    rvec = np.arange(0,nMax)
    if method=='direct-sum':
        # Using actual definition of the correlation coefficient, with a direct sum (expensive)
        rho    = np.zeros(nMax)
        rho[0] =1
        for i, nDelay in enumerate(rvec[1:]):
            rho[i+1] = np.mean( x[0:-nDelay] * x[nDelay:] ) / var

    elif method=='corrcoef':
        rho = np.array([1]+[np.corrcoef(x[:-i], x[i:])[0,1]  for i in range(1, nMax)])

    elif method=='correlate':
        rho = np.correlate(x, x, mode='full')[-n:] / (var * n)
        rho = rho[:nMax]
    else:
        raise NotImplementedError(method)

    tau = rvec*dt
    return rho, tau


def bin_DF(df, colBin, xBins=None, stats='mean', dropObjects=True):
    """ 
    Perform bin averaging of a dataframe
    INPUTS:
      - df   : pandas dataframe
      - colBin: column name (string) of the dataframe, used for binning 
      - xBins: end points delimiting the bins, array of ascending x values
      - stats: statistics to compute (either mean or std)
      - dropObjects: columns that are objects (e.g. strings or timestamps) are discarded

    OUTPUTS:
       binned dataframe, with additional column 'Counts' for the number 
    """
    if colBin not in df.columns.values:
        raise Exception('The column `{}` does not appear to be in the dataframe'.format(colBin))

    if xBins is None:
        xBins = np.linspace(*np.percentile(df[colBin].dropna(), [0.5, 99.5]), 21)
    else:
        xBins = np.asarray(xBins)

    if dropObjects:
        non_numeric_cols = df.select_dtypes(exclude=['float64', 'int64', 'Int64', 'Float64']).columns
        df[non_numeric_cols] = np.nan

    binMid  = (xBins[:-1]+xBins[1:])/2 # Geometric midpoints

    df['Bin'] = pd.cut(df[colBin], bins=xBins, labels=binMid ) # Adding a column that has bin attribute
    if stats=='mean':
        df2       = df.groupby('Bin', observed=False).mean(numeric_only=True) # Average by bin
    elif stats=='std':
        df2       = df.groupby('Bin', observed=False).std(numeric_only=True)  # Std by bin
    else:
        raise ValueError("stats must be 'mean' or 'std'")
    # also counting
    df['Counts'] = 1
    dfCount=df[['Counts','Bin']].groupby('Bin', observed=False).sum()
    df2['Counts'] = dfCount['Counts']
    # Ensure all bins are present, empty ones get NaN
    df2 = df2.reindex(binMid)
    # Add Boundary and center
    df2 = df2.reset_index(names='binMid') # binMid is now a columns index → column
    df2['binStart'] = xBins[:-1]
    df2['binEnd']   = xBins[1:]
    bin_means = df.groupby('Bin', observed=False)[colBin].mean()
    df2['binMean'] = bin_means.reindex(binMid).values # Note: same as df2[colBin] when stats=mean
    # Optional: sort just to be extra safe
    df2 = df2.sort_values('binMid').reset_index(drop=True)
    return df2




if __name__ == '__main__':
    # Small testing code
    import matplotlib.pyplot as plt

    t = np.linspace(0, 40*np.pi, 1000)
    y =  10 * np.sin(2*np.pi*t)
    f, S = psd(t, y)
    plt.loglog(f, S)

    plt.show()
