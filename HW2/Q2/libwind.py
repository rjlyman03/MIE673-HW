import numpy as np
import pandas as pd


def kaimal_f(f, U0, sigma, L):
    r""" 
    Single-sided Kaimal spectrum S_X(f)

    INPUTS:
     - f : frequency array [Hz]
     - U0: reference velocity [m/s]
     - sigma: standard deviation of wind [m/s]
     - L: turbulence length scale [m] 

    OUTPUTS:
      - S_X: Single-sided Kaimal spectrum as function of f
    """
    f = np.asarray(f)
    S = (4 * sigma**2 * L / U0) / (1 + 6 * f * L/U0) ** (5./3) 

    # NOTE: the Kaimal spectrum has a non-zero DC value S[0] = 4 * sigma**2 * L / U0
    #       You need this DC value to get \int S = sigma^2.
    #       But one shouldn't think that the mean of the time series is not zero.
    #       But it doesn't mean that the mean of the time series is not zero.
    # It's possible to artificially account for the mean (when generating a time series)
    # as follows:
    #    S[f == 0] = U0**2 * T # T is the total length in [s] 
    return S

def kaimal_omega(omega, U0, sigma, L):
    r""" 
    Single-sided Kaimal spectrum S_X(\omega)
         Sf(f) df = So(om)dom   =>   So(om) = Sf(f) / (2pi)
    """
    return  kaimal_f(omega/(2*np.pi), U0, sigma, L)/(2*np.pi)

def NTM_TI(WS, turb_class='A'):
    """ Normal Turbulence intensity model according to the IEC standards
    INPUTS: 
       -  WS: array of wind speed [m/s]
       -  turb_class: turbine class, "A", "B" or "C"
    OUTPUT:
       - TI: turbulence intensity at WS
    """
    WS = np.asarray(WS).astype(float)
    Iref={'A+':0.18, 'A':0.16, 'B':0.14, 'C':0.12}[turb_class]
    c    = 2.
    Vave = 10.
    sigma = c * Iref * (0.072 * (Vave / c + 3.) * (WS / c - 4.) + 10.)
    TI = np.full_like(WS, np.nan, dtype=float)
    mask = WS > 0
    TI[mask] = sigma[mask]/WS[mask]
    return TI

if __name__ == '__main__':
    # Small testing code
    import matplotlib.pyplot as plt

    WS = np.linspace(0, 10, 20)
    TI = NTM_TI(WS)

    fig, ax = plt.subplots(1, 1, figsize=(6.4,4.8))
    ax.plot(WS, TI, label='Class A')
    ax.set_xlabel('Wind speed [m/s]')
    ax.set_ylabel('TI [-]')
    plt.show()
