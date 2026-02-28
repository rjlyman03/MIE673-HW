import numpy as np
from libwind import kaimal_omega
from libsignal import psd
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid

def spectrum2TimeSeries(S_Xw, omegas):
    PI = np.pi
    E = np.e
    N = len(omegas)
    dw = (omegas[1] - omegas[0])
    dt = 1/(2*PI*dw*N)
    t = range(0,N)*dt
    f = np.fft.rfftfreq(N, dt)
    N_plus = len(f)
    U = np.random.uniform(size = N_plus + 1)
    phi = np.random.uniform(size = N_plus + 1) * 2 * PI
    phi[0] = 0
    U[0] = 1/E
    X = [None] * (N_plus+1)
    A = [None] * (N_plus+1)
    for i in range(N_plus+1):
        if i==0:
            X[0] = N*np.sqrt(S_Xw[0]*dw)
        A[i] = np.sqrt(2*S_Xw[i]*dw)*np.sqrt(-1*np.log(U[i]))
        X[i] = 0.5*N*A[i]*np.exp(1j*phi[i])
    x = np.fft.irfft(X, N)
    return x, t

U0 = 11
sigma = 3
L = 100
f = 5*np.array(range(30))
omegas = 2*np.pi*f
S_test = kaimal_omega(omegas, U0, sigma, L)
x_test,t_test = spectrum2TimeSeries(S_test, omegas)
sigma_squared = np.var(x_test)
freqs, S_Xf = psd(t_test, x_test)
integral = trapezoid(freqs, S_Xf)
print(freqs)
print("Integral of PSD: ", integral)
print("Variance: ", sigma_squared)
plt.plot(t_test, x_test)
plt.xlabel('Time')
plt.ylabel('Turbulence')
plt.title('Kaimal Process Realization from Spectrum')
plt.show()


