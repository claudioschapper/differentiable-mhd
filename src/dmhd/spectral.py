import numpy as np


def wavenumbers(N):
    """assumes a 2π-periodic domain, even N"""
    return np.fft.fftfreq(N, d=1/N)


def diffspec(v):
    N = len(v)
    assert N % 2 == 0
    wavenum = wavenumbers(N)
    v_hat = np.fft.fft(v)
    w_hat = 1j*wavenum*v_hat
    w_hat[N//2] = 0
    return np.fft.ifft(w_hat).real

def diffspec_gen(v, m):
    N = len(v)
    assert N % 2 == 0
    wavenum = wavenumbers(N)
    v_hat = np.fft.fft(v)
    w_hat = ((1j*wavenum)**m)*v_hat
    # Need to zero nyquist mode for odd m, N/2 and -N/2 wavemodes have opposite derivatives
    if m%2 == 1:
        w_hat[N//2] = 0

    return np.fft.ifft(w_hat).real

