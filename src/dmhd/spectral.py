import numpy as np

"""assumes a 2π-periodic domain, even N"""
def wavenumbers(N):
    return np.fft.fftfreq(N, d=1/N)


def diffspec(v):
    N = len(v)
    assert N % 2 == 0
    wavenum = wavenumbers(N)
    v_hat = np.fft.fft(v)
    w_hat = 1j*wavenum*v_hat
    w_hat[N//2] = 0
    return np.fft.ifft(w_hat).real


