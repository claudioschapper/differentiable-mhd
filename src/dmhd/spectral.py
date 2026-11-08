import numpy as np
from dataclasses import dataclass

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

@dataclass(frozen=True, eq=False)
class Grid:
    L: float
    N: int
    x: np.ndarray
    waveseven: np.ndarray
    wavesodd: np.ndarray
    coords: tuple
    kx_even_der: np.ndarray
    ky_even_der: np.ndarray
    kx_odd_der: np.ndarray
    ky_odd_der: np.ndarray
    mask_twothirds: np.ndarray
    @classmethod 
    def make(cls, L, N):
        assert N%2 == 0
        x = L*np.arange(N)/N
        x.flags.writeable = False
        # need to separate into two wavenumber arrays, as we need to handle the nyquist mode which should be zero for odd derivatives
        waveseven = np.fft.fftfreq(N, d=L/N)*2*np.pi
        waveseven.flags.writeable = False
        wavesodd = waveseven.copy()
        wavesodd[N//2] = 0
        wavesodd.flags.writeable = False
        xc, yc = np.meshgrid(x, x, indexing="ij")
        xc.flags.writeable = False
        yc.flags.writeable = False
        kx_even_der = np.meshgrid(waveseven, np.zeros(N), indexing="ij")[0]
        kx_even_der.flags.writeable = False
        ky_even_der = np.meshgrid(np.zeros(N), waveseven, indexing="ij")[1]
        ky_even_der.flags.writeable = False
        kx_odd_der = np.meshgrid(wavesodd, np.zeros(N), indexing="ij")[0]
        kx_odd_der.flags.writeable = False
        ky_odd_der = np.meshgrid(np.zeros(N), wavesodd, indexing="ij")[1]
        ky_odd_der.flags.writeable = False
        modes = np.fft.fftfreq(N, d=1/N)
        mx = np.meshgrid(modes, np.zeros(N), indexing="ij")[0]
        my = np.meshgrid(np.zeros(N), modes, indexing="ij")[1]
        M = (N-1)//3
        mask_twothirds = np.logical_and(np.abs(mx) <=M, np.abs(my) <= M )
        mask_twothirds.flags.writeable = False
        return cls(L = L, N = N, x = x, waveseven = waveseven, wavesodd = wavesodd, coords = (xc, yc),
                    kx_even_der = kx_even_der, ky_even_der = ky_even_der, kx_odd_der = kx_odd_der, ky_odd_der = ky_odd_der,
                    mask_twothirds = mask_twothirds)

def partial_x(u, grid):
    u_hat = np.fft.fft2(u)
    part_x_u_hat = 1j*grid.kx_odd_der*u_hat
    return np.fft.ifft2(part_x_u_hat)

def partial_y(u, grid):
    u_hat = np.fft.fft2(u)
    part_y_u_hat = 1j*grid.ky_odd_der*u_hat
    return np.fft.ifft2(part_y_u_hat)

def laplacian(u, grid):
    u_hat = np.fft.fft2(u)
    laplacian_u_hat = -(grid.kx_even_der**2 + grid.ky_even_der**2)*u_hat
    return np.fft.ifft2(laplacian_u_hat)

def poisson_inversion(u, grid):
    u_hat = np.fft.fft2(u)
    kx = grid.kx_even_der.copy()
    kx[0][0] = 1
    ky = grid.ky_even_der.copy()
    ky[0][0] = 1
    poi_inv_u_hat = u_hat/(kx**2+ky**2)
    poi_inv_u_hat[0][0] = 0
    return np.fft.ifft2(poi_inv_u_hat)
