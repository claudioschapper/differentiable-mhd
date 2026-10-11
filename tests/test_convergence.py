from dmhd.spectral import diffspec
import numpy as np
import pytest



def test_convergence():
    Ns = { N:0 for N in range(4,50,2)}

    for N in Ns:
        h = 2*np.pi/N
        x = h*np.arange(N)
        v = np.exp(np.sin(x))
        w = diffspec(v)
        Ns[N] = np.max(np.abs(w-np.cos(x)*v))
        if N >= 6 and Ns[N] >= 1e-13:
            assert Ns[N] <= Ns[N-2], f"error rose between N = {N-2} and N = {N}"
        if N >= 24:
            assert Ns[N] < 1e-12
        