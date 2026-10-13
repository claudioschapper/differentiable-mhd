from dmhd.spectral import diffspec, diffspec_gen
import numpy as np
import pytest
@pytest.mark.parametrize("N", [16, 32, 64])
@pytest.mark.parametrize("k", [1, 3, 5])


def test_diff(N, k):
    h = 2*np.pi/N
    x = h*np.arange(N)
    v1 = np.sin(k*x)
    w1 = diffspec(v1)
    assert np.max(np.abs(w1 - diffspec_gen(v1, 1))) < 1e-16 
    assert np.max(np.abs(w1-k*np.cos(k*x))) < 1e-12

@pytest.mark.parametrize("n", [16, 32, 64])
def test_nyquist(n):
    h = 2*np.pi/n
    x = h*np.arange(n)
    v2 = np.cos(n/2*x)
    w2 = diffspec(v2)
    assert np.max(np.abs(w2 - diffspec_gen(v2, 1))) < 1e-16 
    assert np.max(np.abs(w2)) < 1e-12


def test_evenm():
    N = 64
    h = 2*np.pi/N
    x = h*np.arange(N)
    v = (-1)**np.arange(N)
    w = diffspec_gen(v, 2)
    assert np.max(np.abs(w-(-((N/2)**2)*v))) < 1e-12








