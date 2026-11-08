from dmhd.spectral import diffspec, diffspec_gen, Grid, partial_x, partial_y, laplacian, poisson_inversion, jacobian
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

@pytest.mark.parametrize("L", [2*np.pi, 4*np.pi, 8*np.pi])
@pytest.mark.parametrize("N", [28, 30, 32])
def test_2d_derivatives(L, N):
    g = Grid.make(L,N)
    X, Y = g.coords
    f = np.sin(3*X) * np.cos(Y)
    dfdx = partial_x(f, g)
    dfdy = partial_y(f, g)
    lapf = laplacian(f,g)
    assert np.max(abs(3*np.cos(3*X)*np.cos(Y)-dfdx)) < 1e-11
    assert np.max(abs(-np.sin(3*X)*np.sin(Y)-dfdy)) < 1e-11
    assert np.max(abs(-10*np.sin(3*X)*np.cos(Y)-lapf)) < 1e-11


@pytest.mark.parametrize("L", [2*np.pi, 4*np.pi, 8*np.pi])
@pytest.mark.parametrize("N", [28, 30, 32])
def test_poisson_mean(L, N):
    c = 1.7
    g = Grid.make(L,N)
    X, Y = g.coords
    f = c + np.sin(2*np.pi*2*X/L) * np.cos(2*np.pi*3*Y/L) + 0.5*np.cos(2*np.pi*X/L)
    poiss_inv = poisson_inversion(f, g)
    lapf = laplacian(f,g)

    assert np.max(abs(-laplacian(poiss_inv, g)-(f-c))) < 1e-12
    assert np.max(abs(-poisson_inversion(lapf, g)-(f-c))) < 1e-12
    assert abs(np.mean(poiss_inv)) < 1e-12


def test_imaginary():
    L = 2*np.pi
    N = 32
    g = Grid.make(L,N)
    X, Y = g.coords
    f =  np.random.default_rng(0).standard_normal((N, N))
    dfdx = partial_x(f, g)
    dfdy = partial_y(f, g)
    poiss_inv = poisson_inversion(f, g)
    lapf = laplacian(f,g)

    assert np.max(abs(dfdx.imag)) < 1e-12
    assert np.max(abs(dfdy.imag)) < 1e-12
    assert np.max(abs(lapf.imag)) < 1e-12
    assert np.max(abs(poiss_inv.imag)) < 1e-12

@pytest.mark.parametrize("N", [14, 28, 30, 32])
def test_mask_count(N):
    L = 2*np.pi
    m = (N-1)//3
    g = Grid.make(L, N)
    mask = g.mask_twothirds
    count = 0
    for i in range(N):
        for j in range(N):
            if mask[i][j]:
                count += 1
    assert count == (2*m + 1)**2

@pytest.mark.parametrize("N", [14, 28, 30, 32])
def test_mask_symmetry(N):
    L = 2*np.pi
    g = Grid.make(L, N)
    mask = g.mask_twothirds
    sym = mask[np.ix_(-np.arange(N)%N, -np.arange(N)%N )]
    assert np.array_equal(sym, mask)

def test_jacobian():
    L = 2*np.pi
    N = 24
    g = Grid.make(L,N)
    X, Y = g.coords
    psi = np.sin(2*X)*np.cos(Y)
    omega = np.cos(X)*np.sin(3*Y)
    J_ref = 1/8*(5*(np.cos(X-2*Y)+np.cos(3*X+4*Y)+np.cos(3*X-4*Y)+np.cos(X+2*Y))+
                        7*(np.cos(3*X-2*Y)+np.cos(3*X+2*Y)+np.cos(X+4*Y)+np.cos(X-4*Y)))
    J_hat = jacobian(psi, omega, g)


    assert np.max(abs(np.fft.ifft2(J_hat).real-J_ref)) < 1e-12
    assert np.max(abs(np.fft.ifft2(J_hat).imag)) < 1e-13
    assert np.max(abs(J_ref-(6*np.cos(2*X)*np.cos(Y)*np.cos(X)*np.cos(3*Y)-np.sin(2*X)*np.sin(Y)*np.sin(X)*np.sin(3*Y)))) < 1e-13
    assert np.max(abs(J_hat+jacobian(omega, psi, g))) < 1e-13
