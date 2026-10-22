import numpy as np
from dmhd.solvers.burgers import solve_burgers, N, u_exact
from dmhd.spectral import wavenumbers


def test_u_exact():
    t = 0.3
    x = 1
    v = 0.1
    a = 0.5
    h1 = 1e-4
    h2 = 1e-3
    u = u_exact(x, t, v, a)
    u_t = (u_exact(x, t+h1, v, a) - u_exact(x, t-h1, v, a))/(2*h1)
    u_x = (u_exact(x+h1, t, v, a) - u_exact(x-h1, t, v, a))/(2*h1)
    u_xx = (u_exact(x+h2, t, v, a) - 2*u_exact(x, t, v, a) + u_exact(x-h2, t, v, a))/(h2)**2
    res = u_t + u*u_x -v*u_xx
    assert np.abs(res) < 1e-8

def test_N():
    n = 64
    h = 2*np.pi/n
    x = h*np.arange(n)
    test = np.fft.ifft(N(n, np.sin(x)))
    analytic = -1/2*np.sin(2*x)
    assert np.max(np.abs(test.real-analytic)) < 1e-14
    assert np.max(np.abs(test.imag)) < 1e-14

def test_N_mean():
    rng = np.random.default_rng(0)
    n = 64
    u = rng.standard_normal(n) + 1
    assert N(n,u)[0] == 0

def test_exact_spectrum():
    n = 64
    h = 2*np.pi/n
    nu = 0.1
    a = 0.5
    x = h*np.arange(n)
    u_hat = np.fft.fft(u_exact(x, 0, nu, a))
    assert np.abs(u_hat[0]) < 1e-14
    assert np.abs(u_hat[n//3])/np.max(np.abs(u_hat)) < 1e-11

def test_aliasing():
    n = 64
    h = 2*np.pi/n 
    x = h*np.arange(n)
    u = np.cos(21*x)
    u_mask = np.abs(N(n,u))
    wavenums = wavenumbers(n)
    u_no_mask = (-1j*wavenums*np.fft.fft((u.real)**2)/2)
    assert np.max(u_mask) < 1e-10
    assert np.abs(u_no_mask[22]) > 1
    assert np.all(u_mask[22:43] == 0)

def test_dynamics():
    n = 96
    nu = 0.1
    a = 0.5
    T = 1
    n_steps = 200
    h = 2*np.pi/n 
    x = h*np.arange(n)
    u0 = u_exact(x, 0, nu, a)
    u = solve_burgers(u0, nu, T, n_steps)
    assert np.max(np.abs(u-u_exact(x, T, nu, a))) < 1e-13
    u_2 = solve_burgers(u0, nu, T, 2*n_steps)
    assert np.max(np.abs(u_2-u_exact(x, T, nu, a))) < 1e-13

def test_spatial():
    n_steps = 200          # fixed so time error stays ~1e-15
    nu = 0.1
    a = 0.5
    T = 1
    ns = np.arange(8, 104, 8)

    errs = np.zeros(len(ns))
    for i, n in enumerate(ns):
        x = 2 * np.pi / n * np.arange(n)
        u0 = u_exact(x, 0, nu, a)
        u = solve_burgers(u0, nu, T, n_steps)
        errs[i] = np.max(np.abs(u - u_exact(x, T, nu, a)))
    assert errs[ns == 64][0] < 1e-13
    sel = (ns >= 16) & (errs > 1e-12)
    assert sel.sum() >= 4
    rate = -np.polyfit(ns[sel], np.log(errs[sel]), 1)[0]
    y0_T = np.arccosh(np.exp(nu * T) / a)
    # Prediction: if u is analytic in the strip |Im x| < y0, its Fourier
    # coefficients decay like exp(-y0 |k|). Here the singularity is the zero of
    # theta = 1 + a e^{-nu T} cos x at x = pi + i y0, so cosh(y0) = e^{nu T} / a.
    # The 2/3 rule keeps only |k| <= n/3, so the truncation error scales as
    # exp(-y0 n / 3). The error is measured at T, so y0 is evaluated at T
    # (a at t = 0 would give 0.439, 8% off). Measured rate 0.478 vs prediction
    # 0.477; the 5% tolerance allows for the k r^k prefactor and fit window.
    assert np.abs((y0_T/3-rate)/rate) < 5e-2
