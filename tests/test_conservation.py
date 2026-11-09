from dmhd.spectral import Grid, partial_x, partial_y, laplacian, poisson_inversion, jacobian
import numpy as np
import pytest
import dataclasses 

@pytest.mark.parametrize("N", [24, 30, 32])
@pytest.mark.parametrize("psi_seed", [0, 1, 2])
@pytest.mark.parametrize("omega_seed", [3, 4, 5])
def test_jac_conservation(N, psi_seed, omega_seed):
    L = 2*np.pi

    g = Grid.make(L,N)
    psi =  np.random.default_rng(psi_seed).standard_normal((N, N))
    psi_hat = np.fft.fft2(psi)*g.mask_twothirds
    psi = np.fft.ifft2(psi_hat)
    assert np.max(np.abs(psi.imag)) < 1e-13
    psi = psi.real
    rms_psi = np.sqrt(np.mean(psi**2))
    rms_grad_psi = np.sqrt(np.mean(((partial_x(psi,g).real)**2 + (partial_y(psi,g).real)**2)))

    omega =  np.random.default_rng(omega_seed).standard_normal((N, N))
    omega_hat = np.fft.fft2(omega)*g.mask_twothirds
    omega = np.fft.ifft2(omega_hat)
    assert np.max(np.abs(omega.imag)) < 1e-13
    omega = omega.real
    rms_omega = np.sqrt(np.mean(omega**2))
    rms_grad_omega = np.sqrt(np.mean(((partial_x(omega,g).real)**2 + (partial_y(omega,g).real)**2)))

    jac = np.fft.ifft2(jacobian(psi, omega, g)).real
    assert np.max(np.abs(jac.imag)) < 1e-13
    rms_jac = np.sqrt(np.mean(jac**2))

    assert np.abs(np.mean(jac*psi)/(rms_psi*rms_jac)) < 1e-16
    assert np.abs(np.mean(jac*omega)/(rms_omega*rms_jac)) < 1e-16
    assert np.abs(rms_jac/(rms_grad_omega*rms_grad_psi)) > 0.1

@pytest.mark.parametrize("N", [24, 30])
@pytest.mark.parametrize("psi_seed", [0, 1, 2])
@pytest.mark.parametrize("omega_seed", [3, 4, 5])

def test_jac_failure(N, psi_seed, omega_seed):
    L = 2*np.pi

    g = Grid.make(L,N)
    modes = np.fft.fftfreq(N, d=1/N)
    mx = np.meshgrid(modes, np.zeros(N), indexing="ij")[0]
    my = np.meshgrid(np.zeros(N), modes, indexing="ij")[1]
    M = (N)//3
    bad_mask = np.logical_and(np.abs(mx) <=M, np.abs(my) <= M )
    g = dataclasses.replace(g, mask_twothirds=bad_mask)

    psi =  np.random.default_rng(psi_seed).standard_normal((N, N))
    psi_hat = np.fft.fft2(psi)*g.mask_twothirds
    psi = np.fft.ifft2(psi_hat)
    psi = psi.real
    rms_psi = np.sqrt(np.mean(psi**2))
    rms_grad_psi = np.sqrt(np.mean(((partial_x(psi,g).real)**2 + (partial_y(psi,g).real)**2)))

    omega =  np.random.default_rng(omega_seed).standard_normal((N, N))
    omega_hat = np.fft.fft2(omega)*g.mask_twothirds
    omega = np.fft.ifft2(omega_hat)
    omega = omega.real
    rms_omega = np.sqrt(np.mean(omega**2))
    rms_grad_omega = np.sqrt(np.mean(((partial_x(omega,g).real)**2 + (partial_y(omega,g).real)**2)))

    jac = np.fft.ifft2(jacobian(psi, omega, g)).real
    rms_jac = np.sqrt(np.mean(jac**2))

    assert np.abs(np.mean(jac*psi)/(rms_psi*rms_jac)) > 1e-4
    assert np.abs(np.mean(jac*omega)/(rms_omega*rms_jac)) > 1e-4
