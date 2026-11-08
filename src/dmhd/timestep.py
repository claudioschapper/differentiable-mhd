import numpy as np
def rk4_step(f, t, y, h):
    k1 = f(t, y)
    k2 = f(t +h /2, y + h*k1/2)
    k3 = f(t + h/2, y + h*k2/2)
    k4 = f(t + h, y + h*k3)
    return y + h/6*(k1 + 2*k2 +2*k3 + k4)

def if_step(u_hat, L, N, dt):
    """L and N use the right hand convention u^hat_t = Lu + N(u), N takes a physical space field and returns a fourier-space array
    Keep in mind that ifft can return imaginary values (for diagnostics that assume that u is real (consider rfft and irfft))"""
    E = np.exp(L*dt/2)
    # a1 to a4 are fourier space arrays
    a1 = N(np.fft.ifftn(u_hat))
    a2 = N(np.fft.ifftn(E*(u_hat+dt/2*a1)))
    a3 = N(np.fft.ifftn(E*u_hat+dt/2*a2))
    a4 = N(np.fft.ifftn((E**2)*u_hat+dt*E*a3))
    return dt/6*((E**2)*a1+2*E*a2+2*E*a3+a4) + (E**2)*u_hat


