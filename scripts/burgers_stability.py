import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from dmhd.solvers.burgers import solve_burgers, u_exact, N
from dmhd.timestep import rk4_step
from dmhd.spectral import wavenumbers

n_steps = 20   
nu = 0.1
a = 0.5
T = 1
n = 64
x = 2 * np.pi / n * np.arange(n)
dt = T/n_steps
wavenums = wavenumbers(n)

def rhs(t, u_hat):
    L = -nu*(wavenums**2)*u_hat
    return L + N(n, np.fft.ifft(u_hat))
u0 = u_exact(x, 0, nu, a)
u_rk = u0
u_rk_hat = np.fft.fft(u0)
ts = dt * np.arange(n_steps + 1)                 # reach T = 1
max_rk, errs_rk, max_if, errs_if = (np.zeros(n_steps + 1) for _ in range(4))
ks = np.array([n // 2, n // 2 - 1, n // 2 - 2])      # k = 32, 31, 30
modes_rk = np.full((len(ks), n_steps + 1), np.nan)
modes_if = np.full((len(ks), n_steps + 1), np.nan)
blown = False
for i, t in enumerate(ts):
    u_if = solve_burgers(u0, nu, t, i) if t > 0 else u0
    u_if_hat = np.fft.fft(u_if)
    modes_if[:, i] = np.abs(u_if_hat[ks])
    u_exac = u_exact(x, t, nu, a)
    if not blown:
        m = np.max(np.abs(u_rk_hat))
        if m > 1e10:                 # produced by the previous step: discard it
            blown = True
    if blown:
        max_rk[i] = errs_rk[i] = np.nan
    else:
        max_rk[i], errs_rk[i] = m, np.max(np.abs(u_rk - u_exac))
        modes_rk[:, i] = np.abs(u_rk_hat[ks])
    max_if[i], errs_if[i] = np.max(np.abs(u_if_hat)), np.max(np.abs(u_if - u_exac))
    if not blown:
        u_rk_hat = rk4_step(rhs, dt * i, u_rk_hat, dt)
        u_rk = np.fft.ifft(u_rk_hat)


k_max = np.max(np.abs(wavenumbers(n)))
z = -nu * k_max**2 * dt
R = abs(1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24)
print(f"nu k_max^2 dt = {-z:.2f}, predicted |R| per step = {R:.2f}")

# measured growth on the clean linear part: above the roundoff seed,
# below the size where the garbage starts feeding the nonlinearity
grow = np.isfinite(errs_rk) & (errs_rk > 100 * errs_if) & (errs_rk < 1e-2)
idx = np.flatnonzero(grow)
if len(idx) >= 3:
    g = np.exp(np.polyfit(ts[idx], np.log(errs_rk[idx]), 1)[0] * dt)
    print(f"measured RK4 error growth per step: {g:.2f}")

def clean(y, floor=1e-17):
    y = np.asarray(y, float).copy()
    y[~np.isfinite(y)] = np.nan
    return np.maximum(y, floor)  # t = 0 error is 0; keep it on a log axis

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4), sharex=True)

ax1.semilogy(ts, clean(max_if), "o-", label="IF-RK4")
ax1.semilogy(ts, clean(max_rk), "s-", label="plain RK4")
ax1.set_xlabel("t")
ax1.set_ylabel(r"max $|\hat u|$")
ax1.set_title("solution amplitude")
ax1.legend()

ax2.semilogy(ts, clean(errs_if), "o-", label="IF-RK4")
ax2.semilogy(ts, clean(errs_rk), "s-", label="plain RK4")
if len(idx):
    j = idx[0]
    ax2.semilogy(ts[j:], errs_rk[j] * R ** ((ts[j:] - ts[j]) / dt), "k--",
                 label=rf"$|R(z)|^{{t/\Delta t}}$, $|R|={R:.1f}$")
ax2.set_ylim(1e-17, 1e12)
ax2.set_xlabel("t")
ax2.set_ylabel("max error vs exact")
ax2.set_title(rf"Burgers, n={n}, $\nu$={nu}, $\nu k_{{max}}^2\Delta t$={-z:.2f}")
ax2.legend()

out = Path(__file__).resolve().parent.parent / "figures" / "burgers_stability.png"
fig.savefig(out, dpi=150, bbox_inches="tight")

def R_of(z):
    return abs(1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24)

fig, ax = plt.subplots(figsize=(6.5, 4.5))
for j, k in enumerate(ks):
    Rk = R_of(-nu * k**2 * dt)
    ok = np.isfinite(modes_rk[j])
    g = np.exp(np.polyfit(ts[ok], np.log(modes_rk[j][ok]), 1)[0] * dt)
    ax.semilogy(ts, modes_rk[j], "s", color=f"C{j}",
                label=rf"RK4, $k={k}$: measured $\times${g:.2f}, predicted $\times${Rk:.2f}")
    ax.semilogy(ts, modes_rk[j, 0] * Rk ** (ts / dt), "--", color=f"C{j}", lw=1)
    print(f"k={k}: measured growth {g:.3f}, predicted |R| {Rk:.3f}")
ax.semilogy(ts, np.maximum(modes_if.max(axis=0), 1e-19), "o-", color="0.4",
            label=r"IF-RK4, $k=30,31,32$ (max)")
ax.set_ylim(1e-19, 1e6)
ax.set_xlabel("t")
ax.set_ylabel(r"$|\hat u_k|$")
ax.set_title(rf"Grid-scale modes, n={n}, $\nu$={nu}, $\nu k_{{max}}^2\Delta t$={-z:.2f}"
             "\n(dashed: predicted $|R|^{t/\\Delta t}$)", fontsize=10)
ax.legend(fontsize=8, loc="upper left")
out = Path(__file__).resolve().parent.parent / "figures" / "burgers_stability_modes.png"
fig.savefig(out, dpi=150, bbox_inches="tight")