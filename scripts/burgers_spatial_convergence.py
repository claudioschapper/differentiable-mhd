import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from dmhd.solvers.burgers import solve_burgers, u_exact

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

print("n:", ns)
print("errors:", errs)

# fit ln(err) vs n on the straight part: past the first points, above the floor
sel = (ns >= 16) & (errs > 1e-12)
rate = -np.polyfit(ns[sel], np.log(errs[sel]), 1)[0]
y0_T = np.arccosh(np.exp(nu * T) / a)

print(f"fitted decay rate per n (ln): {rate:.3f}   y0(T)/3: {y0_T/3:.3f}   y0(T)/2: {y0_T/2:.3f}")
fig, ax = plt.subplots()

n0, e0 = ns[2], errs[2]
ax.semilogy(ns, e0 * np.exp(-y0_T / 3 * (ns - n0)), "k--", label=r"$e^{-y_0(T)\,n/3}$ (2/3 rule)")
ax.semilogy(ns, e0 * np.exp(-y0_T / 2 * (ns - n0)), "k:", label=r"$e^{-y_0(T)\,n/2}$ (full band)")
ax.semilogy(ns, np.maximum(errs, 1e-17), "o-", label="IF-RK4 + 2/3 rule")
n0, e0 = ns[2], errs[2]
ax.set_ylim(1e-17, 1)
ax.set_xlabel("n (grid points)")
ax.set_ylabel("max error at T")
ax.set_title(f"Burgers spatial convergence (nu={nu}, a={a}, T={T}, {n_steps} steps)")
ax.legend()
out = Path(__file__).resolve().parent.parent / "figures" / "burgers_spatial_convergence.png"
fig.savefig(out, dpi=150, bbox_inches="tight")