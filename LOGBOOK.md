# Logbook

## 2026-10-10 / 2026-10-11
Goal: 1D heat equation, pseudo-spectral, vs analytic solution.
- Spectral first derivative (`diffspec`) done: matches analytic to ~1e-14, Nyquist mode zeroed.
- Bugs hit: `N/2` float used as index, missing factor of i, negative wavenumbers had wrong sign.
- Tests: exactness over N and k, Nyquist, convergence on exp(sin x). All green.
- Convergence plot of exp(sin x) saved to figures/convergence_exp_sin.png.
NEXT: second derivative multiplier (-k^2, decide Nyquist), RK4 in timestep.py, heat equation vs exp(-nu k^2 t) sin(kx).

## 2026-10-13
Goal: 1D heat equation, pseudo-spectral + RK4, verified against the analytic solution.
- `diffspec_gen(v, m)`: general order-m derivative. Nyquist slot zeroed for odd m only. For even m, (ik)^m is even in k, so +N/2 and -N/2 agree and nothing is biased; keeping the mode is correct. Test: (-1)^j under m=2 returns -(N/2)^2 v.
- `rk4_step(f, t, y, h)` in timestep.py. Verified on y' = lam*y (lam = -1..-5): fitted order 4.02-4.12. On y' = t^3 the error is at roundoff (<1e-16), so the stage times are right.
- `solve_heat` (evolves u_hat, diagonal RHS -nu k^2 u_hat): matches e^(-nu t) sin x + 0.5 e^(-9 nu t) sin 3x below 1e-12 (n_steps=1000, N=32). Mean conserved, imaginary part of the final ifft ~5e-17. Fitted time order 4.06.
- Bug hit: order fit gave slope -4 because I regressed against n_steps/N instead of h = T/n_steps. A slope of exactly -p means the x-axis is inverse to the step.
- Tooling lesson: editing files on disk while VS Code held unsaved buffers caused a save conflict and an overwritten LOGBOOK. Save/close tabs first.
- RK4 stability verified on heat (notebooks/2026-10-13-heat-stability.ipynb): measured per-step Nyquist growth = R(z) to 7e-16; h_crit = 2.7853/(nu (N/2)^2) = 0.1088 (N=32, nu=0.1), growth there = 1 to ~1e-14. Blow-up runs match eps*R^n with the same roundoff seed eps ~ 7e-17; h = 0.111 (just above h_crit) only reaches 5e-3 by T=40, so "bounded at finite T" does not prove stability.
- test_timestep now asserts y' = 10y matches e^10; mean-conservation tolerance loosened to 1e-12.
NEXT: integrating-factor RK4 in timestep.py (h_crit ~ 1/(nu N^2) is the stiffness it removes); then 1D Burgers with/without 2/3 dealiasing.

## 2026-10-17
Goal: integrating-factor RK4 (`if_step` in timestep.py), verified before any Burgers work. Source: Trefethen, Spectral Methods in MATLAB, Ch. 10 (Program 27); Kassam & Trefethen 2005 is the ETDRK4 follow-up, not read yet.
- Derivation (interaction picture): for u_hat_t = L u_hat + N(u_hat), set U = e^{-Lt} u_hat. Then U_t = e^{-Lt} N(e^{Lt} U), so L drops out of the linear part and RK4 only sees N. Stages need e^{L dt/2} and e^{L dt}. Convention chosen: L on the right-hand side, E = exp(L dt/2) with entries <= 1 (no positive exponentials, so no overflow for diffusion). Physical meaning: E is the exact viscous damping over half a step, so dt is limited by the nonlinear timescale (advective CFL), not by 1/(nu k_max^2).
- Derivation errors caught on the way: missing dt factors, wrong sign in the exponentials, missing the factor 2 on the middle RK4 weights, reusing a stage value instead of the new one.
- Implementation checks (tests/test_timestep.py):
  - `NoNonlinear`: N = 0, nu k_max^2 dt ~ 1e3 (plain RK4 would blow up) -> matches exp(L dt) u_hat to 1e-14.
  - `NoLinear`: L = 0, N = w^2 -> agrees with `rk4_step` to 1e-14 over 100 steps (this caught a missing N() call on stage a4).
  - `test_order`: Bernoulli y' = -lam y + mu y^2, lam=5, mu=1, y0=1, T=1, exact y = 1/(mu/lam + (1 - mu/lam) e^{lam}) from z = 1/y (z' = lam z - mu, linear). This is the only test that exercises the E-powers on the nonlinear stages, since the other two have E = 1 or no stages. Error ratios per halving give observed order 3.83, 3.92, 3.96, 3.98 (n = 10..160), approaching 4 from below. The 5% tolerance on log-ratio passes at the coarsest pair with little margin (3.83 -> ratio 14.2); do not tighten it.
- Bugs hit: `np.max` on a complex array compares lexicographically, not by magnitude, so it silently ignored most of the error (use `np.abs` first); in the NoLinear loop `u` was overwritten before being compared with the `rk4_step` output (stale-state comparison).
- Not done yet: stiff check in the Bernoulli problem (lam=100, dt=0.05, lam*dt=5 beyond RK4's ~2.8 limit; assert bounded and accurate for `if_step`, unstable for `rk4_step`; no order assertion since order reduction at large lam*dt is expected). This is the miniature of the logbook figure.
- Committed: bf1980a (timestep: integrating-factor RK4 with linear, nonlinear and order tests).
NEXT: stiff Bernoulli check; `burgers.py` (L = -nu k^2 on the layout of `wavenumbers(n)`, dealias the product with the 2/3 mask, take .real before squaring, Cole-Hopf initial condition checked against the setup before running dynamics); dt-convergence (slope 4), then spatial convergence on log-linear axes. No 2D until that plot is exponential.

## 2026-10-19
Goal: close out timestep.py with a stiff check of `if_step` against plain `rk4_step`.
- `test_stability_and_blowup`: Bernoulli y' = -lam y + mu y^2 with lam = 70, mu = 1, dt = 0.05 (n_steps = 20), so lam*dt = 3.5, past RK4's real-axis limit of ~2.8 (|R(-3.5)| ~ 2.7 per step). Asserts `if_step` is within 1e-3 relative error of the exact solution and that `rk4_step` is more than 10x off.
- Measured: `if_step` relative error ~5e-4 at t = 3 dt. My first tolerance of 1e-10 failed. It assumed roundoff-level error, but this is truncation error of the nonlinear stages at lam*dt = 3.5, where the O(dt^4) constant grows with (lam dt)^p. The stepper was fine and the tolerance was wrong, so tolerances should come from a printed measurement, not a guess.
- Comparing against the exact value at T = 1 is meaningless here: y(1) ~ e^-70 ~ 4e-31, so any stepper that returns ~0 passes. Compare relative error at early times, where y is still O(1e-2..1).
- RK4's failure is partly the quadratic term, not only linear instability: once the amplified mode pushes y past lam/mu = 70, mu y^2 takes over. Written as `not (err <= 10)` so a nan from overflow still counts as failure.
- Limits of the test: the loop runs 4 iterations and the `i > 2` guard means both assertions fire once, at t = 3 dt. It shows accuracy and RK4 failure at one point, not boundedness to T = 1. The 1e-3 tolerance is only ~2x above the measured error; it is a regression guard, do not tighten it. n_steps = 20 fixes lam*dt = lam*T/n_steps, so changing it moves the stiff regime.
- 23 tests pass (`uv run pytest`).
NEXT: `burgers.py`. Periodic [0, 2 pi), u_t = -u u_x + nu u_xx, nu ~ 0.1, N ~ 64. Exact Cole-Hopf reference with theta(x, 0) = 1 + a cos x, 0 < a < 1 (start a ~ 0.5); derive u = -2 nu theta_x / theta at general t by hand and check periodicity and the PDE before coding. Check the t = 0 round trip through fft/ifft, plot |u_hat(k)| of the exact IC to confirm it reaches roundoff before the 2/3 cutoff. Then dt-convergence (slope 4), spatial convergence on log-linear axes, and the RK4-vs-IF stability figure on Burgers.

## 2026-10-20 .. 2026-10-22
Goal: 1D viscous Burgers (`solvers/burgers.py`), verified against an exact solution, as the Phase 0 deliverable.

### Advective term
- Conservative form: u u_x = (u^2/2)_x, so N_hat(k) = -i k * FFT[u^2/2](k), 2/3-masked. The product is formed in physical space and the derivative applied afterwards in Fourier space, so `ik` never multiplies `u` itself. This fits the `if_step` contract (N takes a physical field, returns a Fourier array) with one forward FFT per stage; the non-conservative u*u_x would need three extra transforms to get u_x from u alone.
- Mask and derivative are both diagonal in k and commute; the mask goes on the transform of the product. `.real` is taken before squaring so a drifting imaginary part cannot feed back.
- k[0] = 0 exactly, so N[0] = 0 for any field: N is a flux divergence, hence the mean of u is conserved to roundoff. The non-conservative form only conserves it approximately, since aliasing breaks the discrete equivalence.

### Exact reference (Cole-Hopf)
- u = -2 nu theta_x / theta solves Burgers exactly when theta solves the heat equation, so theta(x, 0) = 1 + a cos x evolves as theta = 1 + a cos x e^{-nu t} (mode k = 1, so the decay is e^{-nu k^2 t} = e^{-nu t}), giving u = 2 nu a sin x e^{-nu t} / (1 + a cos x e^{-nu t}).
- 0 < a < 1 keeps theta > 0 (the heat equation only damps), so u never blows up; u is the derivative of a periodic log, so its mean is exactly zero. a is the steepness knob: as a -> 1 theta nearly vanishes at x = pi and a front forms there, with width set by nu. Small-a limit: u ~ 2 nu a sin x e^{-nu t}, satisfying u_t = nu u_xx with the nonlinear term O(a^2).
- theta_x is differentiated BY HAND, not spectrally. If the reference used `diffspec`, an error in the derivative code would appear on both sides of the comparison and cancel. Same principle as z = 1/y in the Bernoulli test. Only u_exact(x, 0) and u_exact(x, T) are needed, not all t.
- The solver stays general (takes an `initial` array); only the reference is specific to theta_0 = 1 + a cos x.

### Diagnostics before any dynamics (tests/test_burgers.py)
- `test_u_exact`: finite-difference residual u_t + u u_x - nu u_xx at one point (h = 1e-4 first derivatives, 1e-3 second). Measured 7.23e-10, tolerance 1e-8. Finite differences, not analytic derivatives, so the test does not reuse the derivation it is checking.
- `test_N`: u = sin x gives -(u^2/2)_x = -1/2 sin 2x. Checks sign, the factor 1/2, the wavenumber layout, and that mode 1 feeds mode 2. Imaginary part also asserted zero.
- `test_N_mean`: field with a nonzero mean, N[0] == 0 exactly (not approximately).
- `test_exact_spectrum`: |u_hat[0]| at roundoff and |u_hat[n//3]|/max = 3.6e-12. The spectrum decays like k r^k with r = (1-sqrt(1-a^2))/a ~ 0.27 for a = 0.5, so at n = 64 the tail is ~1e-11 (NOT roundoff) and n = 96 is needed for the dt study to have no spatial floor.
- `test_aliasing`: u = cos(21x) (the top kept mode at n = 64). The product has content at k = 42, which aliases to -22, outside the mask. Correct framing: the 2/3 rule works on band-limited input, so the test feeds a mode INSIDE the mask and checks its alias is removed, rather than feeding a mode above n/3.
- `test_dynamics`: n = 96, nu = 0.1, a = 0.5, T = 1. Error 2.6e-15 at 200 steps, 1.5e-15 at 400: at the roundoff floor and not time-limited.

### dt-convergence (scripts/burgers_dt_convergence.py, figures/burgers_dt_convergence.png)
- n = 96, nu = 1, a = 0.5, T = 1, n_steps = 4..64. Errors 2.98e-5 .. 9.64e-10; observed orders 3.44, 3.71, 3.85, 3.92; fitted slope on the last three points 3.886. Approaches 4 from below, as in the Bernoulli test, and the finest point is far from roundoff so no plateau contaminates the fit.
- nu = 0.1 is useless as an order test: the amplitude is 2 nu a, so the problem is nearly linear, the integrating factor handles the linear part exactly, and every step count lands on the roundoff floor. nu = 1 raises the amplitude to 1 and is the first real test of the nonlinear stages. Larger errors at nu = 1 are therefore expected, not a regression.
- The spectrum depends only on `a`, not nu, so n = 96 keeps the spatial error negligible in both cases.

### Tolerance mistakes worth not repeating
- Every threshold here was set from a printed measurement, after three failures from guessed ones: 1e-10 on the stiff Bernoulli check (truncation error, not roundoff); 1e-17 then 1e-14 on the exact spectrum; 1e-14 on the aliasing noise (cos(21x) carries ~1e-14 of phase noise at argument ~132, the unnormalised FFT sums 64 of them, and the derivative scales by |k| up to 21, giving ~1e-13).
- Two test-shape bugs caught by inspection: `np.all(u_mask[22:43]) == 0` collapses the slice to one bool before comparing, so it only asserted that at least one slot was zero (it did still fail with the mask removed, but a mask zeroing one slot of 21 would have passed); and the second dynamics run called `solve_burgers` with `n_steps` instead of `2*n_steps`, so the step-doubling check tested nothing.
- Committed: f824e72 (burgers + tests), 05971cc (convergence script + figure).
NEXT: spatial convergence. Fix dt small enough that time error is below 1e-12 (nu = 0.1 with ~200 steps), vary n = 8, 12, 16, .., and plot max error vs n on log-LINEAR axes. Expect a straight falling line of slope ln r ~ -1.3 per mode for a = 0.5, flattening at ~1e-14 around n ~ 60. Then the IF-vs-RK4 stability figure on Burgers, and only then Phase 1 (2D grid, Poisson inversion with k = 0 explicit, Jacobian on an analytic pair, energy/enstrophy conservation at nu = 0).

## 2026-10-22 (spatial convergence)
Goal: spatial convergence of the Burgers solver, checked against a predicted decay rate rather than just "it goes down".
- Setup (scripts/burgers_spatial_convergence.py, figures/burgers_spatial_convergence.png): nu = 0.1, a = 0.5, T = 1, n_steps = 200 fixed so time error stays ~1e-15; n = 8, 16, .., 96. Plotted on log-LINEAR axes: spectral error is exponential in n, not a power of n, so there is no "order", only a decay rate.
- Prediction: if u is analytic in the strip |Im x| < y0, its Fourier coefficients decay like exp(-y0 |k|), and the truncation error at cutoff k_max scales like exp(-y0 k_max). u = -2 nu (ln theta)_x is singular where theta vanishes; with theta = 1 + a' cos x and x = pi + i y, cos(pi + i y) = -cosh y, so the zero is at cosh(y0) = 1/a' with a' = a e^{-nu t}. y0 is also the width of the front in x: a -> 1 sharpens it, y0 -> 0 and the required n diverges.
- Result: fitted decay rate 0.478 per unit n (natural log, fit on n >= 16 with errors above 1e-12). Prediction y0(T)/3 = 0.477. Ruled out: y0(0)/3 = 0.439 (8% off) and the full band y0/2 ~ 0.66-0.72.
- So the 2/3 mask sets the effective cutoff (k_max = n/3, not n/2): the solver resolves exactly the modes the dealiasing rule says it should.
- Correction to my own advice: I said to use y0 at t = 0, as the least smooth moment. Wrong. The error is measured at T, so it is the truncation tail of u(T). Two reasons: viscosity smooths the solution (a' shrinks, the singularity retreats from the real axis, y0 grows from 1.32 to 1.43); and truncation error made early sits at k ~ n/3 ~ 20, where the integrating factor damps it by e^{-nu k^2 t} = e^{-40 t}, so only recent error survives. The error is slaved to the current solution, not the initial condition.
- The third-digit agreement (0.2%) is partly luck: the true tail is ~ k r^k with a polynomial prefactor and the fit window shifts the number. The honest claim is agreement within ~5%, which is the test tolerance. It still rejects the t = 0 value and the full-band alternative.
- Floor at ~2.6e-15 for n >= 72, the same value as test_dynamics: roundoff, no information about the rate. Error at n = 64 is 1.3e-14.
- `test_spatial` (tests/test_burgers.py): asserts error at n = 64 < 1e-13, at least 4 points in the fit window (so polyfit cannot silently fit two), and |y0(T)/3 - rate|/rate < 5e-2. The n = 64 check sits after the loop so it fails loudly if 64 drops out of the grid.
- Committed: 4373725.
NEXT: stability figure (scripts/burgers_stability.py). Same Burgers problem, nu = 0.1, n = 64, so nu k_max^2 = 102.4 with k_max = n/2 (L acts on all modes, only N is masked). Choose dt with nu k_max^2 dt ~ 5 (about 20 steps to T = 1), past RK4's real-axis limit of 2.785. Plain RK4 on the full RHS L u_hat + N via rk4_step should blow up from a roundoff seed at |R(-5)| ~ 13.7 per step; IF-RK4 should stay accurate. Plot max|u_hat| vs t with the predicted growth line overlaid. Then Phase 1.

## 2026-10-22 (stability figure)
Goal: the figure that justifies the integrating factor — plain RK4 vs IF-RK4 on the same Burgers problem at the same step. Phase 0 deliverable complete.
- Setup (scripts/burgers_stability.py): nu = 0.1, a = 0.5, n = 64, T = 1, n_steps = 20 so dt = 0.05. k_max = 32 (L acts on all modes; only N is masked), so nu k_max^2 dt = 5.12, past RK4's real-axis stability limit of 2.785. Plain RK4 steps the full RHS L u_hat + N(ifft(u_hat)) through `rk4_step`; both runs use the SAME masked `N`, so the only difference between them is the treatment of L.
- figures/burgers_stability.png: amplitude and max-error vs t. IF-RK4 error flat at ~3e-11 for the whole run; plain RK4 climbs from roundoff and crosses the solution's size near t = 0.75, reaching 2.8e3 by t = 0.9.
- figures/burgers_stability_modes.png is the stronger result. Tracking |u_hat_k| per step for k = 32, 31, 30 gives measured growth 15.251, 11.460, 8.523 against predicted |R(-nu k^2 dt)| = 15.251, 11.460, 8.523 — three decimals, over 14 decades, with no bend. So the blow-up is EXACTLY the linear amplification of the grid-scale diffusion modes, not the nonlinearity and not aliasing.
- Why those lines never bend: the 2/3 mask kills N for |k| > 21, so modes 30-32 receive no nonlinear forcing and evolve under RK4 on the linear term alone at every amplitude. The solution blows up only because their huge u^2 content feeds back into the unmasked modes through `u.real`.
- The max-norm error fit (12.9-13.7 depending on the window) is NOT the number to quote: it mixes modes with different growth rates plus a pre-asymptotic plateau at RK4's own truncation error (~1e-10, same as IF-RK4's), so the fastest mode only dominates late. Per-mode tracking is the right diagnostic; the max-norm curve is the visual story.
- Mistake worth not repeating: my first blow-up guard tested max|u_hat| BEFORE the step and skipped stepping above 1e10. That cannot work here — the last step jumped 5.6e4 -> 3.8e35 in one step, because once u is large RK4 iterates the quadratic through its stages (k1 ~ u^2, k2 ~ u^4, k3 ~ u^8, k4 ~ u^16), so growth stops being ~15x per step and becomes a power law. A guard must test the RESULT of a step and discard it, not the input. The class of bug: a threshold on the input of an operation that amplifies super-exponentially.
- Two latent bugs found by inspection: the IF recording sat below a `continue`, so after RK4 blew up `errs_if` kept its np.zeros value and would have dropped to the log floor; and `ts = dt*np.arange(n_steps)` stopped at 0.95, so the run never reached T.
- The point, in scaling terms: plain RK4 needs dt < 2.785/(nu k_max^2) = 0.027 here (>= 37 steps), and that bound scales like 1/n^2 — at 256^2 in 2D the diffusion term alone would dictate the step. IF-RK4's limit is advective, dt ~ 1/(|u|_max k_max) ~ 1/n. That gap is why the integrating factor is in the solver.
- Committed: 73eb38b.
NEXT: Phase 0 is done (spectral derivative, RK4, heat, IF-RK4, Burgers with dt- and spatial-convergence and the stability figure). Phase 1: 2D grid and wavenumbers; Poisson inversion psi_hat = omega_hat/k^2 with the k = 0 mode set explicitly to zero; the Jacobian J(psi, omega) tested against an analytic pair; then nu = 0 with dealiasing must conserve energy and enstrophy to roundoff. Taylor-Green, then decaying turbulence with the k^-3 enstrophy cascade. Optional small test: mirror test_stability_and_blowup on the PDE (same dt, assert IF accurate and plain RK4 not, nan-safe).

## 2026-11-05 .. 2026-11-07
Goal: Phase 1, step 1. The 2D spatial operators (grid, d/dx, d/dy, Laplacian, Poisson inversion, 2/3 mask), each verified alone before any time stepping or the Jacobian. A time stepper hides RHS bugs.

### Grid (src/dmhd/spectral.py)
- Wavenumbers are k = 2 pi m / L. Periodicity fixes the spacing 2 pi / L. The grid fixes which N integers m can be told apart: at x_j = jL/N the phase is 2 pi m j / N, L cancels, so m and m + N agree at every grid point. `fftfreq(N, d=L/N)` returns m/L, hence the 2 pi. kx, ky are labels of the fft2 output slots, set by geometry only, never by the data.
- `meshgrid(..., indexing="ij")` puts x on axis 0, matching fft2. meshgrid returns a TUPLE of arrays; the first version stored the whole tuple.
- Odd derivatives use the Nyquist-zeroed array, the Laplacian and Poisson the full one (k^2 is even and real, so the Nyquist mode stays consistent).
- `Grid` is a frozen dataclass (`eq=False`) built by the classmethod `Grid.make(L, N)`; every array has `writeable = False`. Frozen blocks reassigning fields; the flag blocks in-place writes into the arrays, which frozen alone does not.
- Why it ended up frozen: the first version kept the wavenumbers as module-level globals, and `poisson_inversion` did `kx = kx_even_der; kx[0][0] = 1`. That is an alias, not a copy, so the first call silently changed k = 0 for every later `laplacian` call. Read-only arrays turn that class of bug into an exception.
- Other bugs on the way: two functions both named `partial_x` (the second, the y-derivative, silently replaced the first); `.copy` without `()` stored the method; the odd-derivative block assigned to the `*_even_der` names, overwriting them; field declared as `wavesood`, passed as `wavesodd`.
- Poisson inversion solves -lap psi = u (psi_hat = u_hat / k^2), matching omega = -lap psi. k = 0 is set to zero on a local copy, so the output has zero mean (the mean of psi is arbitrary: only grad psi enters the velocity).
- The operators return the complex ifft2 output, no `.real` inside, so tests can see the imaginary part. Open: take `.real` once, checked, on the physical-space nonlinear term.

### Tests (tests/test_spectral.py, 40 passing)
- `test_2d_derivatives`: f = sin(3x) cos(y), analytic d/dx, d/dy, lap to 1e-11 for L in {2pi, 4pi, 8pi}, N in {28, 30, 32}. Not symmetric in x and y, so an axis swap would be an O(1) error.
  - The test function must be periodic on the box. L = 5 or L = 1 cuts sin(3x) mid-period; the periodic extension jumps, coefficients decay like 1/k, and the spectral derivative has O(1) errors that do not shrink with N. Nothing to do with the code.
  - N = 20, L = 8 pi failed with error 5.000: sin(3x) is mode m = 3L/(2 pi) = 12 > N/2 = 10, which aliases to 12 - 20 = -8, i.e. wavenumber 2 instead of 3, and the max difference is 3 + 2. Under-resolution, not a bug. Resolution is a statement about m: need |m| < N/2, i.e. N > 2L/lambda. N alone means nothing until L is attached.
- `test_poisson_mean`: f = c + sin(2 pi 2x/L) cos(2 pi 3y/L) + 0.5 cos(2 pi x/L), c = 1.7. Asserts -lap(P f) = f - c, -P(lap f) = f - c, mean(P f) = 0. Nonzero mean on purpose (a zero-mean field would hide the k = 0 handling); different modes in x and y so nothing cancels by symmetry. First version failed with a CONSTANT error ~2.7 = 2 max|f - c| at every L and N: the test expected the opposite sign. Error/expected = 2 names the bug class (sign) before reading any code. The code was right; the test was wrong.
- `test_imaginary`: random real field (every mode, including Nyquist); imaginary part of all four operators < 1e-12. Smooth test fields never touch Nyquist, so they cannot check the zeroing.
  - Control: pointing `partial_x` at the unzeroed array gave O(1) imaginary parts, 5.1, 2.6, 1.3 at L = 2pi, 4pi, 8pi. Halves with L because k_Nyq = (N/2)(2 pi / L).
  - Why: a real field has f_hat(-m) = conj f_hat(m). The Nyquist slot pairs with itself, so f_hat there is real, and i k_Nyq f_hat is purely imaginary unless k_Nyq = 0. On the grid: the derivative of (-1)^j is -sin(pi j) times a constant, zero at every grid point.
- `test_mask_count`: (2M + 1)^2 kept modes, N in {14, 28, 30, 32}. N = 30 is the case where N/3 is an integer.
- `test_mask_symmetry`: mask[ix_(-i mod N, -j mod N)] == mask.

### 2/3 mask
- Square mask |mx| <= M and |my| <= M on INTEGER mode numbers (not k, so it does not depend on L). The 2D DFT is separable and each axis wraps mod N independently, so the per-axis condition is the right one.
- M = (N - 1)//3, the largest integer strictly below N/3. Two kept modes with |m| <= M multiply into |m| <= 2M, which wraps to 2M - N; that must lie outside the kept band, so N - 2M > M. With "<= N/3" at N = 30: 10 + 10 = 20 wraps to -10, which is kept, so the top mode gets aliasing error. Nyquist is excluded automatically (N/2 > M).
- The square mask keeps corners up to |k| = sqrt(2) N/3, so it is anisotropic. Fine for correctness; E(k) shell sums must stop at N/3.
- Symmetry: a real 0/1 mask with mu(-m) = mu(m) maps conjugate-symmetric spectra to conjugate-symmetric spectra, so real fields stay real (paired terms add to 2 Re(...)). An asymmetric mask leaves a lone e^{i theta} term, which is complex. The index test checks mu(-m) = mu(m) directly, which implies realness, so no separate operational test on the array. That becomes worth having once the Jacobian applies the mask, because then it tests the application code.
- Bugs on the way: `fftfreq(N, d=L/N)` for mode numbers gave m/L, not m (units; `d=1/N` gives integers); Python `and` on arrays raises "ambiguous truth value" (needs `np.logical_and` or `&` with parentheses).
- Test-shape bugs: `(np.arange(N) - N) % N` is the identity, not the mirror; two 1D index arrays pick a diagonal, not a grid (needs `np.ix_`); `sym.all() == mask.all()` is False == False and passes for ANY mask. Same class as `np.all(slice) == 0` in the Burgers entry: reducing to one bool before comparing. A test that has only ever passed proves nothing; break the code on purpose once.
- Committed: 1941caf (branch phase1-2d-operators).
NEXT: (1) Restore the immutability test: `pytest.raises(FrozenInstanceError)` on `g.N = 5`, `pytest.raises(ValueError)` on an array write, in SEPARATE `with` blocks, plus a loop asserting `flags.writeable` is False on every array field including `mask_twothirds`. (2) Jacobian J(psi, omega) = psi_x omega_y - psi_y omega_x: derivatives spectrally, product in physical space, fft, mask the product. Decide the interface (physical or spectral in/out) and where `.real` is taken. Test A: hand-derived J for a NON-proportional pair with mixed x, y dependence and |m| well below N/3 (derive by hand first). Test B: <psi J> = <omega J> = 0 to ~1e-14 for random fields band-limited to M (energy and enstrophy conservation of the nonlinear term). Controls: mask off with content up to N/2 must fail; cutoff "<= N/3" at N = 30 must fail. Not Taylor-Green yet: omega is proportional to psi there, so J vanishes identically.

## 2026-11-08 (stepper transform, 2D)
Goal: make `if_step` safe for 2D before any 2D PDE uses it.
- Problem: `if_step` called `np.fft.ifft`, which transforms the LAST AXIS ONLY. On an N x N array it returns the right shape and a plausible wrong field, so nothing raises. Linear-only tests cannot see it (L acts in Fourier space, N is never called with real content), and a Jacobian test calls N directly, bypassing the stepper. Class of bug: a dimension-specific transform hidden inside a generic routine.
- Fix: `np.fft.ifftn` in all four stages. `ifft2` would have broken the 1D Burgers solver, whose tests are now the regression suite for the stepper; they pass unchanged (7 passed).
- Caveat: `ifftn` transforms ALL axes. Fine for one field. A stacked state (omega, psi_B) of shape (2, N, N) would be transformed across the field axis too. Decide before Phase 2: transform only the spatial axes, or pass the transform in.
- Test (`test_recordingN`, tests/test_timestep.py): N is a spy that appends each argument it receives and returns zeros, so with N = 0 every stage input is known in closed form. Asserts: exactly 4 calls; call 1 = ifft2(u_hat), calls 2 and 3 = ifft2(E u_hat), call 4 = ifft2(E^2 u_hat); and the returned step = E^2 u_hat. Setup: N = 24, L = 2 pi, random N x N field (so a one-axis transform gives an O(1) error), nonuniform L = -(kx^2 + ky^2) from the Grid. The reference uses `ifft2`, not `ifftn`, so a shared bug cannot cancel. All tolerances 1e-12.
- Control: with `ifft` restored in `if_step`, the test fails at call 1 with max error 16.164224170535206 (O(1), not roundoff); the call-count assert passes, so the failure is the transform and not the number of calls. With `ifftn` all 6 tests in test_timestep.py pass. The control was run on the final version of the test (dt = 1e-5). The error is identical to the one at dt = 0.1 because call 1 sees ifft(u_hat) with no E in it, so it does not depend on dt.
- dt: first 0.1, then 1e-5. At 0.1 the top modes have E ~ e^-14, so only low modes distinguish E from E^2. At 1e-5 all modes have E ~ 1, and the E vs E^2 difference is still ~1e-5 or larger, far above the 1e-12 tolerance.
- Also: `Grid.make` now asserts even N (the Nyquist slot N//2 is only the Nyquist mode for even N).
- Committed: (hash)

## 2026-11-08 .. 2026-11-09 (Jacobian)
Goal: the nonlinear term J(psi, omega), verified against a hand-derived result before it enters any solver.

### Convention (derived on paper, no test can check this part)
- v = (psi_y, -psi_x), omega = curl v = d_x v_y - d_y v_x = -lap(psi). With J = psi_x omega_y - psi_y omega_x one gets J = -(v . grad) omega, so the vorticity equation is omega_t = J(psi, omega) + nu lap(omega).
- Physical check: omega = f(psi) gives J = f'(psi) (psi_x psi_y - psi_y psi_x) = 0; streamlines and vorticity contours coincide, so nothing is advected.
- A global sign error in the CONVENTION is invisible to every test here: antisymmetry and conservation hold for either sign. Test A pins J to its DEFINITION only.

### Interface
- `jacobian(psi, omega, grid)` in spectral.py (the MHD bracket reuses it): physical real fields in, masked Fourier array out. Matches the `if_step` contract. Built from the tested `partial_x`/`partial_y`; `.real` taken on the four derivatives before multiplying. The redundant transforms (each derivative does its own fft2/ifft2) are accepted for now: correctness first, and JAX will fuse them.
- Bug on the way: the first version did `temp[grid.mask_twothirds]`, boolean indexing, which SELECTS the kept entries of the PHYSICAL-space product and returns a flat array (225 = (2M+1)^2 for N = 24, visible in the traceback). The mask is indexed by Fourier mode slot, so it must zero entries of the Fourier array, not select from physical space. Class of bug: a mask applied in the wrong space with the wrong operation. The tell was the shape in the error.

### Test A (tests/test_spectral.py::test_jacobian)
- Pair: psi = sin(2x) cos(y), omega = cos(x) sin(3y), N = 24 (M = 7), chosen non-proportional, both depending on x and y, inputs |m| <= 3, so J sits well inside the kept band and the test isolates the formula from aliasing and masking (Test B's job).
- Hand derivation (product-to-sum): J = 5/8 [cos(x-2y) + cos(3x+4y) + cos(3x-4y) + cos(x+2y)] + 7/8 [cos(3x-2y) + cos(3x+2y) + cos(x+4y) + cos(x-4y)]. Highest mode (3, 4); no cancellation between the two products; 16 distinct Fourier modes. J(0,0) = 6, but that only exercises the first product (psi_y = 0 there), so a second reference was added.
- Direct closed form: J = 6 cos2x cos x cos y cos3y - sin2x sin x sin y sin3y. Asserted equal to the mode sum on the whole grid in pure NumPy: this audits the reference itself, independent of any spectral code, so a failure of the main comparison can be blamed on code and not on my algebra.
- Assertions: max |ifft2(J_hat).real - J_ref| < 1e-12; max |ifft2(J_hat).imag| < 1e-13 (the mask application keeps the field real); reference vs direct formula < 1e-13; antisymmetry J(omega, psi) = -J(psi, omega) < 1e-13.
- Control: with the sign flipped (`-jacobian`), the first assertion fails with max error 12.000000000000016 = 2 max|J|, since max|J| = 6. Without the flip, 41 tests pass in test_spectral.py.
- Tolerances: the first thresholds were guesses; the measured error was 7.99e-15 against a threshold of 1e-15. 1e-14 would leave about 25% margin, which is fragile across NumPy and platforms, so the three tight ones were loosened to 1e-13. A real bug here is O(1), so 1e-13 loses nothing.
- Not done: the kx/ky swap control, and the exact-spectrum check (16 nonzero modes at N^2 c / 2). The mask is not exercised by this pair on purpose.
- Committed: (hash)
## 2026-11-09 (Jacobian conservation, Test B)
Goal: the nonlinear term conserves energy and enstrophy exactly on the grid, with controls showing the test can fail.

### What is tested
- E = 1/2 <psi omega>, Z = 1/2 <omega^2>; dE/dt = <psi omega_t> (-lap is self-adjoint), so the nonlinear contribution is <psi J(psi, omega)>, and dZ/dt gives <omega J(psi, omega)>. These are exactly the integrals the test checks. The weights are the functional derivatives, dE/domega = psi, dZ/domega = omega.
- Why they vanish: int f J(g, h) is cyclic in (f, g, h) (write J(g, h) = d_x(g h_y) - d_y(g h_x), integrate by parts, use antisymmetry). If two slots hold the same field, the integral equals its own negative. This holds for ANY pair (psi, omega), so psi and omega are drawn independently and are not tied by omega = -lap(psi). Physically: v . grad(psi) = 0, advection does no work; enstrophy is a Casimir (int F(omega) conserved for every F).
- Why it holds exactly on the grid: mode level, <psi J> = sum over p + q + r = 0 of (p x q) psi_hat(p) omega_hat(q) psi_hat(r); the terms cancel in pairs, but only if p + q + r = 0 exactly. The grid cannot tell k from k +- N, so triples with p + q + r = +-N would be wrongly counted as zero-sum. With |m| <= M = (N - 1)//3 the largest sum is 3M < N, so none exist. That is the 2/3 rule, and it only protects band-limited input.

### Test (tests/test_conservation.py::test_jac_conservation)
- N in {24, 30, 32}, 3 psi seeds x 3 omega seeds = 27 runs. Random real fields, band-limited by zeroing the spectrum outside `g.mask_twothirds`; imaginary part of the result asserted at roundoff before `.real`.
- Asserted: |<psi J>| / (rms psi rms J) < 1e-14 and |<omega J>| / (rms omega rms J) < 1e-14 (measured maximum 8.9e-17, below machine epsilon because the mean of N^2 random-sign terms cancels); and rms(J) / (rms|grad psi| rms|grad omega|) > 0.1 as a guard against a Jacobian that returns zeros (so conservation would hold trivially). By Cauchy-Schwarz the ratio lies in [0, 1].
- Measured guard ratio 0.37 to 0.43. For independent isotropic gradients the unmasked value would be 1/sqrt(2) = 0.71; the mask removes the modes of J above M (J has modes up to 2M), which is why it is lower.

### Controls
- Wrong bracket (psi_x omega_y + psi_y omega_x): all 27 runs fail, relative error 1.2e-3 to 9.1e-2 against a 1e-16 threshold.
- Wrong cutoff (M = N//3 for both the input band-limiting and the Jacobian mask, via `dataclasses.replace` on the frozen Grid): fails at N = 24 and N = 30 (relative error 4.4e-3 to 3.0e-2), passes at N = 32, where N//3 = (N - 1)//3 so the two masks are identical. Kept permanently as `test_jac_failure` (N in {24, 30}, asserts the error is above 1e-4).
- Inputs not band-limited (content to N/2): not run separately. The failure mechanism is the same wrapped triples, and the wrong-cutoff control is the harder case (only a sliver of triples, those with all three modes at exactly +-N/3, can wrap), so bulk aliasing fails by more.
- Flipping the sign of J is NOT a control here: -J conserves energy and enstrophy equally.

### Mistakes worth not repeating
- First guard ratio used mean(|grad|^2) without the square root. Symptom: values 0.01 to 0.03 drifting with N (the mean square grows like N^2), which a scale-free ratio cannot do. Corrected value 0.37.
- `.real` inside partial_x, partial_y, laplacian and poisson_inversion makes `test_imaginary` vacuous (the imaginary part of a real array is exactly 0) and hides the Nyquist-zeroing control (O(1) imaginary parts). Reverted; `.real` is taken at the call sites. Class of bug: a projection moved into a primitive hides the signal a test relied on.
- `dataclasses.replace(g, ...)` returns a NEW frozen grid; the first attempt discarded the return value, so the "bad mask" control changed nothing and 27 tests still passed. A test that cannot fail is not a control. Also, control code was left inside `jacobian` in the source and silently overrode the mask for every caller; removed.
- Tight thresholds set 10 to 25 percent above the measurement are fragile; set them one to two orders above the measured maximum instead.
- Committed: (hash)
## 2026-11-10 (mask inside `jacobian`)
Goal: show that `jacobian` zeroes the Fourier array outside the 2/3 band, with a test that fails both when the mask is missing and when it is too narrow.

### Why a dedicated test
The product of two fields band-limited to M has modes out to 2M. With 3M < N the wrapped modes land outside the band, so the in-band modes are exact and the mask discards the rest. Test A (J modes <= 4 < M) and Test B (inputs band-limited, so <psi J> only reads J in band) cannot see a missing mask. If it were missing, the symptom would appear only in long runs, as spectral pile-up near the cutoff and energy drift from aliased energy folding back.

### Test (tests/test_spectral.py::test_jac_mask)
- N in {24, 28, 32} x 3 psi seeds x 3 omega seeds = 27 runs. N = 28 is the tight case (3M = N - 1).
- Region |m| <= M = (N-1)//3 built from `fftfreq` integer mode numbers, independent of `g.mask_twothirds`. Inputs band-limited by zeroing outside the Grid mask.
- Unmasked reference computed inside the test from partial_x, partial_y and fft2, without calling `jacobian`.
- Asserts: (1) max|J_hat| outside the region == 0 exactly. (2) max|unmasked outside| / max|unmasked inside| > 0.1, so that (1) can fail. (3) max|unmasked * mask - J_hat| < 1e-12.

### Measured
- Ratio in (2): 1.06 to 2.05 over the runs visible in the terminal, no trend with N, so the 0.1 threshold has about a 10x margin. Physically: the product's spectrum is a convolution out to 2M, so the content outside the band is as large as the content inside; the mask throws away a lot.
- First version of (2) was an absolute threshold (> 10) in unnormalised fft2 units; measured values were ~600 to ~3000 and grew like N^2. Replaced by the ratio, which is scale-free.

### Controls
- Mask multiplication removed from `jacobian`: (1) fails 27/27, outside magnitude ~600 to ~3000.
- Grid mask narrowed by one mode (M - 1): (3) fails 27/27, magnitude ~270 to ~1900 (N = 24: ~270-380, N = 28: ~790-1060, N = 32: ~1140-1910), while (1) still passes. So (1) guards a missing or too-wide mask, (3) a too-narrow one.
- Both control edits reverted; `git diff src/` clean before the commit.

### Notes / mistakes worth not repeating
- (3) is nearly tautological: both sides come from the same operations, so the difference is exactly 0 and the 1e-12 carries no margin. Its only value is catching in-band truncation, which is why it needed its own control.
- Absolute thresholds on unnormalised FFT quantities drift with N; use ratios.
- Control edits live in the source file while running; revert and check `git diff src/` every time (second time this has needed saying).
- The test is about `jacobian`'s contract, not conservation, so it lives in test_spectral.py.
Committed: (hash)
## 2026-11-13 .. 2026-11-14
Goal: first validation tests for `nse2d.py`: the linear part, then the nonlinear wiring and sign.

### Conventions fixed
- ω = ∂ₓv_y − ∂_yvₓ, ω = −∇²ψ, v = (ψ_y, −ψ_x), J = ψₓω_y − ψ_yωₓ.
- v·∇ω = −J, so ∂ₜω = +J + ν∇²ω and `nonlinear` returns +J.
- State, input and output are vorticity. Zero mean flow assumed (k = 0 mode of v is conserved and removable by a Galilean shift on a periodic box). Not true for a mean B field in MHD, which changes the physics; ⟨B⟩ = 0 must be assumed there.
- Decision: no immutability test. `Grid` stays frozen with read-only arrays (incl. `mask_twothirds`); the test would only catch a forgotten field later.

### Tests added (tests/test_nse.py)
- `test_linear`: N ≡ 0 via monkeypatch of `dmhd.solvers.nse2d.nonlinear`, calls `solve_nse`, compares with ω₀e^{−13ν(2π/l)²T} for ω₀ = sin(2kx)cos(3ky). Parametrised over N, n_steps (1 to 1000) and l ∈ {2π, 3}. Guard: decay > 0.01.
- `test_nonlinear`: ω = cos(3x+y) + cos(2x+y) (wavevectors non-parallel, |k|² = 10 and 5). Hand-derived J = −(1/20)[cos x − cos(5x+2y)], amplitude independent of l. Checks the full Fourier array against the FFT of the hand field, that the four predicted modes (±1,0), (±5,±2) are nonzero, and that everything else is ~0. All relative to max|ĵ_hand|.

### Measured
- Linear control (drop k_y²): error 0.297 predicted, 0.2967 measured (l = 2π, ν = 0.05); 0.358 predicted, 0.3581 measured (l = 3, N = 24). Error identical for 1 to 1000 steps (integrating factor is exact, so no dt dependence).
- Sign-flip control in `nonlinear`: relative error 1.9999999999999993 to 2.0 for every N and l, which is exactly 2|J|.
- Skip Poisson inversion: J(ω, ω) = 0, max = 0.0.
- x/y swapped in `predicted`: |ĵ| at the swapped modes ~1e-15, against a threshold of N²/80.

### Controls
- Dropped k_y² in `L` (in `nse2d.py`, not in the test): fails 24/24.
- −J, skipped Poisson inversion, swapped mode indexing: each fails, all reverted, `git diff src/` checked.
- Mask removed: not a useful control here. All content (|m| ≤ 6) lies inside the mask and below Nyquist, so the output is unchanged. The mask is covered by the `jacobian` test in `test_spectral.py`.

### Notes / mistakes worth not repeating
- The first `test_linear` rebuilt its own L and time loop, so controls on it proved nothing about `nse2d.py`. Always call the function under test; monkeypatch what you need to switch off.
- The `decay` guard caught ν = 0.2 at l = 3 (decay 1e-5): the factor (2π/l)² makes a box with l < 2π much more diffusive.
- Mode prediction for J: self-interactions vanish (ω_a ∝ ψ_a), and cos × cos gives sum and difference modes. I first listed (6,2), (4,2) and missed (1,0).
- A mode pair is (m, n) and (−m, −n). Writing (±5, 2) is wrong.
- Use a bool array as a mask, not a float one.
- Wavenumber for mode m is 2πm/l, not l m/2π. They agree at l = 2π and hide the bug.
- Unnormalised `fft2` coefficients scale as N², so use ratios. Take the scale from the reference, not from the output under test.

Committed: (hash)
NEXT: (1) Taylor–Green: N(ω₀) at roundoff, then ω₀e^{−2ν(2π/l)²T} for l = 2π and 3, small and large dt. It does not test the N coupling, since N ≈ 0 on that field. (2) dt self-convergence on the `test_nonlinear` field: ratio → 16, fixed n, dt below advective CFL. (3) Energy budget dE/dt = −2νZ with the mask on.