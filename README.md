# differentiable-mhd

> **Status:** work in progress (Phase 0: 1D spectral warm-up)

A differentiable 2D resistive MHD solver (pseudo-spectral, periodic box, JAX) and a Fourier Neural Operator surrogate, with a systematic study of where the surrogate fails.

## Install

```bash
uv sync
uv run nbstripout --install
uv run pytest
```
