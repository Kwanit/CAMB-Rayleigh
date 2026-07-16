# Massive-Neutrino Momentum Kernels

## Summary

CAMB evolves massive-neutrino perturbations at a small number of sampled comoving
momenta `q` and integrates them with precomputed weights. The current optimized
rules live in `fortran/massive_neutrinos.f90`, and the maintained reproduction and
diagnostic script is `scripts/nu_integration_kernels.py`.

The script replaces the old Mathematica notebook reference and is self-documenting;
run `python scripts/nu_integration_kernels.py --help` for the full case list, current
findings, and usage notes.

## Live Sampling Rules

`TNuPerturbations_init` sets `nqmax` from `Accuracy * neutrino_q_boost`:

- `Accuracy <= 1`: use the default 3-node optimized rule.
- `1 < Accuracy <= 2`: use the current 4-node optimized rule.
- `2 < Accuracy <= 3`: use the current 5-node optimized rule.
- `Accuracy > 3`: use `nint(Accuracy * 10)` uniformly spaced midpoint samples.

The current hard-coded low-node rules are:

| `nqmax` | Node choice | Weight choice | Status |
| --- | --- | --- | --- |
| `3` | Free nodes | Exact solution of the six Fermi-Dirac moment conditions `n = -4, -2, -1, 0, 1, 2` | Kept as the fast default; accurate enough for most default runs. |
| `4` | Free nodes | Least-squares fit to Fermi-Dirac moments `n = -4, -2..2` plus representative `v(am/q)` and `1/v(am/q)` kernels | Best tested 4-node rule. |
| `5` | Free nodes with moment constraints | Exact for moments `n = -4, -2..2`; remaining freedom fit to representative `v(am/q)` and `1/v(am/q)` kernels | Best tested 5-node rule. |
| `>5` | Uniform midpoint grid to `qmax = 12 + nqmax / 5` | Direct midpoint evaluation of the Fermi-Dirac kernel, normalized by `fermi_dirac_const` | Used only for higher accuracy/reference runs. |

The original fixed-node 4- and 5-node rules are kept as comments beside the live
rules for traceability. They are reproducible with `scripts/nu_integration_kernels.py`
but were not adopted after compiled-CAMB comparisons.

## Validation

The main comparisons used compiled CAMB spectra, not just static quadrature
errors. The standard validation setup was:

- `AccuracyBoost = 1.2`, `lSampleBoost = 3`, and `lmax = 3000`.
- Massive-neutrino models with `mnu = 0.06`, `0.12`, and `0.3`, using one massive
  eigenstate for the first two and three for `mnu = 0.3`.
- Matter power tests using `pk_kmax = 5`, with focused checks at `z = 0` and
  broader diagnostics at `z = 0`, `0.5`, `1`, and `2`.
- Reference runs at `nqmax = 80`, cross-checked against `nqmax = 160` and spot
  checked against `nqmax = 240`.

`nqmax = 80` is stable enough for ranking the low-node candidates, but it is not
an exact reference. Typical residual reference-level differences are around
`1e-5`, increasing to a few `1e-5` or `1e-4` in the higher-mass matter-power
stress cases.

The current 4-node and 5-node rules are the best compiled-CAMB candidates found.
The 5-node rule is not guaranteed to improve every max-error diagnostic over the
4-node rule, because the rules are not nested and the max matter-power error can
be localized around `k/h ~ 0.007-0.012`. In broader diagnostics the 5-node rule is
generally tied with or slightly better than the 4-node rule in RMS error, and it
remains the preferred high-accuracy low-node rule.

## Alternative Searches

Several alternative searches were tested:

- Static benchmark minimax over `{1, q^-2} x {1, v, 1/v}`.
- A closer static proxy including `v^2` and `q^-2 v^3` kernels used by the
  perturbation equations.
- Dense minimax searches over the proxy kernels.
- Post-hoc minimax fits to dumped evolved q-dependent integrands from compiled
  CAMB runs.

Some minimax candidates substantially improved the static proxy errors, and fits
to dumped evolved integrands improved their offline objective. These improvements
did not translate into better compiled spectra. The evolved perturbation
multipoles are themselves changed when the q nodes are changed, so post-hoc
integration of functions dumped from an `nqmax = 80` run is not a valid adoption
criterion by itself. The dumped-integrand fits were also intentionally diagnostic:
the unconstrained minimax solutions compressed nodes into the mid-q range and
lost the small-q and large-q moment matching that protects the full evolution.

Future changes to these kernels should therefore keep compiled-CAMB spectra as
the decisive objective, preferably with tight moment constraints retained. If
dumped integrands are used again as candidate generators, the fit should restrict
to the semi/non-relativistic range and avoid giving weight to epochs where the
relativistic expansion already makes the q dependence trivial.

A follow-up constrained dump test included the scalar `pinudot` integrand as well
as the `Nu_Integrate_L012` density/flux/pressure/shear integrands. Enforcing the
same analytic moment constraints as the 3-node rule gave physically sensible
nodes again, including high-q tail nodes. The best tested 5-node constrained
dump-fit candidate improved some matter-power max-error diagnostics, but it
worsened lensed CMB max errors and did not improve the broad matter-power RMS
diagnostics. It was therefore not adopted over the current live 5-node rule.

## Reproducing

Useful script entry points:

```bash
python scripts/nu_integration_kernels.py
python scripts/nu_integration_kernels.py --case camb-compare
python scripts/nu_integration_kernels.py --case camb-compare --reference-nqmax 160
python scripts/nu_integration_kernels.py --case dump-kernels
python scripts/nu_integration_kernels.py --case minimax-compare --minimax-global-search
```

Run `python scripts/nu_integration_kernels.py --help` for the current script behavior
and interpretation notes. Testing a new candidate in compiled CAMB still requires
temporarily editing `fortran/massive_neutrinos.f90`
and rebuilding before running the compiled comparison.
