# CorrFuncFullSky optimization and d4m4 recursion fix

Date: 2026-06-12

## d4m4 Wigner-d recursion fix (correctness)

The `d4m4` recursion used by the full curved-sky lensing method (lensing_method=4)
and by `camb.correlations` (both `lensed_cls`/`lensed_correlations` and
`lensed_cl_derivatives`) used `(l+2)(l-1)` in the coefficients where the exact
identity requires `l(l+1)`:

```
d4m4 = (-7/5 (l(l+1)-6) d2m2 + 12/5 (-l(l+1) + (9x+26)/(1-x)) d3m3) / (l(l+1)-12)
```

Checked against exact sympy `Rotation.d(l,4,-4,beta)` with exact `d2m2`/`d3m3`
inputs: the `l(l+1)` form (already used by the long-standing apodized Fortran
method 1) is exact to ~1e-13; the `(l+2)(l-1)` form has relative error ~2/l^2
(factor ~1e3 at l=5, 2.5e-3 at l=30). The term only enters the second-order
C_gl,2^2 contribution to the Q-U correlation, so the physical effect is tiny:
lensed BB changes by at most ~2e-8 (relative) at default settings. Fixed
consistently in `fortran/lensing.f90` and `camb/correlations.py` (two places),
preserving the Python/Fortran parity test at rtol 1e-6.

## CorrFuncFullSky optimization (lensing_method=4, default for AccurateBB=T)

`CorrFuncFullSkyPointAccumulate` was inlined into `CorrFuncFullSky` (matching the
structure of the other correlation routines) and restructured:

- The two corr accumulation loops were merged, and the per-l
  `exp(-l(l+1)sigma^2/2)` (previously evaluated twice per l per point) replaced
  by a multiplicative recurrence: one `exp` per integration point, relative
  error growing only like l*epsilon (~1e-12 at l=8000).
- All x-independent square roots, reciprocals and the `theta_cut` stability
  threshold are precomputed once per call instead of per (l, point); the builds
  use plain `-O3` (no fast-math), so per-l divisions/sqrts were full cost.
- The blanket zero-initialization of 16 work arrays per point was removed
  (every element read is written first).
- The high-m d's (`d1m2`...`d4m4`) are computed as scalars inside the merged
  loop instead of stored arrays, halving per-thread stack and memory traffic;
  only `P, dP, d11, dm11, d20, d22, d2m2` are stored, as in the apodized method.

Timings for the relensing step alone (`get_lensed_cls_with_spectrum`,
lmax=2500, lens_potential_accuracy=1, 20-core machine):

| Path                         | Before | After |
|------------------------------|--------|-------|
| method 4 AccurateBB=T (1 thr)| 1.00 s | 0.33 s|
| method 4 AccurateBB=T        | 245 ms | 86 ms |
| method 4 AccurateBB=F        | 11.5 ms| 3.4 ms|
| method 1 (default)           | unchanged (bit-identical) |

Method 4 spectra change by <2e-10 (TT/EE/TE) and <2e-8 (BB) relative, from the
exp recurrence and the d4m4 fix.

## Shared setup and other tidying

- The duplicated setup in the two curved-sky routines (lmax_lensed
  determination/allocation, weighted spectra, high-L template tail,
  ALens_Fiducial override, amplitude sanity check) is factored into
  `InitLensedClArrays` and `PrepareLensedCLSpectra`. The low-L EE taper is now
  an explicit argument: applied for the apodized scheme (unchanged) and now
  also for the Gauss-Legendre method in its truncated (AccurateBB=F) mode;
  the full-range AccurateBB=T reference remains untapered. Reason: the
  untapered truncated integral is window-sensitive in the lowest-L BB for high
  tau (errors growing from -2.3e-2 at L=2 to -6.6e-2 as LensingBoost widens
  theta_max from pi/32 to pi/8 for tau=0.11, i.e. not converging to the
  full-range result), whereas the tapered version is stable across the same
  sweep at the ~-2e-2 plateau set by the deterministic omission of the
  reionization-EE term, matching the long-standing method-1 behavior.
  `camb.correlations.lensed_cls`/`lensed_correlations` gained a matching
  `low_l_ee_taper` option (kernel-only, smoothstep zero at l=2 to one at l=20)
  so the python/Fortran short-range parity test compares like with like.
- Physical note: the reionization-bump EE lensing term kept by the untapered
  full-range reference is computed with the z~1100 deflection kernel, but only
  ~55-75% of the deflection power at the relevant L<~60 is sourced below z~7
  (checked with a `GaussianSourceWindow(source_type="lensing")` at z=7). The
  full-range untapered reference therefore overestimates the term by ~25-45%
  of its size (term is +1.7e-2 relative in BB at L=2 for tau=0.11, +3.5e-3 for
  tau=0.054), while tapering removes ~all of it: untapered full range remains
  the more accurate reference by roughly a factor of 2, and getting low-L BB
  right at the few-1e-3 level for high tau would need a multi-source kernel
  treatment.
- Python optimization: the per-point evaluation of the Legendre function table
  P_l(x), P_l'(x) in `legendre_funcs` (distinct from the Gauss-Legendre
  nodes/weights, which were already cached and computed by the Fortran
  `gauss_legendre`) now also uses the Fortran library: new
  `MathUtils.Legendre_Table` exposed as `camb.mathutils.legendre_polynomials`,
  with the previous scipy `legendre_p_all` retained as fallback when the
  compiled library is unavailable. Agrees with scipy to ~1e-12 and benchmarks
  ~3x faster per point (36us vs 116us at lmax=4400); end-to-end
  `lensed_cls` is ~23% faster for the default short range (60->46 ms) and ~9%
  for the full range (1726->1565 ms), the remainder being the per-point numpy
  d-function/correlation algebra. A blocked vectorization of that remaining
  numpy work (BLAS matvec reductions over chunks of points) was also tried and
  validated to 1e-14 but benchmarked neutral-to-slower (the ~lmax-sized
  per-point ops are bandwidth-bound and cache-resident; blocking pushed the
  working set out of L2), so the simple per-point loop was kept.
- The misplaced high-L-template normalization check was moved out of the tail
  fill loop (it tested the same element every iteration).
- `GetCachedGaussLegendre` now documents that the module-level cache is not
  safe for concurrent use with different npoints.

`camb.tests.camb_test` passes (18 tests), including the Python/Fortran method-4
parity checks.
