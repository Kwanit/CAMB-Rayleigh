# Thermal Neutrino Background Fit Validation

`scripts/nu_background_quadrature.py` is the standalone validation script for
the thermal massive-neutrino background fits in
`fortran/massive_neutrinos.f90`. It numerically reproduces the coefficient
sets, branch cuts, mappings, and Horner evaluations in:

- `ThermalNuBackground_rho_P`;
- `ThermalNuBackground_rho`;
- `ThermalNuBackground_drho`.

The active density and pressure implementation uses low-`am` polynomials for
`am <= 0.42`, direct polynomials in a mapped coordinate for `0.42 < am < 70`,
and large-`am` asymptotic expansions thereafter. `drho` has its own low branch
below `am = 0.3`, a fit to `D = rho - 3P` up to `am = 70`, and the corresponding
large-`am` expansion. The script reproduces `rho` separately in its public API,
although it has the same numerical fit as the density returned by `rho_P`.

## Reference and reported quantities

The reference values are adaptive Fermi-Dirac integrals with `epsabs = epsrel =
1e-13`. The script evaluates `rho`, `P`, and `drho / adotoa`, where the latter
is the exact `a * d rho / da`. The default command evaluates 801 logarithmically
spaced points over `1e-4 <= am <= 1e4`, so it tests both the fit interiors and
their low- and high-`am` branches.

It also rebuilds the removed 1400-point, 100-step quadrature/spline table as a
historical baseline. This comparison is limited to its old table range,
`0.3 <= am <= 70`; it is not part of the live CAMB implementation.

With the current coefficients and default sample grids, the reported relative
errors are:

| Quantity | Live fit: max (RMS) | Removed spline: max (RMS) |
| --- | --- | --- |
| `rho` | `2.335e-5` (`5.901e-6`) | `9.377e-6` (`6.865e-6`) |
| `P` | `4.640e-5` (`1.826e-5`) | `3.996e-5` (`4.823e-6`) |
| `drho / adotoa` | `3.545e-4` (`4.605e-5`) | `6.878e-4` (`9.675e-5`) |

The live density and pressure maxima occur at `am = 0.60256` and `7.94328`.
The `drho` maximum occurs at `am = 0.29512`, immediately below its low-branch
cut at `am = 0.3`; the `D`-fit comment in the Fortran source applies to the
intermediate direct-fit branch rather than this boundary-adjacent low expansion.

## Usage

```bash
python scripts/nu_background_quadrature.py
python scripts/nu_background_quadrature.py --samples 1601
python scripts/nu_background_quadrature.py --am-min 1e-5 --am-max 1e5
```

The script makes no changes to the Fortran implementation. When altering the
live coefficients or branch points, update the mirrored constants and rerun
this script to record the new reference errors.
