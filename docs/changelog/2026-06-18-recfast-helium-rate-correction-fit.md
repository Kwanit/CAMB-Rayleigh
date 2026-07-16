# RECFAST helium rate correction fit

- Date: 2026-06-18
- Scope: optional RECFAST helium recombination correction, compared to direct CosmoRec vX histories.

## Motivation

The existing seven-parameter `recfast_cosmorec` approximation matches the overall CosmoRec history well, but the HeI recombination part still leaves a coherent residual around the He recombination feature. A first pass fitted a free multiplier to the RECFAST He RHS, but a log/exp shape with explicit tapers was unnecessarily expensive and visually arbitrary.

## Adopted trial form

The current trial uses a cheap algebraic multiplier on the He equation,

```text
F_He(z) = 1 + (a0 + a1*u + a2*u^2)/(1 + d2*u^2 + d4*u^4),
u = (z_scale - z0)/width,
z_scale = (T_CMB / COBE_CMBTemp) * (1 + z) - 1.
```

Here `F_He` is a direct multiplicative correction to the RECFAST HeI ionization-fraction derivative `dx_He/dz` after the singlet and triplet He terms have been summed. It does not alter the individual atomic rates or escape probabilities. Using `z_scale` makes the fitted feature depend on photon temperature rather than raw redshift, matching the convention used by the existing RECFAST hydrogen fudge.

The fitted constants are:

```text
a0    =  0.07805480599148856
a1    =  0.899368710079757
a2    = -3.063896868095767
d2    =  3.752520187028219
d4    = 18.126758791787946
z0    = 1968.1134791219417
width =  774.0312440582285
```

It is applied only when `RECFAST_He_rate_correction` is true, `Heflag == 6`, and `1500 < z_scale < 3000`. The constants are internal to `recfast.f90`; only the boolean switch is exposed through the Python/ini interface.

The fit target was total `x_e(z)` from `camb.get_background` against direct CosmoRec vX, not the pointwise inferred `dx_He/dz`. This matters because small high-redshift RHS errors accumulate into visible history shifts.

## Current comparison

Using the standard test cosmology in the diagnostic scripts:

- Baseline RECFAST vs CosmoRec, `z > 800`: RMS `4.009e-4`, max `1.391e-3`.
- Rational He correction, `z > 800`: RMS `5.386e-5`, max `1.738e-4`.
- Baseline RECFAST in the He window: RMS `6.424e-4`, max `1.391e-3`.
- Rational He correction in the He window: RMS `3.509e-5`, max `7.749e-5`.

The remaining `z > 800` maximum is at `z ~= 1328`, i.e. outside the He correction target region.

Refitting the existing seven RECFAST nuisance parameters on top of the new He-rate correction moved the hydrogen Gaussian parameters only modestly and did not further improve the He-window residual. The joint fit used `800 < z < 3000`, terminated after 6 evaluations, and found:

| Parameter | Start | Joint fit | Delta |
| --- | ---: | ---: | ---: |
| `RECFAST_fudge_He` | `0.8472384364` | `0.8472367976735188` | `-1.6387264811790203e-6` |
| `AGauss1` | `-0.1421246011` | `-0.13952724827828716` | `0.002597352821712834` |
| `AGauss2` | `0.0747533979` | `0.07298919518696294` | `-0.0017642027130370663` |
| `zGauss1` | `7.2841557503` | `7.281306128192896` | `-0.002849622107103933` |
| `zGauss2` | `6.7306248459` | `6.766703867907752` | `0.03607902200775115` |
| `wGauss1` | `0.1741263787` | `0.16389664101046253` | `-0.010229737689537471` |
| `wGauss2` | `0.3349828575` | `0.2785834127132674` | `-0.05639944478673259` |

Joint-refit `x_e` residuals:

- `z > 800`: RMS `4.211e-5`, max `1.141e-4`.
- He window: RMS `3.509e-5`, max `7.749e-5`.

CMB lensed `lmax=3000` residuals:

- TT RMS: baseline `5.143e-5`, corrected `5.141e-5`.
- EE RMS: baseline `7.126e-5`, corrected `6.655e-5`.
- TE scaled RMS: baseline `3.058e-5`, corrected `2.832e-5`.

Timing was checked with the dedicated Fortran `InitVars` timing executable:

```text
make -C fortran initvars_timing
fortran/initvars_timing 6000 5 baseline
fortran/initvars_timing 6000 5 he_rate
```

Two orderings were run, each with 5 warmup calls discarded and 6000 timed calls:

| Mode | Run 1 avg CPU seconds | Run 2 avg CPU seconds | Mean |
| --- | ---: | ---: | ---: |
| baseline | `2.483264e-2` | `2.565465e-2` | `2.524365e-2` |
| He-rate correction | `2.543300e-2` | `2.500467e-2` | `2.521884e-2` |

The ordering scatter is larger than the mean difference, so this timing does not show a measurable `InitVars` overhead. The correction itself is only multiplies/adds/divides and one branch in the He-active region.

## Outputs

Updated plots and summaries:

- `docs/changelog/recfast_cosmorec_head_fits/xe_recfast_he_rate_correction_vs_cosmorec.png`
- `docs/changelog/recfast_cosmorec_head_fits/cmb_recfast_he_rate_correction_vs_cosmorec.png`
- `docs/changelog/recfast_cosmorec_head_fits/helium_rhs_fudge_multiplier_fit.png`
- `docs/changelog/recfast_cosmorec_head_fits/helium_rhs_cosmorec_vs_recfast_fudge_fit.png`
- `docs/changelog/recfast_cosmorec_head_fits/helium_xhe_evolution_fudge_fit.png`
- `docs/changelog/recfast_cosmorec_head_fits/recfast_he_rate_correction_trial_summary.json`
- `docs/changelog/recfast_cosmorec_head_fits/helium_rhs_fudge_fit_summary.json`

## Validation

- Rebuilt CAMB with `RECOMBINATION_FILES="recfast hyrec cosmorec"` and `COSMOREC_PATH` pointing to `external/CosmoRec.vX`.
- Regenerated the diagnostic plots and JSON summaries after the z-gated internal-constant implementation.
- Built and ran `fortran/initvars_timing` for baseline and He-rate correction modes.
