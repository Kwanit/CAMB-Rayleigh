# RECFAST CosmoRec/HyRec HEAD Fudge Refit

- Date: 2026-06-10
- Source HEAD: `a8c0c32ff590067577f4c2c7b38b756416a2aca7`
- Temp worktree: `/tmp/camb-recfast-fit-head-NofB30`
- Build: `RECOMBINATION_FILES="recfast cosmorec hyrec"` with local `external/CosmoRec` and `external/HYREC-2`

The fit used the current HEAD RECFAST implementation with default CAMB timesteps and default RECFAST numerical settings.  CosmoRec was patched only in the temp worktree to expose shell/absorption overrides through the CAMB wrapper.  Each CosmoRec setting was evaluated in a fresh Python process because native CosmoRec caches global initialization.

Fit objective:

- fit window: `z = 800..2500`
- observables fitted: fractional `x_e`
- residuals: uniform fractional residual plus a modest visibility-weighted term, with a small prior to keep the parameters close to the starting values
- fixed unless stated: `RECFAST_fudge = 1.125`
- six-parameter fit: Gaussian H parameters only, `RECFAST_fudge_He = 0.86`
- seven-parameter fit: same plus free `RECFAST_fudge_He`

This is **not** the old note's final `visibility_regularized_plus_uniform_minimax` metric.  The temp script used:

- `frac = (x_e^RECFAST - x_e^ref) / x_e^ref`
- residual block 1: `frac`
- residual block 2: `0.35 * sqrt(visibility / mean(visibility)) * frac`
- residual block 3: `0.02 * (params - start) / prior_scale`

So the fit below is a conservative mixed uniform/visibility fit for the current HEAD RECFAST parameter surface.  It is closer in spirit to the old note's "do not chase pure visibility weighting too hard" conclusion than to the named `visibility_regularized_plus_uniform_minimax` fit, but it is not numerically the same objective.

Highest valid shell reference tested:

- `cosmorec_accuracy = 6`
- `cosmorec_nshells = 10`
- `cosmorec_n_s_2gamma = 8`
- `cosmorec_n_s_raman = 7`
- `cosmorec_diff_iteration_max = 2`

Reproduction artifacts for this reference:

- native-style CosmoRec settings file: `docs/changelog/recfast_cosmorec_head_fits/cosmorec_accuracy6_h10_reference_runfile.txt`
- exact CAMB-wrapper `x_e(z)` dump, with `T_b` included as the third column: `docs/changelog/recfast_cosmorec_head_fits/cosmorec_accuracy6_h10_reference_xe.txt`

The exact reference used in the fits is the CAMB wrapper call, which passes CAMB's sampled `H(z)` directly to CosmoRec.  The native CosmoRec runfile records the equivalent final shell/radiative-transfer settings, but unmodified native CosmoRec does not expose every batch-wrapper detail in `parameters.dat`.

In a direct native-console check, the closest output to the CAMB wrapper reference is now the native final `.trans` file, with fractional `x_e` RMS `4.12e-5` and max `1.31e-4` over `z = 800..2500`.  The intermediate one-diffusion-iteration file `CosmoRec...Xe.it_1.accuracy6_h10.dat` is farther away: RMS `2.84e-4`, max `7.18e-4` over the same range.

Native-console difference plot and metrics:

- `docs/changelog/recfast_cosmorec_head_fits/xe_native_cosmorec_trans_vs_camb_accuracy6_h10.png`
- `docs/changelog/recfast_cosmorec_head_fits/xe_native_cosmorec_trans_vs_camb_accuracy6_h10_summary.json`

The CAMB wrapper exposes this source hook as `CosmoRec.diff_iteration_max` / `cosmorec_diff_iteration_max`.  With the companion CosmoRec patch that reads `runpars[4]`, all diffusion-on CosmoRec references in this note set `diff_iteration_max = 2`, returning the final `.trans`-equivalent table programmatically through CAMB's usual recombination spline without reading the native output files.

The local CosmoRec build only accepts hydrogen `nShells` in `{3, 4, 5, 10}`; trying `nShells = 12` exits natively.

## CosmoRec Convergence

Fractional `x_e` differences over `z = 400..3000`, relative to the highest valid shell reference above:

| Reference | RMS | Max |
| --- | ---: | ---: |
| `cosmorec_accuracy6` | `6.07e-5` | `1.76e-4` |
| `cosmorec_h4_he3_diffusion` | `6.94e-5` | `1.76e-4` |
| `cosmorec_he_diffusion` | `3.25e-4` | `7.41e-4` |
| `cosmorec_h10_he5` | `3.31e-4` | `1.37e-3` |
| `cosmorec_h4` | `4.55e-4` | `1.38e-3` |
| `cosmorec_default` | `5.31e-4` | `1.38e-3` |
| `cosmorec_h10_he3` | `5.55e-4` | `1.46e-3` |
| `cosmorec_hyrec_like` | `9.99e-4` | `4.13e-3` |
| `hyrec` | `1.13e-3` | `4.51e-3` |

For these RECFAST refits, `cosmorec_accuracy = 6` is already close to the explicit H10/ntg8/nR7 reference.  Pushing the hydrogen shell count to 10 after the accuracy-6 preset changes `x_e` by only `6e-5` RMS in this range.

The convergence plots compare direct CosmoRec settings to the highest valid `accuracy6_h10` reference.  Current default RECFAST is included as a dashed black reference curve.  Each CosmoRec setting was evaluated in a separate Python process so the native CosmoRec global initialization could not reuse the first setting.

CosmoRec convergence plot files:

- `docs/changelog/recfast_cosmorec_head_fits/cosmorec_xe_convergence_vs_accuracy6_h10.png`
- `docs/changelog/recfast_cosmorec_head_fits/cosmorec_xe_convergence_vs_accuracy6_h10_with_hyrec.png`
- `docs/changelog/recfast_cosmorec_head_fits/cosmorec_cmb_convergence_lmax4000_vs_accuracy6_h10.png`
- `docs/changelog/recfast_cosmorec_head_fits/cosmorec_convergence_lmax4000_summary.json`
- `docs/changelog/recfast_cosmorec_head_fits/cosmorec_xe_convergence_with_hyrec_summary.json`

The `with_hyrec` version of the `x_e` plot uses one shared legend below both panels and adds direct HyRec as a dotted independent-reference curve.

Selected lensed-scalar CMB convergence metrics against direct CosmoRec `accuracy6_h10`, using default CAMB timesteps and `lmax = 4000`:

| Model | Max `|Delta TT / TT|`, `2000 < ell <= 4000` | Max `|Delta EE / EE|`, `2000 < ell <= 4000` | Max `|Delta TE / sqrt(TT EE)|`, `2000 < ell <= 4000` |
| --- | ---: | ---: | ---: |
| Current default RECFAST | `5.82e-4 @ 2009` | `6.68e-4 @ 2173` | `3.92e-4 @ 2578` |
| `cosmorec_default` | `1.64e-3 @ 3990` | `1.86e-3 @ 3862` | `5.44e-4 @ 3288` |
| `cosmorec_h4` | `4.09e-4 @ 3973` | `4.55e-4 @ 3853` | `1.40e-4 @ 3278` |
| `cosmorec_h4_he3_diffusion` | `3.59e-4 @ 3958` | `4.12e-4 @ 3846` | `1.18e-4 @ 3906` |
| `cosmorec_accuracy6` | `3.27e-4 @ 3955` | `3.78e-4 @ 3843` | `1.09e-4 @ 3270` |

Full-range SO `Delta chi2` values for the same direct CosmoRec convergence models, against direct CosmoRec `accuracy6_h10`, using `camb.check_accuracy.calculate_cmb_delta_chi2` with `T,E,B`, `L = 2..4000`, and `get_total_cls(..., CMB_unit="muK")`:

| Model vs direct CosmoRec `accuracy6_h10` | SO `Delta chi2`, `L=2..4000` | Largest per-`L` contribution |
| --- | ---: | ---: |
| Current default RECFAST | `0.636` | `5.03e-4 @ 2149` |
| `cosmorec_default` | `7.054` | `3.99e-3 @ 3355` |
| `cosmorec_h4` | `4.34e-1` | `2.46e-4 @ 3322` |
| `cosmorec_h10_he3` | `2.09e-2` | `1.28e-5 @ 3401` |
| `cosmorec_h10_he5` | `1.34e-1` | `9.23e-5 @ 3578` |
| `cosmorec_he_diffusion` | `7.643` | `4.37e-3 @ 3357` |
| `cosmorec_h4_he3_diffusion` | `3.41e-1` | `1.92e-4 @ 3331` |
| `cosmorec_accuracy6` | `2.92e-1` | `1.62e-4 @ 3326` |

The CMB convergence plot shows the same pattern as the `x_e` table: the default CosmoRec-level settings are not converged at high ell relative to `accuracy6_h10`, while `accuracy6` is already close to the explicit H10/ntg8/nR7 reference.

## Recommended HEAD Fit Against Highest CosmoRec Reference

Six-parameter Gaussian-only fit:

- `AGauss1 = -0.1409170221`
- `AGauss2 = 0.0755105195`
- `zGauss1 = 7.2836736667`
- `zGauss2 = 6.7301232681`
- `wGauss1 = 0.1751063244`
- `wGauss2 = 0.3340950663`
- `RECFAST_fudge = 1.125`
- `RECFAST_fudge_He = 0.86`

Seven-parameter fit including helium fudge:

- `AGauss1 = -0.1421246011`
- `AGauss2 = 0.0747533979`
- `zGauss1 = 7.2841557503`
- `zGauss2 = 6.7306248459`
- `wGauss1 = 0.1741263787`
- `wGauss2 = 0.3349828575`
- `RECFAST_fudge = 1.125`
- `RECFAST_fudge_He = 0.8472384364`

Fit-window metrics against the highest CosmoRec reference:

| Model | RMS | Visibility-weighted RMS | `theta_*` fraction |
| --- | ---: | ---: | ---: |
| Current defaults | `1.10e-3` | `6.38e-4` | `-1.13e-5` |
| Six-parameter fit | `1.05e-3` | `8.15e-5` | `+2.28e-6` |
| Seven-parameter fit | `4.91e-4` | `8.09e-5` | `+2.67e-6` |

Over the wider plotted range `z = 400..3000`, relative to direct CosmoRec with the same `accuracy6_h10` settings:

| Model | RMS `Delta x_e / x_e` | Max `|Delta x_e / x_e|` |
| --- | ---: | ---: |
| Current default RECFAST | `1.14e-3` | `3.31e-3` |
| Six-parameter fit | `1.08e-3` | `3.31e-3` |
| Seven-parameter He-fudge fit | `7.64e-4` | `2.24e-3` |

The six-parameter fit mainly improves the visibility-weighted part of the recombination history and `theta_*`; it does not reduce the broad `z = 400..3000` max because the largest residual is in the high-redshift/helium tail.  Freeing `RECFAST_fudge_He` is what materially improves the broad `x_e` plot.

## Spectra Against Direct CosmoRec

The same models were compared to direct CosmoRec with `cosmorec_accuracy = 6`, `cosmorec_nshells = 10`, `cosmorec_n_s_2gamma = 8`, `cosmorec_n_s_raman = 7`, and `cosmorec_diff_iteration_max = 2`, using default CAMB accuracy settings and lensed scalar spectra to `lmax = 3000`.

RMS differences for `ell >= 30`:

| Model | `Delta TT / TT` | `Delta EE / EE` | `Delta TE / sqrt(TT EE)` |
| --- | ---: | ---: | ---: |
| Current default RECFAST | `3.31e-4` | `3.40e-4` | `1.86e-4` |
| Six-parameter fit | `1.44e-4` | `1.90e-4` | `7.10e-5` |
| Seven-parameter He-fudge fit | `8.63e-5` | `7.51e-5` | `6.41e-5` |

RMS differences for `ell >= 500`:

| Model | `Delta TT / TT` | `Delta EE / EE` | `Delta TE / sqrt(TT EE)` |
| --- | ---: | ---: | ---: |
| Current default RECFAST | `3.54e-4` | `2.93e-4` | `1.88e-4` |
| Six-parameter fit | `1.51e-4` | `1.93e-4` | `4.74e-5` |
| Seven-parameter He-fudge fit | `8.43e-5` | `5.03e-5` | `3.38e-5` |

Plots and the diagnostic JSON are saved in:

- `docs/changelog/recfast_cosmorec_head_fits/xe_fractional_difference_vs_cosmorec_accuracy6_h10.png`
- `docs/changelog/recfast_cosmorec_head_fits/xe_abs_fractional_difference_vs_cosmorec_accuracy6_h10.png`
- `docs/changelog/recfast_cosmorec_head_fits/cmb_spectra_difference_vs_cosmorec_accuracy6_h10.png`
- `docs/changelog/recfast_cosmorec_head_fits/diagnostics_summary.json`

## Old-Note-Style Visibility-Regularized Minimax Check

The old named `visibility_regularized_plus_uniform_minimax` script is not in this repository, so this recheck used an explicit implementation of the idea rather than the exact old script:

- visibility term: `sqrt(visibility / mean(visibility)) * sign(frac) * max(|frac| - 1e-4, 0)`
- uniform soft-minimax term: a cubic residual in `frac`, so the least-squares cost applies a smooth high-power penalty to the largest uniform residuals
- small stabilization prior: `0.006 * (params - start) / prior_scale`

This is the relevant comparison for `check_accuracy`-style high-ell residuals.  Spectra were recomputed to `lmax = 4000`; CAMB returns spectra through `ell = 4050` after the usual padding, and the max tables below use `ell > 2000`.

Old-note-style CosmoRec-targeted seven-parameter fit:

- `AGauss1 = -0.1409344303`
- `AGauss2 = 0.0748146374`
- `zGauss1 = 7.2836018868`
- `zGauss2 = 6.7308403696`
- `wGauss1 = 0.1747358203`
- `wGauss2 = 0.3356371126`
- `RECFAST_fudge = 1.125`
- `RECFAST_fudge_He = 0.8444375595`

Old-note-style HyRec-targeted seven-parameter fit:

- `AGauss1 = -0.1407228723`
- `AGauss2 = 0.0730048298`
- `zGauss1 = 7.2834789409`
- `zGauss2 = 6.7334927006`
- `wGauss1 = 0.1744109693`
- `wGauss2 = 0.3395288309`
- `RECFAST_fudge = 1.125`
- `RECFAST_fudge_He = 0.8706252645`

Against direct CosmoRec `accuracy6_h10`, the old-note-style fit is **not a high-ell spectra win** over the mixed fit when both are trained against the two-diffusion-iteration reference:

| Model vs direct CosmoRec | Max `|Delta TT / TT|`, `ell > 2000` | Max `|Delta EE / EE|`, `ell > 2000` | Max `|Delta TE / sqrt(TT EE)|`, `ell > 2000` |
| --- | ---: | ---: | ---: |
| Current default RECFAST | `5.82e-4 @ 2009` | `6.68e-4 @ 2173` | `3.92e-4 @ 2578` |
| Mixed CosmoRec fit with He | `1.36e-4 @ 4005` | `7.45e-5 @ 3584` | `5.19e-5 @ 3939` |
| Vreg/minimax CosmoRec fit with He | `1.74e-4 @ 4008` | `9.76e-5 @ 3603` | `6.72e-5 @ 3633` |
| Direct HyRec | `1.31e-3 @ 4009` | `1.39e-3 @ 4050` | `4.12e-4 @ 3610` |

So for a CosmoRec-targeted replacement, the mixed fit remains preferable by this high-ell max-error criterion.

Against direct HyRec, however, the old-note-style HyRec fit is clearly better than the previous mixed HyRec fit on high-ell maxima:

| Model vs direct HyRec | Max `|Delta TT / TT|`, `ell > 2000` | Max `|Delta EE / EE|`, `ell > 2000` | Max `|Delta TE / sqrt(TT EE)|`, `ell > 2000` |
| --- | ---: | ---: | ---: |
| Current default RECFAST | `1.37e-3 @ 2648` | `1.52e-3 @ 3472` | `6.76e-4 @ 2593` |
| Mixed HyRec fit with He | `4.76e-4 @ 4016` | `4.20e-4 @ 3933` | `1.55e-4 @ 3941` |
| Vreg/minimax HyRec fit with He | `2.31e-4 @ 4003` | `1.73e-4 @ 3586` | `1.00e-4 @ 2679` |
| Direct CosmoRec `accuracy6_h10` | `1.30e-3 @ 4009` | `1.39e-3 @ 4050` | `4.12e-4 @ 3610` |

The broad `z = 400..3000` `x_e` plot metrics are also essentially tied between mixed and vreg/minimax for the matched target:

| Reference | Model | RMS `Delta x_e / x_e` | Max `|Delta x_e / x_e|` |
| --- | --- | ---: | ---: |
| Direct CosmoRec `accuracy6_h10` | Mixed CosmoRec fit with He | `7.64e-4` | `2.24e-3` |
| Direct CosmoRec `accuracy6_h10` | Vreg/minimax CosmoRec fit with He | `7.44e-4` | `2.25e-3` |
| Direct HyRec | Mixed HyRec fit with He | `9.95e-4` | `2.63e-3` |
| Direct HyRec | Vreg/minimax HyRec fit with He | `9.92e-4` | `2.62e-3` |

## SO CMB Delta Chi2

The CMB `Delta chi2` numbers below use `camb.check_accuracy.calculate_cmb_delta_chi2` directly, with the built-in SO noise settings:

- fields: `T,E,B`
- `noise_muK_arcmin_T = 2.6`
- `noise_muK_arcmin_P = 3.67`
- `fwhm_arcmin = 1.4`
- `fsky = 0.6`
- spectra: `get_total_cls(..., CMB_unit="muK")`

The main column is the check-accuracy-style total over `L = 2..4000`.  The `L > 2000` column repeats the same calculation with `lmin = 2001` to isolate the high-ell contribution.

| Model vs direct CosmoRec `accuracy6_h10` | SO `Delta chi2`, `L=2..4000` | SO `Delta chi2`, `L=2001..4000` | Largest per-`L` contribution |
| --- | ---: | ---: | ---: |
| Current default RECFAST | `0.636` | `0.322` | `5.03e-4 @ 2149` |
| Mixed CosmoRec fit with He | `4.28e-2` | `2.88e-2` | `2.25e-5 @ 3043` |
| Vreg/minimax CosmoRec fit with He | `6.44e-2` | `4.68e-2` | `3.29e-5 @ 3055` |
| Direct HyRec | `3.620` | `3.294` | `2.20e-3 @ 3558` |

For the CosmoRec target, the mixed fit is still better by this SO metric, but both fitted RECFAST histories are well below current default RECFAST.

| Model vs direct HyRec | SO `Delta chi2`, `L=2..4000` | SO `Delta chi2`, `L=2001..4000` | Largest per-`L` contribution |
| --- | ---: | ---: | ---: |
| Current default RECFAST | `5.962` | `4.879` | `3.14e-3 @ 2641` |
| Mixed HyRec fit with He | `4.36e-1` | `3.79e-1` | `2.56e-4 @ 3061` |
| Vreg/minimax HyRec fit with He | `1.16e-1` | `8.87e-2` | `6.82e-5 @ 3043` |
| Direct CosmoRec `accuracy6_h10` | `3.622` | `3.296` | `2.20e-3 @ 3558` |

For the HyRec target, the vreg/minimax fit is preferred by both the high-ell max residuals and the SO `Delta chi2` metric.

## HyRec vs CosmoRec Parameter-Bias Projection

As a rough cosmological-impact estimate, I projected the direct HyRec versus direct CosmoRec `accuracy6` spectral residual onto `N_eff` and `n_s` templates using `camb.check_accuracy.calculate_cmb_delta_chi2`, SO noise, `T,E,B`, and `L = 2..4000`.

This is a fixed-background Fisher/parabolic estimate: all parameters other than the listed parameter(s) are held fixed, including `H0`, `ombh2`, `omch2`, `tau`, and `A_s`.  A full cosmological fit with acoustic-scale and matter-density degeneracies marginalized would be needed for final parameter biases.

For the concrete case where the fitted model uses HyRec but the "data" are direct CosmoRec `accuracy6`, the baseline mismatch is `Delta chi2 = 1.891`.  The effective fixed-background shifts are:

| Fit template | Best-fit shift | Fixed-parameter SO sigma | Shift/sigma | Residual `Delta chi2` |
| --- | ---: | ---: | ---: | ---: |
| `N_eff` only | `-6.13e-4` | `2.04e-3` | `-0.30` | `1.89` |
| `n_s` only | `-4.27e-4` | `3.18e-4` | `-1.34` | `0.089` |
| joint `N_eff,n_s` | `Delta N_eff = -1.64e-4`, `Delta n_s = -4.40e-4` | `sigma(N_eff)=2.07e-3`, `sigma(n_s)=3.22e-4` | `-0.08`, `-1.37` | `~0` |

The reverse direction, fitting direct CosmoRec `accuracy6` to HyRec "data", gives approximately the opposite tilt shift: one-parameter `Delta n_s = +4.29e-4`, joint `Delta n_s = +4.17e-4`.  The `N_eff` shift is smaller and more likelihood-direction dependent: one-parameter `Delta N_eff = +2.39e-4`, joint `Delta N_eff = -2.14e-4`.

So in this restricted SO projection, the HyRec/CosmoRec `accuracy6` difference still looks more like a small tilt bias, `|Delta n_s| ~= 4.3e-4`, than an `N_eff` bias.  In the fixed-background SO metric used here this is about `1.3 sigma` because only `n_s` is allowed to move; it should not be read as a marginalized cosmological constraint.

If TT information is removed, or mostly downweighted, the same fixed-background projection gives a smaller effective tilt bias:

| SO variant, HyRec model vs CosmoRec `accuracy6` data | Baseline `Delta chi2` | `N_eff`-only shift | `n_s`-only shift | Joint shift | Joint residual `Delta chi2` |
| --- | ---: | ---: | ---: | ---: | ---: |
| `E,B` only | `0.274` | `Delta N_eff = -2.79e-4` | `Delta n_s = -3.70e-4` | `Delta N_eff = -1.47e-4`, `Delta n_s = -3.73e-4` | `0.008` |
| `T,E,B`, 10x T noise | `0.469` | `Delta N_eff = -2.77e-4` | `Delta n_s = -3.58e-4` | `Delta N_eff = -1.60e-4`, `Delta n_s = -3.60e-4` | `0.021` |

So excluding or downweighting TT reduces the effective tilt bias modestly, from `|Delta n_s| ~= 4.3e-4` to about `3.6e-4..3.7e-4`.  The `N_eff` projection remains weak: the `N_eff`-only fit leaves most of the residual `Delta chi2` in both temperature-light variants.

Varying only `YHe` is also an efficient one-parameter absorber of the HyRec/CosmoRec `accuracy6` difference:

| SO variant, HyRec model vs CosmoRec `accuracy6` data | Baseline `Delta chi2` | Best-fit `Delta YHe` | Fixed-parameter SO sigma | Shift/sigma | Exact `Delta chi2` at shifted `YHe` |
| --- | ---: | ---: | ---: | ---: | ---: |
| `T,E,B` | `1.891` | `+3.75e-4` | `2.84e-4` | `+1.32` | `0.078` |
| `E,B` only | `0.274` | `+3.53e-4` | `7.39e-4` | `+0.48` | `0.038` |
| `T,E,B`, 10x T noise | `0.469` | `+3.43e-4` | `5.48e-4` | `+0.63` | `0.063` |

The reverse direction gives approximately the opposite sign: fitting CosmoRec `accuracy6` to HyRec "data" prefers `Delta YHe = -3.85e-4` for full SO, `-3.62e-4` for `E,B` only, and `-3.53e-4` for 10x T noise.

The new plot and diagnostic files include both direct HyRec and the HyRec-targeted RECFAST fits.  The `x_e` comparison plot also includes the native CosmoRec console `.trans` output as a black dot-dashed curve, using the same native history checked above:

- `docs/changelog/recfast_cosmorec_head_fits/xe_fractional_difference_cosmorec_hyrec_refs.png`
- `docs/changelog/recfast_cosmorec_head_fits/xe_fractional_difference_cosmorec_hyrec_refs_with_native_trans_summary.json`
- `docs/changelog/recfast_cosmorec_head_fits/cmb_spectra_lmax4000_vs_direct_cosmorec_accuracy6_h10.png`
- `docs/changelog/recfast_cosmorec_head_fits/cmb_spectra_lmax4000_vs_direct_hyrec.png`
- `docs/changelog/recfast_cosmorec_head_fits/vreg_lmax4000_summary.json`
- `docs/changelog/recfast_cosmorec_head_fits/so_chi2_lmax4000_summary.json`
- `docs/changelog/recfast_cosmorec_head_fits/hyrec_cosmorec_accuracy6_param_bias_so_lmax4000.json`
- `docs/changelog/recfast_cosmorec_head_fits/hyrec_cosmorec_accuracy6_param_bias_so_temperature_variants_lmax4000.json`
- `docs/changelog/recfast_cosmorec_head_fits/hyrec_cosmorec_accuracy6_yhe_bias_so_lmax4000.json`

## Direct HyRec/CosmoRec and Hybrid-History Check

I also made a focused comparison of HyRec against CosmoRec `accuracy6_h10` to `lmax = 4000`.  This plot uses default CAMB timesteps and lensed scalar spectra, with CosmoRec `accuracy=6`, `nShells=10`, `n_s_2gamma=8`, `n_s_raman=7`, and `diff_iteration_max=2` as the zero line:

- `docs/changelog/recfast_cosmorec_head_fits/cmb_spectra_cosmorec_hyrec_hybrid_he_switch_lensed_lmax4000.png`

For `2000 < ell <= 4000`, the maxima are:

| Model vs direct CosmoRec `accuracy6_h10` | Max `|Delta TT / TT|` | Max `|Delta EE / EE|` | Max `|Delta TE / sqrt(TT EE)|` |
| --- | ---: | ---: | ---: |
| Direct HyRec | `1.31e-3 @ ell=4000` | `1.38e-3 @ ell=4000` | `4.12e-4 @ ell=3610` |
| CosmoRec He history, then HyRec H evolution | `1.91e-4 @ ell=2648` | `2.29e-4 @ ell=2825` | `8.53e-5 @ ell=2600` |

The corresponding `x_e` plot includes a diagnostic hybrid history: CosmoRec above `z = 1600`, then HyRec's lower-redshift hydrogen-evolution shape below that switch, with a continuity offset applied at the switch.  The vertical red line marks the switch redshift:

- `docs/changelog/recfast_cosmorec_head_fits/xe_cosmorec_hyrec_hybrid_he_switch.png`
- `docs/changelog/recfast_cosmorec_head_fits/hybrid_cosmorec_he_hyrec_h_validated_summary.json`

Over `z = 400..3000`, relative to direct CosmoRec `accuracy6_h10`, the direct HyRec history has RMS `1.13e-3` and max `4.51e-3`, while the hybrid history has RMS `1.54e-4` and max `3.92e-4`.

The hybrid spectra use a native temp-worktree recombination class, not a precomputed recombination table.  The class initializes direct CosmoRec and direct HyRec with the same CAMB state, applies the switch at `z = 1600`, and uses the small continuity offsets `Delta x_e = 7.56e-8` and `Delta T_m = 3.25e-6` at the switch.

## Freeing `RECFAST_fudge`

The fits above keep the hydrogen `RECFAST_fudge` fixed at `1.125`.  I checked whether freeing it helps, especially in the lower-redshift hydrogen tail.  This check uses only the direct CosmoRec high-shell reference with `diff_iteration_max = 2`, matching the convention used elsewhere in this note.

- `fixed_h_mixed7`: the previous seven-parameter mixed fit, with `RECFAST_fudge = 1.125`
- `h_free_mixed`: refit the mixed objective with `RECFAST_fudge` and `RECFAST_fudge_He` both free
- `h_scan_best_z400_3000`: keep the previous seven parameters fixed and scan only `RECFAST_fudge`; the broad `z = 400..3000` RMS optimum is `RECFAST_fudge = 1.12725`

The simple one-dimensional scan confirms that `RECFAST_fudge` is the right lever for the lower-redshift residual: increasing it from `1.125` to `1.12725` reduces the `z = 400..800` RMS by about a factor of `3.5` and the broad `z = 400..3000` RMS by about `40%`.  The cost is worse high-ell spectra, especially EE/TE.  Letting all parameters move in the mixed objective gives only a small fudge shift (`1.12532`), improves `x_e` less, and still worsens high-ell spectra.

| Reference | Model | `RECFAST_fudge` | `z400..800` RMS / max | `z400..3000` RMS / max | Max high-ell TT / EE / TE |
| --- | --- | ---: | ---: | ---: | ---: |
| diffusion iter2 | fixed mixed | `1.12500` | `1.65e-3` / `2.24e-3` | `7.56e-4` / `2.24e-3` | `1.74e-4` / `1.12e-4` / `6.28e-5` |
| diffusion iter2 | H-only scan | `1.12725` | `4.73e-4` / `8.03e-4` | `4.44e-4` / `1.30e-3` | `1.83e-4` / `3.03e-4` / `1.15e-4` |
| diffusion iter2 | H-free mixed | `1.12532` | `1.44e-3` / `2.00e-3` | `6.68e-4` / `2.00e-3` | `3.16e-4` / `2.61e-4` / `1.05e-4` |

So freeing or retuning `RECFAST_fudge` is worthwhile only if the priority is the low-redshift `x_e` tail itself.  For the CMB spectra/check-accuracy target, keeping `RECFAST_fudge = 1.125` remains the better compromise.

I also tried a broad joint refit that includes `z = 400..2500`, but downweights the low-redshift tail where the visibility is small.  The objective used the correct cached direct-CosmoRec `plot_xe` history below `z = 800` rather than extrapolating the shorter `fit_history` grid:

- uniform fractional residual for `z >= 800`
- a full-window visibility term `0.55 * sqrt(visibility / max(visibility)) * Delta x_e / x_e`
- a small low-redshift guard term, either `0`, `0.05`, or `0.10` times the uniform fractional residual for `z < 800`
- the same small parameter prior as the other hydrogen-fudge checks

This is a saner version of "include the low-z tail, but do not let it dominate."  It improves the plotted tail modestly, but it is not a spectra win:

| Model | Low-z guard | `RECFAST_fudge` | `RECFAST_fudge_He` | `z400..800` RMS / max | Max high-ell TT / EE / TE |
| --- | ---: | ---: | ---: | ---: | ---: |
| fixed mixed | -- | `1.12500` | `0.84662` | `1.65e-3` / `2.24e-3` | `1.74e-4` / `1.12e-4` / `6.28e-5` |
| joint broad | `0.00` | `1.12519` | `0.84370` | `1.55e-3` / `2.12e-3` | `4.48e-4` / `4.07e-4` / `1.48e-4` |
| joint broad | `0.05` | `1.12543` | `0.84370` | `1.40e-3` / `1.95e-3` | `4.40e-4` / `4.01e-4` / `1.47e-4` |
| joint broad | `0.10` | `1.12592` | `0.84370` | `1.10e-3` / `1.61e-3` | `4.25e-4` / `3.88e-4` / `1.44e-4` |

So the visibility-downweighted joint fit confirms the same conclusion: once the low-z tail is prevented from dominating, the optimizer only moves `RECFAST_fudge` slightly, and the accompanying helium/Gaussian shifts degrade the `ell > 2000` spectra more than they help the broad `x_e` plot.  I would keep the fixed-hydrogen mixed fit as the CosmoRec-targeted CMB compromise.

Diagnostics:

- `docs/changelog/recfast_cosmorec_head_fits/xe_recfast_h_fudge_free_lowz_check.png`
- `docs/changelog/recfast_cosmorec_head_fits/recfast_h_fudge_free_lowz_check_lmax4000.json`
- `docs/changelog/recfast_cosmorec_head_fits/xe_recfast_joint_visibility_downweighted.png`
- `docs/changelog/recfast_cosmorec_head_fits/recfast_joint_visibility_downweighted_lmax4000.json`

## Early-Dark-Energy Stress Test

I also checked the fitted `recfast_cosmorec` parameters on a substantially non-standard background, using CAMB's `EarlyQuintessence` model with:

- `thetastar = 0.0104085`
- `n = 3`, `theta_i = 3.1`
- `fde_zc = 0.1`
- two peak-redshift cases: requested `zc = 3500` and `zc = 2000`
- solved values: `zc = 3499.9383`, `fde_zc = 0.09999999`; and `zc = 1999.9620`, `fde_zc = 0.10000009`

The reference is direct CosmoRec with the same settings used above: `accuracy = 6`, `n_shells = 10`, `n_shells_hei = 5`, `n_s_2gamma = 8`, `n_s_raman = 7`, `flag_hi_absorption = 3`, and `diff_iteration_max = 2`.

The `x_e` plot compares this direct CosmoRec reference to the `recfast_cosmorec` fudge fit, with `thetastar` fixed separately for each recombination model.  For context it also includes current default RECFAST as a dashed curve, also with the same fixed `thetastar` input.  The blue dash-dot curve shows the direct CosmoRec EDE-vs-LCDM change at fixed `thetastar`, `(x_e^{EDE} - x_e^{LCDM}) / x_e^{EDE}`, to show the size of the actual non-standard-background signal.

- `docs/changelog/recfast_cosmorec_head_fits/xe_ede_planck_theta_recfast_cosmorec_vs_direct_cosmorec.png`
- `docs/changelog/recfast_cosmorec_head_fits/xe_ede_planck_theta_recfast_cosmorec_vs_direct_cosmorec_summary.json`

Relative to direct CosmoRec over `z = 400..5000`:

| EDE case | Model | RMS `Delta x_e/x_e` | Max `|Delta x_e/x_e|` | `z` at max |
| --- | --- | ---: | ---: | ---: |
| `zc ~= 3500` | `recfast_cosmorec`, same `theta_*` input | `5.87e-4` | `2.31e-3` | `572.6` |
| `zc ~= 3500` | Current default RECFAST, same `theta_*` input | `8.54e-4` | `3.15e-3` | `1923.1` |
| `zc ~= 3500` | Direct CosmoRec EDE minus LCDM, same `theta_*` input | `3.34e-3` | `9.02e-3` | `845.0` |
| `zc ~= 2000` | `recfast_cosmorec`, same `theta_*` input | `5.99e-4` | `2.34e-3` | `572.6` |
| `zc ~= 2000` | Current default RECFAST, same `theta_*` input | `8.50e-4` | `3.05e-3` | `1919.3` |
| `zc ~= 2000` | Direct CosmoRec EDE minus LCDM, same `theta_*` input | `5.73e-3` | `1.59e-2` | `1010.0` |

Restricting to the hydrogen recombination window `z = 800..1600`, the fitted `recfast_cosmorec` history has RMS/max `1.60e-4`/`3.24e-4` for `zc ~= 3500`, and `2.80e-4`/`5.55e-4` for `zc ~= 2000`.  The corresponding direct CosmoRec EDE-vs-LCDM changes are much larger: RMS/max `5.82e-3`/`9.02e-3` and `1.01e-2`/`1.59e-2`.  This suggests the CosmoRec-targeted RECFAST fudge fit extrapolates reasonably well to early-dark-energy models with a 10% peak pre-recombination density fraction, including a lower peak redshift close to recombination.

## Other Reference Fits

Seven-parameter fits, useful as approximate RECFAST settings for each reference:

| Reference | `AGauss1` | `AGauss2` | `zGauss1` | `zGauss2` | `wGauss1` | `wGauss2` | `RECFAST_fudge_He` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `cosmorec_default` | `-0.1354282830` | `0.0714804967` | `7.2856502616` | `6.7265922818` | `0.1753514305` | `0.3316332603` | `0.8520965722` |
| `cosmorec_hyrec_like` | `-0.1404711869` | `0.0741166653` | `7.2844210127` | `6.7299819125` | `0.1744142667` | `0.3346766357` | `0.8635739044` |
| `cosmorec_h4` | `-0.1405247251` | `0.0741190677` | `7.2844578347` | `6.7299708000` | `0.1744042801` | `0.3346438448` | `0.8520965685` |
| `cosmorec_h10_he3` | `-0.1421189108` | `0.0747553747` | `7.2841486601` | `6.7306304256` | `0.1741289943` | `0.3349892250` | `0.8503361937` |
| `cosmorec_h10_he5` | `-0.1421145992` | `0.0747538699` | `7.2841496543` | `6.7306258279` | `0.1741330708` | `0.3349812749` | `0.8526476548` |
| `cosmorec_he_diffusion` | `-0.1354310099` | `0.0714768465` | `7.2856554215` | `6.7265896976` | `0.1753539240` | `0.3316282642` | `0.8494038671` |
| `cosmorec_h4_he3_diffusion` | `-0.1405339308` | `0.0741140807` | `7.2844604008` | `6.7299765780` | `0.1744007653` | `0.3346566698` | `0.8476706714` |
| `cosmorec_accuracy6` | `-0.1405330191` | `0.0741184911` | `7.2844615447` | `6.7299724247` | `0.1744032874` | `0.3346467657` | `0.8472384664` |
| `cosmorec_accuracy6_h10` | `-0.1421246011` | `0.0747533979` | `7.2841557503` | `6.7306248459` | `0.1741263787` | `0.3349828575` | `0.8472384364` |
| `hyrec` | `-0.1417072428` | `0.0733978462` | `7.2841809376` | `6.7316088397` | `0.1734480042` | `0.3373321220` | `0.8658491370` |

Full JSON outputs are in `/tmp/camb-recfast-fit-head-NofB30/build/recfast_fudge_fits/summary.json`.
