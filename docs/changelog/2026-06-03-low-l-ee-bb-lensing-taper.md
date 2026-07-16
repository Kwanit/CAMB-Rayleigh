# low-l EE lensing taper: regularize fast large-scale polarization

- Date / issue id: 2026-06-03, modern accuracy `params_tau_0p11_lmax3000.ini` lensed scalar BB.
- Symptom & reproducer: `python -m camb.check_accuracy fortran/testfiles/modern_accuracy_inis/params_tau_0p11_lmax3000.ini` failed `BB 2 <= L <= 3000`, max abs `0.0441 > 0.02` at `L=2`.
- Root cause: the optimized full-sky short-range lensing correction reuses the low-L reionization-bump EE spectrum in the BB correction; that physical contribution is negligible for large-scale lensed BB but makes the truncated BB integral sensitive at the lowest multipoles.
- Full physical extent: scalar polarization lensing corrections in the fast `CorrFuncFullSkyImpl` path when `AccurateBB=F`. The stored unlensed EE spectrum added back to lensed EE is unchanged. `AccurateBB=T` uses the direct full-range path and is left untapered.
- Locality / rejected alternatives: broad `LensRangeBoost` floors fixed the target but added runtime. Tapering TT/TE/EE together did not help the high-L EE neutrino rows and did not replace the support boosts. Keeping a separate BB-only EE copy fixed the symptom but split the polarization correction convention unnecessarily. The accepted fix tapers only the EE source used by the fast polarization correction.
- Change: `fortran/lensing.f90` adds a shared smoothstep `LowLEELensingTaper(l)` that is zero at `L=2` and reaches unity at `L=20`; the fast path multiplies the computed low-L EE lensing-correction source by this taper. The high-L template extension is left untapered.
- Before/after: the high-tau BB row changed from `0.0441 FAIL at L=2` to `0.003671 OK at L=5`.
- Acceptance matrix: rebuilt the main tree with `python setup.py make`. Targeted `check_accuracy` runs passed for `params_tau_0p11_lmax3000.ini`, `params_tau_0p03_lmax3000.ini`, `params_base_lmax4000.ini`, and the high-L massive-neutrino guard `params_three_deg_omnuh2_0p02400_lmax4000_lpa6.ini`.
- AccurateBB check: no modern accuracy ini sets `accurate_BB=T`; the target checks still compare the fast standard run against the direct AccurateBB reference. The direct path is intentionally unchanged, because its full angular range does not have the truncation instability.
- Runtime impact: no sampling/range boost added. Representative high-tau main-tree run after the final fix used `standard cpu=6.77s`, `boosted cpu=37.72s`.
- Checker verdict: local audit accepted the EE-only fast-path scope; sub-agent checker was not used in this session.
