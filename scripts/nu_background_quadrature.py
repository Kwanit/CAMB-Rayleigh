"""Validate the live CAMB thermal-neutrino background fits against direct quadrature.

This is a line-for-line numerical reproduction of the branches and coefficients
in ``fortran/massive_neutrinos.f90`` for ``ThermalNuBackground_rho_P``,
``ThermalNuBackground_rho``, and ``ThermalNuBackground_drho``. It compares
them with converged Fermi-Dirac integrals and retains the removed spline-table
path only as an accuracy baseline.
"""

import argparse
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.special import zeta

FERMI_DIRAC_CONST = 7.0 * np.pi**4 / 120.0
CONST2 = 5.0 / (7.0 * np.pi**2)
ZETA3 = float(zeta(3.0))
ZETA5 = float(zeta(5.0))
ZETA7 = float(zeta(7.0))

NU_RHOP_AM_MIN = 0.42
NU_RHOP_AM_MAX = 70.0
NU_RHOP_FIT_C = np.sqrt(NU_RHOP_AM_MIN * NU_RHOP_AM_MAX)
NU_RHOP_FIT_INV_ZMAX = (NU_RHOP_AM_MAX + NU_RHOP_FIT_C) / (NU_RHOP_AM_MAX - NU_RHOP_FIT_C)
NU_FIT_RHO_SCALE = 3.0 * ZETA3 / (2.0 * FERMI_DIRAC_CONST)
NU_FIT_P_DENOM_SCALE = FERMI_DIRAC_CONST / ((900.0 / 120.0) * ZETA5)

NU_FIT_RHO_C = np.array(
    [
        7.499710425926954249e-01,
        1.188066888581748443e-01,
        1.545077526625578401e-01,
        -8.378996019212282820e-02,
        2.123075360029347963e-02,
        -2.544964203886362908e-03,
    ]
)
NU_FIT_P_C = np.array(
    [
        9.283417043601330798e-01,
        3.156134695375645838e-01,
        -2.435843189264263742e-01,
        -2.509793530357270000e-02,
        4.067958090718486880e-02,
        2.378742143669706002e-03,
        -1.965863329520473827e-03,
    ]
)
NU_LOW_RHO_C = np.array([-1.872929199948838302e-02, 1.831728681837861694e-02])
NU_LOW_P_C = np.array([1.667044572156534468e-02, -2.277471498716537868e-02])

DRHO_AM_MIN = 0.3
DRHO_AM_MAX = 70.0
DRHO_FIT_C = np.sqrt(DRHO_AM_MIN * DRHO_AM_MAX)
DRHO_FIT_INV_ZMAX = (DRHO_AM_MAX + DRHO_FIT_C) / (DRHO_AM_MAX - DRHO_FIT_C)
NU_FIT_D_C = np.array(
    [
        5.800497880734221123e-01,
        1.801974470841860576e-01,
        -1.602709049602861757e-01,
        3.167230341650149189e-02,
        1.056591590995091534e-02,
        -3.877662402660121861e-03,
    ]
)

LEGACY_AM_MIN = 0.3
LEGACY_AM_MAX = 70.0
LEGACY_NRHOPN = 1400
LEGACY_AM_MINP = LEGACY_AM_MIN + LEGACY_AM_MAX / (LEGACY_NRHOPN - 1) * 1.01
LEGACY_AM_MAXP = LEGACY_AM_MAX * 0.9
LEGACY_QMAX = 30.0
LEGACY_NQ = 100


@dataclass(frozen=True)
class ErrorSummary:
    max_relative_error: float
    rms_relative_error: float
    am_at_max_error: float


@dataclass(frozen=True)
class LegacySplineTable:
    am_grid: np.ndarray
    rhonu: np.ndarray
    pnu: np.ndarray
    drho: np.ndarray
    dp: np.ndarray
    ddrho: np.ndarray


def horner_ascending(x: float, coefficients: np.ndarray) -> float:
    """Evaluate a polynomial with coefficients in increasing order."""
    value = 0.0
    for coefficient in coefficients[::-1]:
        value = coefficient + x * value
    return float(value)


def fit_coordinate(am: float, fit_c: float, fit_inv_zmax: float) -> float:
    return ((am - fit_c) / (am + fit_c)) * fit_inv_zmax


def thermal_nu_background_rho_p(am: float) -> tuple[float, float]:
    """Reproduce Fortran ``ThermalNuBackground_rho_P`` at ``am`` exactly."""
    if am <= NU_RHOP_AM_MIN:
        am2 = am**2
        am4 = am2**2
        rhonu = 1.0 + CONST2 * am2 + am4 * (NU_LOW_RHO_C[0] + am2 * NU_LOW_RHO_C[1])
        pnu = (1.0 - CONST2 * am2) / 3.0 + am4 * (NU_LOW_P_C[0] + am2 * NU_LOW_P_C[1])
        return rhonu, pnu
    if am >= NU_RHOP_AM_MAX:
        rhonu = (
            3.0
            / (2.0 * FERMI_DIRAC_CONST)
            * (ZETA3 * am + ((15.0 * ZETA5) / 2.0 - 945.0 * ZETA7 / (16.0 * am**2)) / am)
        )
        pnu = (900.0 / 120.0) / FERMI_DIRAC_CONST * (ZETA5 - 63.0 * ZETA7 / (4.0 * am**2)) / am
        return rhonu, pnu

    z = fit_coordinate(am, NU_RHOP_FIT_C, NU_RHOP_FIT_INV_ZMAX)
    rhonu = (1.0 + NU_FIT_RHO_SCALE * am) * horner_ascending(z, NU_FIT_RHO_C)
    pnu = horner_ascending(z, NU_FIT_P_C) / (1.0 + NU_FIT_P_DENOM_SCALE * am)
    return rhonu, pnu


def thermal_nu_background_rho(am: float) -> float:
    """Reproduce Fortran ``ThermalNuBackground_rho`` at ``am`` exactly."""
    if am <= NU_RHOP_AM_MIN:
        am2 = am**2
        am4 = am2**2
        return 1.0 + CONST2 * am2 + am4 * (NU_LOW_RHO_C[0] + am2 * NU_LOW_RHO_C[1])
    if am >= NU_RHOP_AM_MAX:
        return (
            3.0
            / (2.0 * FERMI_DIRAC_CONST)
            * (ZETA3 * am + ((15.0 * ZETA5) / 2.0 - 945.0 * ZETA7 / (16.0 * am**2)) / am)
        )

    z = fit_coordinate(am, NU_RHOP_FIT_C, NU_RHOP_FIT_INV_ZMAX)
    return (1.0 + NU_FIT_RHO_SCALE * am) * horner_ascending(z, NU_FIT_RHO_C)


def thermal_nu_background_drho(am: float, adotoa: float = 1.0) -> float:
    """Reproduce Fortran ``ThermalNuBackground_drho`` at ``am`` exactly."""
    if am <= DRHO_AM_MIN:
        am2 = am**2
        return am2 * (2.0 * CONST2 + am2 * (0.04399706676 * np.log(am) - 0.002970400378 - 0.029331377855 * am)) * adotoa
    if am >= DRHO_AM_MAX:
        return (
            3.0
            / (2.0 * FERMI_DIRAC_CONST)
            * (ZETA3 * am + (-(15.0 * ZETA5) / 2.0 + 2835.0 * ZETA7 / (16.0 * am**2)) / am)
            * adotoa
        )

    z = fit_coordinate(am, DRHO_FIT_C, DRHO_FIT_INV_ZMAX)
    d_fit = horner_ascending(z, NU_FIT_D_C)
    return am**2 / (1.0 + 2.0 * am) * d_fit * adotoa


def fermi_dirac_weight(q: float) -> float:
    return np.exp(-q) if q > 40.0 else 1.0 / (np.exp(q) + 1.0)


def exact_rho_p(am: float) -> tuple[float, float]:
    """Compute converged reference density and pressure integrals."""
    rho, _ = quad(
        lambda q: q**3 * fermi_dirac_weight(q) * np.sqrt(1.0 + (am / q) ** 2),
        0.0,
        np.inf,
        epsabs=1e-13,
        epsrel=1e-13,
        limit=300,
    )
    pressure, _ = quad(
        lambda q: q**3 * fermi_dirac_weight(q) / np.sqrt(1.0 + (am / q) ** 2),
        0.0,
        np.inf,
        epsabs=1e-13,
        epsrel=1e-13,
        limit=300,
    )
    return rho / FERMI_DIRAC_CONST, pressure / (3.0 * FERMI_DIRAC_CONST)


def exact_drho(am: float, adotoa: float = 1.0) -> float:
    """Compute the exact ``a * d rho / da * adotoa`` reference integral."""
    value, _ = quad(
        lambda q: am**2 * q**2 * fermi_dirac_weight(q) / np.sqrt(q**2 + am**2),
        0.0,
        np.inf,
        epsabs=1e-13,
        epsrel=1e-13,
        limit=300,
    )
    return value * adotoa / FERMI_DIRAC_CONST


def splint_uniform(values: np.ndarray) -> float:
    dyn = (11.0 * values[-1] - 18.0 * values[-2] + 9.0 * values[-3] - 2.0 * values[-4]) / 6.0
    return float(0.5 * (values[0] + values[-1]) - dyn / 12.0 + np.sum(values[1:-1]))


def splini(size: int) -> np.ndarray:
    g = np.empty(size)
    g[0] = 0.0
    for index in range(1, size):
        g[index] = 1.0 / (4.0 - g[index - 1])
    return g


def splder(values: np.ndarray, g: np.ndarray) -> np.ndarray:
    size = len(values)
    f = np.empty(size)
    f[0] = (-10.0 * values[0] + 15.0 * values[1] - 6.0 * values[2] + values[3]) / 6.0
    f[-1] = (10.0 * values[-1] - 15.0 * values[-2] + 6.0 * values[-3] - values[-4]) / 6.0
    for index in range(1, size - 1):
        f[index] = g[index] * (3.0 * (values[index + 1] - values[index - 1]) - f[index - 1])
    derivative = np.empty(size)
    derivative[-1] = f[-1]
    for index in range(size - 2, -1, -1):
        derivative[index] = f[index] - g[index] * derivative[index + 1]
    return derivative


def legacy_nu_rho_pres(am: float) -> tuple[float, float]:
    """Reproduce the removed 100-step ``nuRhoPres`` quadrature helper."""
    adq = LEGACY_QMAX / LEGACY_NQ
    rho_values = np.zeros(LEGACY_NQ + 1)
    pressure_values = np.zeros(LEGACY_NQ + 1)
    for index in range(1, LEGACY_NQ + 1):
        q = index * adq
        velocity = 1.0 / np.sqrt(1.0 + (am / q) ** 2)
        weight = adq * q**3 / (np.exp(q) + 1.0)
        rho_values[index] = weight / velocity
        pressure_values[index] = weight * velocity
    rhonu = (splint_uniform(rho_values) + rho_values[-1] / adq) / FERMI_DIRAC_CONST
    pnu = (splint_uniform(pressure_values) + pressure_values[-1] / adq) / (3.0 * FERMI_DIRAC_CONST)
    return rhonu, pnu


def legacy_low_rho_p(am: float) -> tuple[float, float]:
    am2 = am**2
    rhonu = 1.0 + am2 * (CONST2 + am2 * (0.01099926669 * np.log(am) - 0.003492416767 - 0.005866275571 * am))
    pnu = (1.0 + am2 * (-CONST2 + am2 * (-0.03299780009 * np.log(am) - 0.0005219952794 + 0.02346510229 * am))) / 3.0
    return rhonu, pnu


def legacy_high_rho_p(am: float) -> tuple[float, float]:
    rhonu = (
        3.0 / (2.0 * FERMI_DIRAC_CONST) * (ZETA3 * am + ((15.0 * ZETA5) / 2.0 - 945.0 * ZETA7 / (16.0 * am**2)) / am)
    )
    pnu = (900.0 / 120.0) / FERMI_DIRAC_CONST * (ZETA5 - 63.0 * ZETA7 / (4.0 * am**2)) / am
    return rhonu, pnu


def legacy_low_drho(am: float) -> float:
    am2 = am**2
    return am2 * (2.0 * CONST2 + am2 * (0.04399706676 * np.log(am) - 0.002970400378 - 0.029331377855 * am))


def legacy_high_drho(am: float) -> float:
    return (
        3.0 / (2.0 * FERMI_DIRAC_CONST) * (ZETA3 * am + (-(15.0 * ZETA5) / 2.0 + 2835.0 * ZETA7 / (16.0 * am**2)) / am)
    )


def build_legacy_spline_table() -> LegacySplineTable:
    am_grid = LEGACY_AM_MIN + np.arange(LEGACY_NRHOPN) * (LEGACY_AM_MAX - LEGACY_AM_MIN) / (LEGACY_NRHOPN - 1)
    rhonu = np.empty(LEGACY_NRHOPN)
    pnu = np.empty(LEGACY_NRHOPN)
    for index, am in enumerate(am_grid):
        rhonu[index], pnu[index] = legacy_nu_rho_pres(float(am))
    g = splini(LEGACY_NRHOPN)
    drho = splder(rhonu, g)
    dp = splder(pnu, g)
    return LegacySplineTable(am_grid, rhonu, pnu, drho, dp, splder(drho, g))


def legacy_spline_rho_p(table: LegacySplineTable, am: float) -> tuple[float, float]:
    if am <= LEGACY_AM_MINP:
        return legacy_low_rho_p(am)
    if am >= LEGACY_AM_MAXP:
        return legacy_high_rho_p(am)
    dam = table.am_grid[1] - table.am_grid[0]
    index = int((am - LEGACY_AM_MIN) / dam + 1.0) - 1
    fraction = (am - LEGACY_AM_MIN) / dam + 1.0 - (index + 1)

    def interpolate(values: np.ndarray, derivatives: np.ndarray) -> float:
        left = values[index]
        right = values[index + 1]
        dleft = derivatives[index]
        dright = derivatives[index + 1]
        return left + fraction * (
            dleft
            + fraction
            * (3.0 * (right - left) - 2.0 * dleft - dright + fraction * (dleft + dright + 2.0 * (left - right)))
        )

    return interpolate(table.rhonu, table.drho), interpolate(table.pnu, table.dp)


def legacy_spline_drho(table: LegacySplineTable, am: float) -> float:
    if am < LEGACY_AM_MINP:
        return legacy_low_drho(am)
    if am > LEGACY_AM_MAXP:
        return legacy_high_drho(am)
    dam = table.am_grid[1] - table.am_grid[0]
    index = int((am - LEGACY_AM_MIN) / dam + 1.0) - 1
    fraction = (am - LEGACY_AM_MIN) / dam + 1.0 - (index + 1)
    left = table.drho[index]
    right = table.drho[index + 1]
    dleft = table.ddrho[index]
    dright = table.ddrho[index + 1]
    derivative = left + fraction * (
        dleft
        + fraction * (3.0 * (right - left) - 2.0 * dleft - dright + fraction * (dleft + dright + 2.0 * (left - right)))
    )
    return am * derivative / dam


def summarize_error(values: np.ndarray, reference: np.ndarray, ams: np.ndarray) -> ErrorSummary:
    relative_error = np.abs(values - reference) / np.abs(reference)
    index = int(np.argmax(relative_error))
    return ErrorSummary(
        max_relative_error=float(relative_error[index]),
        rms_relative_error=float(np.sqrt(np.mean(np.square(relative_error)))),
        am_at_max_error=float(ams[index]),
    )


def print_error(label: str, summary: ErrorSummary) -> None:
    print(
        f"{label:<30} max = {summary.max_relative_error:.3e} at am = {summary.am_at_max_error:.6g}; "
        f"RMS = {summary.rms_relative_error:.3e}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate CAMB thermal-neutrino background fits")
    parser.add_argument("--samples", type=int, default=801, help="Number of logarithmically spaced samples per range")
    parser.add_argument("--am-min", type=float, default=1e-4, help="Minimum am for the direct-fit validation range")
    parser.add_argument("--am-max", type=float, default=1e4, help="Maximum am for the direct-fit validation range")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.samples < 2 or args.am_min <= 0.0 or args.am_max <= args.am_min:
        raise ValueError("Require at least two samples and 0 < am_min < am_max")

    fit_ams = np.geomspace(args.am_min, args.am_max, args.samples)
    fit_rho_p = np.array([thermal_nu_background_rho_p(float(am)) for am in fit_ams])
    fit_rho = np.array([thermal_nu_background_rho(float(am)) for am in fit_ams])
    reference_rho_p = np.array([exact_rho_p(float(am)) for am in fit_ams])
    fit_drho = np.array([thermal_nu_background_drho(float(am)) for am in fit_ams])
    reference_drho = np.array([exact_drho(float(am)) for am in fit_ams])
    if not np.array_equal(fit_rho, fit_rho_p[:, 0]):
        raise RuntimeError("ThermalNuBackground_rho no longer matches the rho_P density branch")

    print(f"Live Fortran fit reproduction over {args.am_min:g} <= am <= {args.am_max:g} ({args.samples} samples)")
    print_error("rho", summarize_error(fit_rho, reference_rho_p[:, 0], fit_ams))
    print_error("P", summarize_error(fit_rho_p[:, 1], reference_rho_p[:, 1], fit_ams))
    print_error("drho / adotoa", summarize_error(fit_drho, reference_drho, fit_ams))

    legacy_ams = np.geomspace(LEGACY_AM_MIN, LEGACY_AM_MAX, args.samples)
    legacy_table = build_legacy_spline_table()
    legacy_rho_p = np.array([legacy_spline_rho_p(legacy_table, float(am)) for am in legacy_ams])
    legacy_drho = np.array([legacy_spline_drho(legacy_table, float(am)) for am in legacy_ams])
    legacy_reference_rho_p = np.array([exact_rho_p(float(am)) for am in legacy_ams])
    legacy_reference_drho = np.array([exact_drho(float(am)) for am in legacy_ams])

    print(
        f"\nRemoved legacy spline baseline over {LEGACY_AM_MIN:g} <= am <= {LEGACY_AM_MAX:g} ({args.samples} samples)"
    )
    print_error("rho", summarize_error(legacy_rho_p[:, 0], legacy_reference_rho_p[:, 0], legacy_ams))
    print_error("P", summarize_error(legacy_rho_p[:, 1], legacy_reference_rho_p[:, 1], legacy_ams))
    print_error("drho / adotoa", summarize_error(legacy_drho, legacy_reference_drho, legacy_ams))


if __name__ == "__main__":
    main()
