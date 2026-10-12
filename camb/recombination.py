from ctypes import c_bool, c_double, c_int

from .baseconfig import CAMBValueError, F2003Class, fortran_class, optional_fortran_class

recfast_planck = "planck"
recfast_cosmorec = "cosmorec"
recfast_hyrec = "hyrec"

recfast_default = recfast_cosmorec

# ##################################################################
# ######### feature added for clumping / PMF recombination #########
# ########## Three-zone phenomenological model of PMF-modified
# ########## recombination (see "Fresh Branch Tests/6 - Clumpy
# ########## Recombination.ipynb"). Enable with Recfast.clumping and
# ########## set the clumping factor Recfast.clumping_b >= 0.
# ##################################################################
clumping_model_M1 = 1  # low-density zone overdensity Delta_1 = 0.1
clumping_model_M2 = 2  # low-density zone overdensity Delta_1 = 0.3
# ##################################################################
# ################ end of feature #################################
# ##################################################################

recfast_approx_model_params = {
    recfast_planck: {
        "RECFAST_fudge": 1.125,
        "RECFAST_fudge_He": 0.86,
        "RECFAST_Hswitch": True,
        "RECFAST_He_rate_correction": False,
        "AGauss1": -0.14,
        "AGauss2": 0.079,
        "zGauss1": 7.28,
        "zGauss2": 6.73,
        "wGauss1": 0.18,
        "wGauss2": 0.33,
    },
    recfast_cosmorec: {
        "RECFAST_fudge": 1.125,
        "RECFAST_fudge_He": 0.8472367977,
        "RECFAST_Hswitch": True,
        "RECFAST_He_rate_correction": True,
        "AGauss1": -0.1395272483,
        "AGauss2": 0.0729891952,
        "zGauss1": 7.2813061282,
        "zGauss2": 6.7667038679,
        "wGauss1": 0.1638966410,
        "wGauss2": 0.2785834127,
    },
    recfast_hyrec: {
        "RECFAST_fudge": 1.125,
        "RECFAST_fudge_He": 0.8658491370,
        "RECFAST_Hswitch": True,
        "RECFAST_He_rate_correction": False,
        "AGauss1": -0.1417072428,
        "AGauss2": 0.0733978462,
        "zGauss1": 7.2841809376,
        "zGauss2": 6.7316088397,
        "wGauss1": 0.1734480042,
        "wGauss2": 0.3373321220,
    },
}


class RecombinationModel(F2003Class):
    """
    Abstract base class for recombination models
    """

    _fields_ = (
        (
            "min_a_evolve_Tm",
            c_double,
            "minimum scale factor at which to solve matter temperature "
            "perturbation if evolving sound speed or ionization fraction perturbations",
        ),
    )

    def write_ini(self, state) -> None:
        state.set("recombination_model", self.__class__.__name__)

    def set_params(self, **kwargs):
        pass


@fortran_class
class Recfast(RecombinationModel):
    """
    RECFAST recombination model (see recfast source for details).

    """

    _fields_ = (
        ("RECFAST_fudge", c_double),
        ("RECFAST_fudge_He", c_double),
        ("RECFAST_Hswitch", c_bool),
        ("RECFAST_He_rate_correction", c_bool),
        ("AGauss1", c_double),
        ("AGauss2", c_double),
        ("zGauss1", c_double),
        ("zGauss2", c_double),
        ("wGauss1", c_double),
        ("wGauss2", c_double),
        ("Nz", c_int),
        ("use_rosenbrock", c_bool),
        ("rosenbrock_handoff_xH", c_double),
        ("rosenbrock_tol", c_double),
        # ##### feature added for clumping / PMF recombination #####
        # ##### must match the field order of type TRecfast in fortran/recfast.f90 #####
        ("clumping", c_bool, "enable the three-zone PMF-clumping recombination model"),
        ("clumping_model", c_int, "clumping zone model: 1 = M1 (Delta_1=0.1), 2 = M2 (Delta_1=0.3)"),
        ("clumping_b", c_double, "clumping factor b = <Delta^2> - 1 (>= 0; 0 reproduces standard recombination)"),
        # ##### end of feature #####
    )

    _fortran_class_module_ = "Recombination"
    _fortran_class_name_ = "TRecfast"

    def set_params(self, recfast_approx_model=None, clumping=None, clumping_b=None, clumping_model=None):
        """
        Set RECFAST approximation parameters.

        :param recfast_approx_model: optional named approximation parameter set. One of

            - ``planck``: CAMB v1.x/RECFAST 1.5.2 Planck-era fudge parameters.
            - ``cosmorec``: RECFAST fit with an He rate correction to direct CosmoRec ``accuracy=6`` with
              10 H shells.
            - ``hyrec``: seven-parameter RECFAST fit to direct HyRec-2.
        :param clumping: (feature added for clumping) enable the three-zone PMF-clumping recombination model
        :param clumping_b: (feature added for clumping) clumping factor b = <Delta^2> - 1, must be >= 0
        :param clumping_model: (feature added for clumping) 1 = M1 (Delta_1=0.1) or 2 = M2 (Delta_1=0.3)

        :return: self
        """
        if recfast_approx_model is not None:
            try:
                params = recfast_approx_model_params[recfast_approx_model]
            except KeyError as err:
                allowed = ", ".join(sorted(recfast_approx_model_params))
                raise CAMBValueError(
                    f"Unknown RECFAST recfast_approx_model '{recfast_approx_model}'. Use one of: {allowed}"
                ) from err
            for name, value in params.items():
                setattr(self, name, value)

        # ##### feature added for clumping / PMF recombination #####
        if clumping_model is not None:
            self.clumping_model = int(clumping_model)
        if clumping_b is not None:
            if clumping_b < 0:
                raise CAMBValueError("clumping factor b must be non-negative")
            self.clumping_b = float(clumping_b)
        if clumping is not None:
            self.clumping = bool(clumping)
        # ##### end of feature #####

        return self

    def write_ini(self, state) -> None:
        super().write_ini(state)
        state.set("RECFAST_H_fudge", self.RECFAST_fudge)
        state.write_fields(
            self,
            names=(
                "RECFAST_fudge_He",
                "RECFAST_Hswitch",
                "RECFAST_He_rate_correction",
                "AGauss1",
                "AGauss2",
                "zGauss1",
                "zGauss2",
                "wGauss1",
                "wGauss2",
                "Nz",
                "use_rosenbrock",
                "rosenbrock_handoff_xH",
                "rosenbrock_tol",
                # ##### feature added for clumping / PMF recombination: ini write side, mirrors TRecfast_ReadParams #####
                "clumping",
                "clumping_model",
                "clumping_b",
                # ##### end of feature #####
            ),
            rename={
                "Nz": "RECFAST_nz",
                "use_rosenbrock": "RECFAST_use_rosenbrock",
                "rosenbrock_handoff_xH": "RECFAST_rosenbrock_handoff_xH",
                "rosenbrock_tol": "RECFAST_rosenbrock_tol",
            },
        )


@optional_fortran_class
class CosmoRec(RecombinationModel):
    """
    `CosmoRec <https://www.jb.man.ac.uk/~jchluba/Science/CosmoRec/CosmoRec.html>`_ recombination model.
    To use this, CAMB must be built with CosmoRec installed from `GitHub <https://github.com/cmbant/CosmoRec>`_
    and RECOMBINATION_FILES set to include cosmorec in the CAMB Makefile.

    """

    _fortran_class_module_ = "CosmoRec"
    _fortran_class_name_ = "TCosmoRec"

    _fields_ = (
        (
            "runmode",
            c_int,
            "Default 0, with diffusion; 1: without diffusion; 2: RECFAST++, 3: RECFAST++ run with correction",
        ),
        ("fdm", c_double, "Dark matter annihilation efficiency"),
        ("accuracy", c_double, "CosmoRec accuracy preset; see CosmoRec source for supported values"),
        ("diff_iteration_max", c_int, "Override CosmoRec diffusion iteration count; negative uses batch default (2)"),
        ("n_shells", c_int, "Override CosmoRec hydrogen shells; negative uses selected accuracy preset"),
        ("n_shells_hei", c_int, "Override CosmoRec helium shells; negative uses selected accuracy preset"),
        ("flag_hi_absorption", c_int, "Override CosmoRec HI absorption flag; negative uses selected accuracy preset"),
        ("n_s_2gamma", c_int, "Override CosmoRec two-photon shell limit; negative uses selected accuracy preset"),
        ("n_s_raman", c_int, "Override CosmoRec Raman shell limit; negative uses selected accuracy preset"),
    )

    def write_ini(self, state) -> None:
        super().write_ini(state)
        state.write_fields(
            self,
            names=(
                "runmode",
                "accuracy",
                "fdm",
                "diff_iteration_max",
                "n_shells",
                "n_shells_hei",
                "flag_hi_absorption",
                "n_s_2gamma",
                "n_s_raman",
            ),
            rename={
                "runmode": "cosmorec_runmode",
                "accuracy": "cosmorec_accuracy",
                "fdm": "cosmorec_fdm",
                "diff_iteration_max": "cosmorec_diff_iteration_max",
                "n_shells": "cosmorec_n_shells",
                "n_shells_hei": "cosmorec_n_shells_hei",
                "flag_hi_absorption": "cosmorec_flag_hi_absorption",
                "n_s_2gamma": "cosmorec_n_s_2gamma",
                "n_s_raman": "cosmorec_n_s_raman",
            },
        )


@optional_fortran_class
class HyRec(RecombinationModel):
    r"""
    `HyRec <https://github.com/nanoomlee/HYREC-2>`_ recombination model.
    To use this, the library must be built with HyRec installed and RECOMBINATION_FILES including hyrec in the Makefile.

    """

    _fortran_class_module_ = "HyRec"
    _fortran_class_name_ = "THyRec"
