from .baseconfig import F2003Class, fortran_class, f_pointer
from ctypes import c_bool, c_double, POINTER, byref, c_void_p


class ReionizationModel(F2003Class):
    """
    Abstract base class for reionization models.
    """
    _fields_ = [
        ("Reionization", c_bool, "Is there reionization? (can be off for matter power, which is independent of it)")]


class BaseTauWithHeReionization(ReionizationModel):
    """
    Abstract class for models that map z_re to tau, and include second reionization of Helium
    """
    _fields_ = [
        ("use_optical_depth", c_bool, "Whether to use the optical depth or redshift paramters"),
        ("redshift", c_double, "Reionization redshift (xe=0.5) if use_optical_depth=False"),
        ("z_end", c_double, "Reionization endpoint if use_optical_depth=False"),
        ("dz", c_double, "Reionization duration (zre-zend)"),
        ("optical_depth", c_double, "Optical depth if use_optical_depth=True"),
        ("fraction", c_double,
         "Reionization fraction when complete, or -1 for full ionization of hydrogen and first ionization of helium."),
        ("include_helium_fullreion", c_bool, "Whether to include second reionization of helium"),
        ("helium_redshift", c_double, "Redshift for second reionization of helium"),
        ("helium_delta_redshift", c_double, "Width in redshift for second reionization of helium"),
        ("helium_redshiftstart", c_double, "Include second helium reionization below this redshift"),
        ("tau_solve_accuracy_boost", c_double, "Accuracy boosting parameter for solving for z_re from tau"),
        ("timestep_boost", c_double,
         "Accuracy boosting parameter for the minimum number of time sampling steps through reionization"),
        ("max_redshift", c_double, "Maxmimum redshift allowed when mapping tau into reionization redshift"),
        ("__min_redshift", c_double, "Minimim redshift allowed when mapping tau into reionization redshift"),
        ("__fHe", c_double, "Helium fraction"),
        ("__state", f_pointer)
    ]

    _fortran_class_module_ = 'Reionization'
    _fortran_class_name_ = 'TBaseTauWithHeReionization'

    _methods_ = [
        ('GetZreFromTau', [c_void_p, POINTER(c_double)], c_double, {"nopass": True}),
        ('GetTauFromXe', [c_void_p, POINTER(c_double)], c_double, {"nopass": True})
        ]
    
    def get_zre(self, params, tau=None):
        """
        Get the midpoint redshift of reionization.

        :param params: :class:`.model.CAMBparams` instance with cosmological parameters
        :param tau: if set, calculate the redshift for optical depth tau, otherwise uses curently set parameters
        :return: reionization mid-point redshift
        """
        if self.use_optical_depth or tau:
            from .camb import CAMBparams
            assert isinstance(params, CAMBparams)
            return self.f_GetZreFromTau(byref(params), c_double(tau or self.optical_depth))
        else:
            return self.redshift

    def set_zrei(self, zrei, dz=None):
        """
        Set the mid-point reionization redshift

        :param zrei: mid-point redshift
        :param dz:  reion duration
        :return:  self
        """
        self.use_optical_depth = False
        self.redshift = zrei
        if dz is not None:
            self.dz = dz
            self.z_end = self.redshift - self.dz
        return self

    def get_tau(self, params, zre=None, dz=None):
        """
        Get the reionization optical depth.

        :param params: :class:`.model.CAMBparams` instance with cosmological parameters
        :param zre: if set, calculate the optical depth for zre, otherwise uses curently set parameters
        :param dz: if set, calculate the optical depth for dz, otherwise uses curently set parameters
        :return: reionization optical depth
        """
        if self.use_optical_depth:
            return self.optical_depth
        else:
            from .camb import CAMBparams
            assert isinstance(params, CAMBparams)
            if dz is None:
                dz = self.dz
            if zre is None:
                zre = self.redshift
            param_reion = [zre, dz]
            param_reion_ptr = (c_double * 2)(*param_reion)
            return self.f_GetTauFromXe(byref(params), param_reion_ptr)

    def set_tau(self, tau):
        """
        Set the optical depth

        :param tau: optical depth
        :return: self
        """
        self.use_optical_depth = True
        self.optical_depth = tau
        return self


@fortran_class
class TanhReionization(BaseTauWithHeReionization):
    """
    This default (unphysical) tanh x_e parameterization is described in
    Appendix B of `arXiv:0804.3865 <https://arxiv.org/abs/0804.3865>`_
    """

    _fortran_class_name_ = 'TTanhReionization'

    def set_extra_params(self):
        
        self.z_end = self.redshift - self.dz


@fortran_class
class ExpReionization(BaseTauWithHeReionization):
    """
        An ionization fraction that decreases exponentially at high z, saturating to fully inionized at fixed redshift.
        This model has a minimum non-zero tau around 0.04 for reion_redshift_complete=6.1.
        Similar to e.g. arXiv:1509.02785, arXiv:2006.16828, but not attempting to fit shape near x_e~1 at z<6.1
    """
    _fields_ = [
        # ("reion_redshift_complete", c_double, "end of reionization"),
        ("reion_exp_smooth_width", c_double, "redshift scale to smooth exponential"),
        ("reion_exp_power", c_double, "power in exponential, exp(-lambda(z-redshift_complete)^exp_power)")]

    _fortran_class_name_ = 'TExpReionization'

    def set_extra_params(self, reion_exp_power=None, reion_exp_smooth_width=None):
        """
        Set extra parameters (not tau, or zrei)

        :param reion_redshift_complete: redshift at which reionization complete (e.g. around 6)
        :param reion_exp_power: power in exponential decay with redshift
        :param reion_exp_smooth_width: smoothing parameter to keep derivative smooth
        """
        self.z_end = self.redshift - self.dz
        if reion_exp_power is not None:
            self.reion_exp_power = reion_exp_power
        if reion_exp_smooth_width is not None:
            self.reion_exp_smooth_width = reion_exp_smooth_width
