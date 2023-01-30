from .baseconfig import F2003Class, fortran_class
from ctypes import c_bool, c_double, POINTER, byref, c_void_p


class ReionizationModel(F2003Class):
    """
    Abstract base class for reionization models.
    """
    _fields_ = [
        ("Reionization", c_bool, "Is there reionization? (can be off for matter power which is independent of it)")]


@fortran_class
class TanhReionization(ReionizationModel):
    """
    This default (unphysical) tanh x_e parameterization is described in
    Appendix B of `arXiv:0804.3865 <https://arxiv.org/abs/0804.3865>`_
    """
    _fields_ = [
        ("use_optical_depth", c_bool, "Whether to use the optical depth or redshift paramters"),
        ("redshift", c_double, "Reionization redshift if use_optical_depth-False"),
        ("optical_depth", c_double, "Optical depth if use_optical_depth=True"),
        ("fraction", c_double, "Reionization fraction when complete, or -1 for full ionization of hydrogen and first ionization of helium."),
        ("include_helium_fullreion", c_bool, "Whether to include second reionization of helium"),
        ("helium_redshift", c_double, "Redshift for second reionization of helium"),
        ("helium_delta_redshift", c_double, "Width in redshift for second reionization of helium"),
        ("helium_redshiftstart", c_double, "Include second helium reionizatio below this redshift"),
        ("tau_solve_accuracy_boost", c_double, "Accuracy boosting parameter for solving for z_re from tau"),
        ("timestep_boost", c_double, "Accuracy boosting parameter for the minimum number of time sampling steps through reionization"),
        ("z_end", c_double, "Reionization endpoint if use_optical_depth=False"),
        ("z_early", c_double, "Maxmimum redshift allowed when mapping tau into reionization redshift")]

    _fortran_class_module_ = 'Reionization'
    _fortran_class_name_ = 'TTanhReionization'

    _methods_ = [
        ('GetZreFromTau', [c_void_p, POINTER(c_double)], c_double, {"nopass": True}),
        ('GetTauFromXe', [c_void_p, POINTER(c_double)], c_double, {"nopass": True})
        ]

    def set_zrei(self, zrei, zend=None):
        """
        Set the mid-point reionization redshift

        :param zrei: mid-point redshift
        :param zend:  end-point redshif
        :return:  self
        """
        self.use_optical_depth = False
        self.redshift = zrei
        if z_end is not None:
            self.z_end = zend
        return self

    def set_tau(self, tau):
        """
        Set the optical depth

        :param tau: optical depth
        :return: self
        """
        self.use_optical_depth = True
        self.optical_depth = tau
        return self

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

    def get_tau(self, params, zre=None, zend=None):
        """
        Get the reionization optical depth.

        :param params: :class:`.model.CAMBparams` instance with cosmological parameters
        :param zre: if set, calculate the optical depth for zre, otherwise uses curently set parameters
        :param zend: if set, calculate the optical depth for zend, otherwise uses curently set parameters
        :return: reionization optical depth
        """
        if self.use_optical_depth:
            return self.optical_depth
        else:
            from .camb import CAMBparams
            assert isinstance(params, CAMBparams)
            return self.f_GetTauFromXe(byref(params), c_double(zre or self.redshift), c_double(zend or self.z_end))

