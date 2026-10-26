"""Depth-resolved bacterial inactivation in the biofilm layer.

A Weibull (Mafart-type) dose-response is applied to the local lethal dose

    D(z) = int c_H2O2 dt  +  w_HNO2 * f_HNO2(pH) * int c_NO2- dt     [mM s]

so that ``log10 N/N0 (z) = -(D/delta)^p``. Membrane permeabilisation (the quantity
seen by impedance) is allowed to precede loss of culturability: the permeabilised
fraction follows the same law with a sensitivity multiplier ``kappa >= 1``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .plasma_chemistry import hno2_fraction


@dataclass
class WeibullInactivation:
    delta: float = 60.0       # scale dose for the first decimal reduction, mM s
    p: float = 0.85           # Weibull shape (<1: tailing, typical of biofilms)
    w_hno2: float = 60.0      # relative lethality of HNO2 vs H2O2
    kappa: float = 1.0        # >1: permeabilisation precedes loss of culturability
    desiccation_lr: float = 0.04   # log reduction caused by the gas flow alone (300 s)

    def local_dose(self, exp_h2o2, exp_no2, pH):
        return exp_h2o2 + self.w_hno2 * hno2_fraction(pH) * exp_no2

    def log_reduction(self, dose):
        return (np.clip(dose, 0, None) / self.delta) ** self.p


def biofilm_response(delivery, pH_profile, model: WeibullInactivation,
                     t_gas: float = 0.0, log10_N0: float = 8.0):
    """Summarise the biofilm state after treatment.

    Returns layer-resolved viable/permeabilised fractions and the
    volume-averaged log10 CFU per cm^2.
    """
    m = delivery["in_bio"]
    dz = delivery["dz"][m]
    dose = model.local_dose(delivery["exposure"]["H2O2"][m],
                            delivery["exposure"]["NO2"][m], pH_profile[m])
    lr = model.log_reduction(dose) + model.desiccation_lr * (t_gas / 300.0)
    viable = 10 ** (-lr)
    perm = 1.0 - 10 ** (-model.kappa * lr)
    viable_mean = np.average(viable, weights=dz)
    perm_mean = np.average(perm, weights=dz)
    return {
        "dose": dose,
        "viable_profile": viable,
        "perm_profile": perm,
        "viable_fraction": viable_mean,
        "permeabilised_fraction": perm_mean,
        "log_reduction": -np.log10(viable_mean),
        "log10_cfu": log10_N0 + np.log10(viable_mean),
    }
