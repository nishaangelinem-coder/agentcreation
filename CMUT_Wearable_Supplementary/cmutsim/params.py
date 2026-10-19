"""Design parameters taken from the manuscript (Tables I-III) plus the
assumptions the manuscript leaves unstated (marked ASSUMED)."""
import numpy as np

# ---------------- Table I : array and CMUT design -----------------
F0 = 1.5e6            # operating frequency (Hz)
C0_TISSUE = 1540.0    # beamformer sound speed (m/s)
RHO0 = 1050.0         # reference tissue density (kg/m^3)
LAM = C0_TISSUE / F0  # wavelength (m)
N_ELEM = 64
PITCH = 0.50e-3
WIDTH = 0.45e-3
KERF = PITCH - WIDTH
APERTURE = N_ELEM * PITCH
FOCUS_Z = 20e-3
P_TARGET = 12e3       # focal peak pressure set-point (Pa)

# CMUT cell (silicon-nitride plate, vacuum gap)
CELL_RADIUS = 30e-6
PLATE_T = 1.5e-6
GAP = 0.20e-6
E_PLATE = 250e9       # ASSUMED Si3N4 Young's modulus (Pa)
NU_PLATE = 0.23       # ASSUMED Poisson ratio
RHO_PLATE = 3100.0    # ASSUMED Si3N4 density (kg/m^3)
T_INSULATOR = 0.15e-6 # ASSUMED SiO2 insulation layer under top electrode
EPS_INS = 3.9
EPS0 = 8.854e-12
BW_FRAC = 0.518
Q_LOADED = 3.5

# Element geometry not given in the manuscript
ELEV_HEIGHT = 6.0e-3  # ASSUMED element elevation height (m)
CELL_PITCH = 66e-6    # ASSUMED cell pitch (2a + 6 um anchor)

# ---------------- Table II : layered tissue model -----------------
# name, z_start (m), z_end (m), c (m/s), rho (kg/m^3), alpha (dB/cm/MHz), G (Pa)
LAYERS = [
    ("Skin",   0.0e-3,  2.0e-3, 1610.0, 1090.0, 0.60, 25.0e3),
    ("Fat",    2.0e-3, 10.0e-3, 1450.0,  950.0, 0.48,  2.0e3),
    ("Muscle",10.0e-3, 18.0e-3, 1580.0, 1050.0, 0.57, 12.0e3),
    ("Tumour",18.0e-3, 28.0e-3, 1560.0, 1050.0, 0.70, 30.0e3),
]
# below 28 mm assume tumour/muscle-like continuation
DEEP = ("Deep", 28.0e-3, 60e-3, 1570.0, 1050.0, 0.60, 12.0e3)

# Thermal properties (IT'IS database v4.1, representative values)
THERMAL = {  # k (W/m/K), rho*c (J/m^3/K), perfusion w_b (1/s)
    "Skin":   (0.37, 1109 * 3391, 1.77e-3),
    "Fat":    (0.21,  911 * 2348, 0.50e-3),
    "Muscle": (0.49, 1090 * 3421, 0.62e-3),
    "Tumour": (0.51, 1050 * 3600, 1.00e-3),
    "Deep":   (0.49, 1090 * 3421, 0.62e-3),
}
RHO_C_BLOOD = 1050 * 3617

DB_PER_NP = 20 / np.log(10)   # 8.686


def layer_at(z):
    """Return layer tuple at depth z (m)."""
    for L in LAYERS:
        if L[1] <= z < L[2]:
            return L
    return DEEP


def alpha_np(alpha_db_cm_mhz, f=F0):
    """Attenuation in Np/m at frequency f (linear frequency law)."""
    return alpha_db_cm_mhz * (f / 1e6) * 100.0 / DB_PER_NP
