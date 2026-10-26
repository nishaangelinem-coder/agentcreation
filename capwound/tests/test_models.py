import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from capwound import electrochem as ec                         # noqa: E402
from capwound.dielectric import Bacterium, biofilm_permittivity, debye_medium  # noqa: E402
from capwound.fitting import dielectric_model, fit_dielectric, fit_randles, lin_kk, randles_model  # noqa: E402
from capwound.impedance import IDESensor, frequency_grid      # noqa: E402
from capwound.plasma_chemistry import Column, PlasmaSource, buffered_ph, simulate_delivery  # noqa: E402

F = frequency_grid()
W = 2 * np.pi * F


def test_mass_conservation_without_sink():
    col = Column(k_cat=0.0)
    d = simulate_delivery(60, PlasmaSource(), col, t_hold=0, dt=0.5)
    total = np.sum(d["conc"]["H2O2"] * d["dz"])
    assert np.isclose(total, PlasmaSource().flux_h2o2 * 60, rtol=1e-6)


def test_buffer_ph_monotonic_and_neutral_start():
    ph = buffered_ph([0.0, 1.0, 5.0, 20.0])
    assert abs(ph[0] - 7.4) < 1e-6 and np.all(np.diff(ph) < 0)


def test_single_shell_limits():
    b = Bacterium()
    e = b.permittivity(np.array([2 * np.pi * 1e9]))
    # far above the membrane relaxation the shell is short-circuited: sigma -> close to cytoplasm
    sig = -e.imag[0] * 2 * np.pi * 1e9 * 8.854e-12
    assert 0.7 * b.sigma_cyt < sig < b.sigma_cyt
    # far below it the intact membrane insulates the cell
    e0 = b.permittivity(np.array([2 * np.pi * 1e2]))
    assert -e0.imag[0] * 2 * np.pi * 1e2 * 8.854e-12 < 1e-3


def test_permeabilisation_reduces_dispersion():
    eps = [biofilm_permittivity(W, 0.85, 0.3, p, Bacterium()).real[100] for p in (0.0, 1.0)]
    assert eps[0] > eps[1]


def test_dielectric_fit_recovers_parameters():
    p = (4e-7, 0.82, 120.0, 180.0, 1.6e7, 0.05, 0.6)
    z = dielectric_model(W, 600.0, *p)
    r = fit_dielectric(F, z, 600.0)
    for k, v in zip(("Q_ep", "n_ep", "eps_inf", "d_eps", "f_c", "alpha", "sigma"), p):
        assert np.isclose(r[k], v, rtol=2e-3), k


def test_randles_fit_recovers_rct():
    f = np.logspace(-1, 5, 61)
    z = randles_model(2 * np.pi * f, 180, 2e-6, 0.9, 1500, 400)
    assert np.isclose(fit_randles(f, z)["Rct"], 1500, rtol=1e-3)


def test_linkk_passes_on_causal_data():
    z = IDESensor().impedance(W, biofilm_permittivity(W, 0.85, 0.3, 0.2, Bacterium()), debye_medium(W, 1.0), 100e-6)
    M, mu, rr, ri = lin_kk(F, z)
    assert np.sqrt(np.mean(rr ** 2 + ri ** 2)) < 2e-3


def test_cv_reversible_limit():
    e = ec.SPCE(k0=1.0, Ru=0.0, C_dl_area=0.0)
    E, i = ec.cyclic_voltammogram(e, 0.0, 0.5)
    m = ec.cv_metrics(E, i)
    assert 0.055 < m["dEp"] < 0.070          # ~59 mV Nernstian separation
    # Randles-Sevcik (25 degC, C in mol/cm^3): starting at -0.2 V converts the probe to the
    # reduced form near the electrode, so ipa lies between the 2.5 mM and 5 mM limits
    rs = lambda c: 2.69e5 * e.area * 1e4 * np.sqrt(e.D_R * 1e4 * 0.05) * c * 1e-6
    assert rs(2.5) < m["ipa"] < rs(5.0)
