"""Lightweight regression tests for siwlab (run: python -m pytest tests -q)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from siwlab import siw
from siwlab.fdtd import graded_axis, FDTD, gaussian_modulated, spectrum
from siwlab.mesh import plate_mesh, rwg_edges
from siwlab.mom import RWGModel, characteristic_modes, aperture_admittance_matrix
from siwlab.constants import c0


def test_siw_equations():
    a_eff = siw.effective_width(5.6e-3, 0.4e-3, 0.8e-3)
    assert abs(a_eff - 5.389e-3) < 1e-5
    fc = siw.cutoff_frequency(a_eff, 2.2)
    assert abs(fc / 1e9 - 18.75) < 0.05
    assert siw.beta(fc * 0.9, a_eff, 2.2) == 0.0


def test_graded_axis_monotonic():
    x = graded_axis([(0, 1, 0.1), (1, 2, 0.1, 0.3), (2, 3, 0.3)])
    assert np.all(np.diff(x) > 0) and abs(x[-1] - 3) < 1e-12


def test_mesh_and_rwg():
    nodes, tris = plate_mesh(-1, 1, -1, 1, [(-0.2, 0.2, -0.5, 0.5)], 0.25)
    r = rwg_edges(nodes, tris)
    assert r['nE'] > 0 and np.all(r['area'] > 0)


def test_strip_dipole_resonance():
    nodes, tris = plate_mesh(-5e-3, 5e-3, -0.25e-3, 0.25e-3, [], 0.5e-3)
    m = RWGModel(nodes, tris)
    lam = [characteristic_modes(m.impedance_matrix(f), 2)[0][0] for f in (13e9, 15e9)]
    assert lam[0] < 0 < lam[1]      # resonance between 13 and 15 GHz (L ~ 0.46 lambda)


def test_aperture_resonance_in_air():
    nodes, tris = plate_mesh(-0.2e-3, 0.2e-3, -2e-3, 2e-3, [], 0.2e-3)
    m = RWGModel(nodes, tris)
    lam = [characteristic_modes(aperture_admittance_matrix(m, f, (1.0, 1.0)), 2)[0][0] for f in (32e9, 37e9)]
    assert lam[0] * lam[1] < 0       # 4 mm slot resonates near 34.8 GHz (L ~ 0.46 lambda)


def test_fdtd_small_domain_stable():
    x = graded_axis([(0, 2e-3, 0.2e-3)]); y = graded_axis([(0, 3e-3, 0.2e-3)]); z = graded_axis([(0, 1e-3, 0.2e-3)])
    sim = FDTD(x, y, z, npml=(3, 3, 3, 3, 0, 3))
    sim.build()
    port = sim.add_siw_port('P', 4, 8, 1e-3, 1.5e-3, 0, 3)
    out = sim.run(300, port, gaussian_modulated(28e9, 6e9), [port])
    assert np.all(np.isfinite(out['rec']['P']['V']))
