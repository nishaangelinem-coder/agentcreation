"""Parametric geometry builder for SIW slot-array FDTD models (dimensions in mm)."""
import numpy as np
from dataclasses import dataclass, field
from .fdtd import FDTD, graded_axis
from .siw import effective_width
from .constants import eps0

MM = 1e-3


@dataclass
class SIWParams:
    er: float = 2.2
    tand: float = 0.0009
    h: float = 0.508          # substrate thickness (mm)
    d: float = 0.4            # via diameter (mm)
    p: float = 0.8            # via pitch (mm)
    a: float = 5.6            # via-row centre spacing (mm)
    f0: float = 28e9          # reference frequency for loss modelling

    @property
    def a_eff(self):
        return effective_width(self.a, self.d, self.p)

    @property
    def sigma_d(self):
        return 2 * np.pi * self.f0 * eps0 * self.er * self.tand


@dataclass
class Layout:
    """Everything in mm. One SIW per entry of xc. Slots: list of lists per SIW of
    (yc, xoff, L, w). Slits: list of (xc, yc, L, w) apertures in the common top wall."""
    xc: list = field(default_factory=lambda: [0.0])
    y_short: float = None        # y of shorting wall (None = through line into +y PML)
    slots: list = field(default_factory=list)
    slits: list = field(default_factory=list)
    walls_y: list = field(default_factory=list)   # extra PEC walls across the substrate: (y, x0, x1) in mm
    y_src: float = -26.0
    y_probe1: float = -22.0
    y_probe2: float = 10.0       # used for through-line / single slot characterisation
    y_min: float = -31.0
    y_max: float = 14.0
    x_margin: float = 3.5        # margin from outer via row to PML (fine->coarse)
    y_fine: tuple = (-21.0, 1.0) # fine region in y
    dfine: float = 0.1
    dcoarse: float = 0.2
    z_air: float = 6.0
    npml: int = 10


def make_grid(lay: Layout, prm: SIWParams):
    xs = np.asarray(lay.xc)
    xo = xs.max() + prm.a / 2 + 1.0
    xi = xs.min() - prm.a / 2 - 1.0
    pmlL = lay.npml * lay.dcoarse
    x = graded_axis([(xi - lay.x_margin - pmlL, xi - lay.x_margin, lay.dcoarse),
                     (xi - lay.x_margin, xi, lay.dcoarse, lay.dfine),
                     (xi, xo, lay.dfine),
                     (xo, xo + lay.x_margin, lay.dfine, lay.dcoarse),
                     (xo + lay.x_margin, xo + lay.x_margin + pmlL, lay.dcoarse)])
    yf0, yf1 = lay.y_fine
    y = graded_axis([(lay.y_min - pmlL, lay.y_min, lay.dcoarse),
                     (lay.y_min, yf0 - 1.0, lay.dcoarse),
                     (yf0 - 1.0, yf0, lay.dcoarse, lay.dfine),
                     (yf0, yf1, lay.dfine),
                     (yf1, yf1 + 1.0, lay.dfine, lay.dcoarse),
                     (yf1 + 1.0, lay.y_max, lay.dcoarse),
                     (lay.y_max, lay.y_max + pmlL, lay.dcoarse)])
    nz_sub = 5
    zpml = 0.25
    z = graded_axis([(0.0, prm.h, prm.h / nz_sub),
                     (prm.h, prm.h + lay.z_air, prm.h / nz_sub, zpml),
                     (prm.h + lay.z_air, prm.h + lay.z_air + lay.npml * zpml, zpml)])
    return x * MM, y * MM, z * MM, nz_sub


def build(lay: Layout, prm: SIWParams, verbose=True):
    x, y, z, nz_sub = make_grid(lay, prm)
    sim = FDTD(x, y, z, npml=(lay.npml, lay.npml, lay.npml, lay.npml, 0, lay.npml))
    kh = nz_sub  # node plane of the top metal
    # substrate everywhere below the top metal
    sim.set_box_material(-1, 1, -1, 1, -1, prm.h * MM * 0.999, prm.er, prm.sigma_d)
    # apertures in the top wall
    aps = []
    for s_list, xc in zip(lay.slots, lay.xc):
        for (yc, xoff, L, w) in s_list:
            aps.append(((xc + xoff - w / 2) * MM, (xc + xoff + w / 2) * MM, (yc - L / 2) * MM, (yc + L / 2) * MM))
    for (xc, yc, L, w) in lay.slits:
        aps.append(((xc - w / 2) * MM, (xc + w / 2) * MM, (yc - L / 2) * MM, (yc + L / 2) * MM))
    sim.add_pec_sheet_z(kh, aps)
    # via rows: from y_min (inside PML) to y_short (or y_max)
    y_end = lay.y_short if lay.y_short is not None else lay.y_max + lay.npml * lay.dcoarse
    y_start = lay.y_min - lay.npml * lay.dcoarse
    yv = np.arange(y_end, y_start, -prm.p)
    for xc in lay.xc:
        for xr in (xc - prm.a / 2, xc + prm.a / 2):
            for yc in yv:
                sim.add_pec_cylinder_z(xr * MM, yc * MM, prm.d / 2 * MM, 0, kh)
        if lay.y_short is not None:
            j = sim.idx('y', lay.y_short * MM)
            sim.add_pec_wall_y(j, (xc - prm.a / 2) * MM, (xc + prm.a / 2) * MM, 0, kh)
    for (yw, x0, x1) in lay.walls_y:
        sim.add_pec_wall_y(sim.idx('y', yw * MM), x0 * MM, x1 * MM, 0, kh)
    ports = []
    for n, xc in enumerate(lay.xc):
        jsrc = sim.idx('y', lay.y_src * MM); jpr = sim.idx('y', lay.y_probe1 * MM)
        ports.append(sim.add_siw_port(f'P{n+1}', jsrc, jpr, xc * MM, prm.a_eff * MM, 0, kh))
    if lay.y_short is None:
        for n, xc in enumerate(lay.xc):
            jpr = sim.idx('y', lay.y_probe2 * MM)
            ports.append(sim.add_siw_port(f'P{n+1}out', jsrc, jpr, xc * MM, prm.a_eff * MM, 0, kh))
    sim.build()
    # aperture DFT nodes (Ex, Ey on plane kh, non-PEC)
    ex = np.where(~sim.pecx[:, :, kh]); ey = np.where(~sim.pecy[:, :, kh])
    dft_nodes = dict(ex=(ex[0], ex[1], np.full(len(ex[0]), kh)), ey=(ey[0], ey[1], np.full(len(ey[0]), kh)))
    ap_geom = dict(ex_xy=(sim.xh[ex[0]], sim.y[ex[1]]), ex_dS=sim.dx[ex[0]] * sim.yd[ex[1]],
                   ey_xy=(sim.x[ey[0]], sim.yh[ey[1]]), ey_dS=sim.xd[ey[0]] * sim.dy[ey[1]])
    if verbose:
        print(f"grid {sim.Nx} x {sim.Ny} x {sim.Nz} = {sim.Nx*sim.Ny*sim.Nz/1e6:.2f} Mcells, dt = {sim.dt*1e12:.4f} ps, "
              f"aperture nodes: {len(ex[0])} Ex, {len(ey[0])} Ey")
    return sim, ports, dft_nodes, ap_geom
