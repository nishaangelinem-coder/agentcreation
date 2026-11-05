"""siwlab - open-source design environment for low-profile SIW slot arrays
with characteristic-mode-guided decoupling (28 GHz 5G mmWave).

Modules
-------
siw       : SIW design equations (effective width, cutoff, dispersion, loss)
slots     : Elliott/Stegen standing-wave slot-array synthesis
mesh      : structured triangular mesher for metal plates with apertures
mom       : RWG EFIE Method-of-Moments and characteristic-mode solver
fdtd      : Numba-accelerated 3-D FDTD (nonuniform Yee grid, CPML, SIW ports,
            aperture near-to-far-field transformation)
geometry  : parametric builders for the 28 GHz dual-port SIW slot array
"""
__version__ = "0.1.0"
