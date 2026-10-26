"""AC grid model: measurement function h(x), Jacobian, Newton-Raphson power flow
and WLS state estimation built on pandapower's internal (ppc) representation."""
import warnings
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

warnings.filterwarnings("ignore")
import pandapower as pp
import pandapower.networks as pn
from pandapower.pypower.makeYbus import makeYbus
from pandapower.pypower.dSbus_dV import dSbus_dV
from pandapower.pypower.dSbr_dV import dSbr_dV

CASES = {"ieee14": pn.case14, "ieee118": pn.case118}


class Grid:
    def __init__(self, name):
        net = CASES[name]()
        pp.runpp(net, numba=False)
        ppc = net._ppc
        self.name = name
        self.base = ppc["baseMVA"]
        self.bus = ppc["bus"].copy()
        self.branch = ppc["branch"].copy()
        self.gen = ppc["gen"].copy()
        self.n = self.bus.shape[0]
        self.m = self.branch.shape[0]
        self.Ybus, self.Yf, self.Yt = makeYbus(self.base, self.bus, self.branch)
        self.Ybus, self.Yf, self.Yt = (sp.csr_matrix(M) for M in (self.Ybus, self.Yf, self.Yt))
        self.f = self.branch[:, 0].real.astype(int)
        self.t = self.branch[:, 1].real.astype(int)
        btype = self.bus[:, 1].astype(int)
        self.ref = int(np.where(btype == 3)[0][0])
        self.pv = np.where(btype == 2)[0]
        self.pq = np.where(btype == 1)[0]
        self.Pd0 = self.bus[:, 2].real / self.base
        self.Qd0 = self.bus[:, 3].real / self.base
        gbus = self.gen[:, 0].real.astype(int)
        self.gbus = gbus
        self.Pg0 = self.gen[:, 1].real / self.base
        self.Vset = np.ones(self.n)
        self.Vset[gbus] = self.gen[:, 5].real
        self.nonref = np.array([i for i in range(self.n) if i != self.ref])
        # adjacency (undirected, unique pairs)
        A = np.zeros((self.n, self.n))
        A[self.f, self.t] = 1
        A[self.t, self.f] = 1
        self.adj = A
        self.nz = 3 * self.n + 2 * self.m  # |V|, P_inj, Q_inj, P_f, Q_f

    # ---------------- measurement model ----------------
    def h(self, Va, Vm):
        V = Vm * np.exp(1j * Va)
        S = V * np.conj(self.Ybus @ V)
        Sf = V[self.f] * np.conj(self.Yf @ V)
        return np.concatenate([Vm, S.real, S.imag, Sf.real, Sf.imag])

    def H(self, Va, Vm):
        """Jacobian wrt state [Va(nonref), Vm(all)]."""
        V = Vm * np.exp(1j * Va)
        dS_dVm, dS_dVa = dSbus_dV(self.Ybus, V)
        dSf_dVa, dSf_dVm, _, _, _, _ = dSbr_dV(self.branch, self.Yf, self.Yt, V)
        n = self.n
        I = sp.eye(n, format="csr")
        Z = sp.csr_matrix((n, n))
        top = sp.hstack([Z, I])
        rows = [top,
                sp.hstack([dS_dVa.real, dS_dVm.real]),
                sp.hstack([dS_dVa.imag, dS_dVm.imag]),
                sp.hstack([dSf_dVa.real, dSf_dVm.real]),
                sp.hstack([dSf_dVa.imag, dSf_dVm.imag])]
        Hf = sp.vstack(rows).tocsc()
        cols = np.concatenate([self.nonref, n + np.arange(n)])
        return Hf[:, cols]

    def sigma(self, s_v=0.004, s_inj=0.01, s_flow=0.008):
        n, m = self.n, self.m
        return np.concatenate([np.full(n, s_v), np.full(2 * n, s_inj), np.full(2 * m, s_flow)])

    # ---------------- power flow ----------------
    def powerflow(self, Pd, Qd, Pg, tol=1e-8, maxit=20):
        """Polar Newton-Raphson. Pd,Qd per bus (pu); Pg per gen row (pu, ref ignored)."""
        n = self.n
        Sbus = -(Pd + 1j * Qd).astype(complex)
        np.add.at(Sbus, self.gbus, Pg)  # Q of PV/ref buses is free
        Va = np.zeros(n)
        Vm = np.ones(n)
        Vm[self.pv] = self.Vset[self.pv]
        Vm[self.ref] = self.Vset[self.ref]
        pvpq = np.r_[self.pv, self.pq]
        pq = self.pq
        for _ in range(maxit):
            V = Vm * np.exp(1j * Va)
            mis = V * np.conj(self.Ybus @ V) - Sbus
            F = np.r_[mis[pvpq].real, mis[pq].imag]
            if np.max(np.abs(F)) < tol:
                return Va, Vm, True
            dS_dVm, dS_dVa = dSbus_dV(self.Ybus, V)
            J = sp.vstack([sp.hstack([dS_dVa[pvpq][:, pvpq].real, dS_dVm[pvpq][:, pq].real]),
                           sp.hstack([dS_dVa[pq][:, pvpq].imag, dS_dVm[pq][:, pq].imag])]).tocsc()
            dx = spsolve(J, -F)
            Va[pvpq] += dx[: len(pvpq)]
            Vm[pq] += dx[len(pvpq):]
        return Va, Vm, False

    # ---------------- WLS state estimation ----------------
    def wls(self, z, sig, Va0=None, Vm0=None, tol=1e-6, maxit=15):
        n = self.n
        Va = np.zeros(n) if Va0 is None else Va0.copy()
        Vm = np.ones(n) if Vm0 is None else Vm0.copy()
        W = sp.diags(1.0 / sig ** 2)
        for it in range(maxit):
            r = z - self.h(Va, Vm)
            Hm = self.H(Va, Vm)
            G = (Hm.T @ W @ Hm).tocsc()
            dx = spsolve(G, Hm.T @ (W @ r))
            Va[self.nonref] += dx[: n - 1]
            Vm += dx[n - 1:]
            if np.max(np.abs(dx)) < tol:
                break
        r = z - self.h(Va, Vm)
        J = float(np.sum((r / sig) ** 2))
        return Va, Vm, J, r
