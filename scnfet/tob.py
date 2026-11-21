"""Ballistic top-of-barrier (ToB) FET model after Rahman, Guo, Datta & Lundstrom, IEEE TED 50, 1853 (2003).

The conventions are pre-registered in 00_Master_Notebook.md ("Pre-registered device conventions"):
  * +k states are filled from mu_S and -k states from mu_D = mu_S - qV_DS. Time-reversal symmetry is assumed,
    so each half of the DOS is filled by one contact;
  * U_scf = -alpha_G V_GS - alpha_D V_DS + q (dN - dP) / C_Sigma, in eV, with ideal alpha_G = 1, alpha_D = 0 and
    C_Sigma = C_G;
  * C_G = eps0 * eps_ox * W / EOT per unit length, with eps_ox = 3.9 and W edge to edge;
  * transmission is mode counting, T(E) = M(E): conduction bands only for the n-FET, valence bands only for the
    p-FET;
  * mu_S comes from the S/D doping: 0.01 extra carriers per atom of the channel material, edge H not counted
    (signed off 2026-11-11; it reproduces Matić & Poljak's 2.7e13 cm^-2 for HfS2). I_OFF is pinned by a rigid
    gate work-function shift: I(V_GS = 0, V_DS = VDD) = I_OFF. Under this pin, with ideal electrostatics and a
    fixed C_G, the doping only translates I_D(V_GS), so I_ON does not depend on it.

Units: energies in eV, lengths in m internally, line densities in 1/m, currents in A. I/W is reported in mA/um
(1 A/m = 1e-3 mA/um).
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit

Q = 1.602176634e-19        # C
H_PLANCK = 6.62607015e-34  # J s
KB = 8.617333262e-5        # eV/K
EPS0 = 8.8541878128e-12    # F/m
A_PER_M_TO_MA_PER_UM = 1e-3


def fermi(E, mu, kT):
    return expit(-(E - mu) / kT)


@dataclass
class Channel:
    """Band structure of a 1D channel on a uniform periodic k-grid.

    E: (nk, nb) sorted eigenvalues (eV); a: transport period (m); n_val: number of valence bands in E.
    """
    E: np.ndarray
    a: float
    n_val: int
    spin: int = 2
    dE: float = 5e-4
    window: float = 2.0                  # eV above the CBM / below the VBM kept for the current integral
    _modes: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        self.Ev = self.E[:, :self.n_val].max()
        self.Ec = self.E[:, self.n_val:].min()
        if self.Ec <= self.Ev:
            raise ValueError(f"no gap between band {self.n_val} and {self.n_val + 1}: Ev={self.Ev}, Ec={self.Ec}")

    @property
    def Eg(self):
        return self.Ec - self.Ev

    def bands(self, carrier):
        return self.E[:, self.n_val:] if carrier == "n" else self.E[:, :self.n_val]

    def modes(self, carrier):
        """M(E) on an energy grid: half the number of crossings of level E over the periodic k-grid.

        This needs no band labels. Sorted bands cannot cross one level in opposite directions within one
        k-segment, so summing crossings over sorted bands counts the true crossings exactly.
        """
        if carrier in self._modes:
            return self._modes[carrier]
        if carrier == "n":
            grid = np.arange(self.Ec - 0.05, self.Ec + self.window, self.dE)
        else:
            grid = np.arange(self.Ev - self.window, self.Ev + 0.05, self.dE)
        B = self.bands(carrier)
        lo = np.minimum(B, np.roll(B, -1, axis=0)).ravel()   # each segment k_i -> k_{i+1}, periodic
        hi = np.maximum(B, np.roll(B, -1, axis=0)).ravel()
        diff = np.zeros(grid.size + 1)
        np.add.at(diff, np.searchsorted(grid, lo, side="right"), 1)
        np.add.at(diff, np.searchsorted(grid, hi, side="right"), -1)
        M = 0.5 * np.cumsum(diff)[:-1]
        self._modes[carrier] = (grid, M)
        return grid, M

    def density(self, U, muS, muD, kT, carrier):
        """Carriers per unit length at the top of the barrier (electrons in the CB, or holes in the VB)."""
        B = self.bands(carrier) + U
        if carrier == "n":
            occ = 0.5 * (fermi(B, muS, kT) + fermi(B, muD, kT))
        else:
            occ = 0.5 * (fermi(-B, -muS, kT) + fermi(-B, -muD, kT))     # 1 - f, written stably
        return self.spin * occ.sum() / (B.shape[0] * self.a)

    def _segments(self, carrier):
        key = "seg_" + carrier
        if key not in self._modes:
            B = self.bands(carrier)
            nxt = np.roll(B, -1, axis=0)
            self._modes[key] = (np.minimum(B, nxt).ravel(), np.maximum(B, nxt).ravel())
        return self._modes[key]

    def current(self, U, muS, muD, kT, carrier):
        """Landauer current (A): I = (spin q/h) ∫ M(E) [f_S - f_D] dE, signed (negative for holes S->D).

        M(E) = (1/2) * #{k-segments whose energy range contains E}. That makes the integral exact for linearly
        interpolated bands: (1/2) sum over segments of ∫_lo^hi (f_S - f_D) dE, done in closed form. Electrons use
        F(E) = kT ln(1 + e^{(mu-E)/kT}). Holes use the 1 - f antiderivative, kT ln(1 + e^{(E-mu)/kT}), so that
        deep-subthreshold hole currents are not differences of numbers close to 1.
        """
        lo, hi = self._segments(carrier)
        lo, hi = lo + U, hi + U
        if carrier == "n":
            F = lambda E, mu: kT * np.logaddexp(0.0, (mu - E) / kT)         # ∫_E^inf f
            integral = (F(lo, muS) - F(hi, muS)) - (F(lo, muD) - F(hi, muD))
        else:
            Hh = lambda E, mu: kT * np.logaddexp(0.0, (E - mu) / kT)        # ∫_-inf^E (1 - f)
            integral = (Hh(hi, muD) - Hh(lo, muD)) - (Hh(hi, muS) - Hh(lo, muS))
        return self.spin * Q / H_PLANCK * 0.5 * integral.sum() * Q          # eV -> J


@dataclass
class ToBFET:
    ch: Channel
    W: float                 # m, edge to edge
    EOT: float               # m
    carrier: str = "n"       # "n" or "p"
    n_atoms_per_cell: float = 1.0   # atoms of the channel material per transport cell (edge H not counted)
    doping: float = 0.01            # extra carriers per atom in the S/D
    T: float = 300.0
    eps_ox: float = 3.9
    alpha_G: float = 1.0
    alpha_D: float = 0.0
    C_G: float = None        # F/m; default parallel plate
    muS: float = None

    def __post_init__(self):
        self.kT = KB * self.T
        if self.C_G is None:
            self.C_G = EPS0 * self.eps_ox * self.W / self.EOT
        self.n_SD = self.doping * self.n_atoms_per_cell / self.ch.a
        if self.muS is None:
            self.muS = self._mu_from_doping()
        self.N0 = self.ch.density(0.0, self.muS, self.muS, self.kT, self.carrier)

    def _mu_from_doping(self):
        f = lambda mu: self.ch.density(0.0, mu, mu, self.kT, self.carrier) - self.n_SD
        lo, hi = self.ch.Ev - 3.0, self.ch.Ec + 3.0
        return brentq(f, lo, hi, xtol=1e-12)

    def solve(self, VGS, VDS):
        """Self-consistent U, carrier density and current at one bias point. Returns (U, n, I)."""
        muD = self.muS - VDS
        UL = -self.alpha_G * VGS - self.alpha_D * VDS
        sign = 1.0 if self.carrier == "n" else -1.0           # electrons raise U, holes lower it
        Cs = self.C_G

        def F(U):
            dn = self.ch.density(U, self.muS, muD, self.kT, self.carrier) - self.N0
            return U - UL - sign * Q * dn / Cs

        U = brentq(F, UL - 5.0, UL + 5.0, xtol=1e-10)
        n = self.ch.density(U, self.muS, muD, self.kT, self.carrier)
        I = self.ch.current(U, self.muS, muD, self.kT, self.carrier)
        return U, n, I

    def transfer(self, VGS, VDS, shift=0.0):
        """I_D(V_GS) at fixed V_DS (A). shift = gate work-function shift (V), applied as V_GS + shift."""
        return np.array([self.solve(v + shift, VDS)[2] for v in np.atleast_1d(VGS)])

    def metrics(self, VDD=0.7, I_OFF=10e-9 / 1e-6, n_sweep=141, span=1.5):
        """Pin I_OFF (A/m) at V_GS = 0 and read I_ON at V_GS = VDD. p-FET: biases are negative, |I| used.

        span: the gate sweep (±span V) that brackets the I_OFF point; it is a search range, not a convention. A
        doped ScN ribbon needs ≈ 2 V to deplete at EOT 1 nm (q·n_SD/C_G), so callers there pass a wider span with
        n_sweep scaled to keep the ≈ 21 mV step."""
        s = 1.0 if self.carrier == "n" else -1.0
        VD = s * VDD
        target = I_OFF * self.W
        vg = s * np.linspace(-span, span, n_sweep)
        I = np.abs(self.transfer(vg, VD))
        order = np.argsort(I)
        logI = np.log10(np.maximum(I[order], 1e-300))
        if not (logI[0] < np.log10(target) < logI[-1]):
            raise RuntimeError("I_OFF not bracketed by the gate sweep; widen it")
        v_off = np.interp(np.log10(target), logI, vg[order])  # gate voltage where |I| = I_OFF
        shift = v_off                                          # rigid work-function shift
        I_off = abs(self.solve(shift, VD)[2])
        U_on, n_on, I_on = self.solve(shift + s * VDD, VD)
        # subthreshold swing from the unpinned sweep: steepest decade/volt below I_OFF
        dlog = np.abs(np.gradient(np.log10(np.maximum(I, 1e-300)), vg))
        sub = (I < target) & (I > 1e-12 * target)
        SS = 1e3 / dlog[sub].max() if sub.any() else np.nan
        return dict(shift=shift, I_off=I_off / self.W * A_PER_M_TO_MA_PER_UM,
                    I_on=abs(I_on) / self.W * A_PER_M_TO_MA_PER_UM, SS_mV_dec=SS,
                    n_on=n_on, v_inj=abs(I_on) / (Q * n_on) if n_on > 0 else np.nan,
                    C_G=self.C_G, muS=self.muS, Eg=self.ch.Eg)
