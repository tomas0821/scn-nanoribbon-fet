"""Tunnelling guard (plan D6, heavy-plan M16; settings fits/PREREG_corr.md).

Effective-mass, single-band, OFF state of the L = 15 nm channel. Energies in eV with the source band edge at 0; the
carrier is treated as electron-like (holes mirrored), so one code serves n and p.

  * Barrier: the natural-length profile ψ'' = (ψ − U_g)/λ_s² on [0, L], ψ(0) = 0, ψ(L) = −V_DD, flat leads. U_g is
    set so that max ψ = δ_S + E_bF (δ_S = source Fermi level above its band edge; E_bF = OFF-state barrier above the
    source Fermi level, both from the ToB pin).
  * Source–drain tunnelling: R = I(T_qm) / I(T_classical) with Landauer occupations f_S − f_D. T_qm by a transfer
    matrix through the piecewise-constant profile (exact within the effective-mass model); WKB alongside.
  * Interband (GIDL-type) window on the drain side: electrons in the channel's opposite band under the barrier
    tunnel into empty drain states. Upper bound I ≤ (2q/h) ∫ T_K(E) (1 − f_D) dE, with T_K from the two-band (Kane)
    WKB decay κ = √(2 m_r (E − E_v)(E_c − E) / E_g)/ħ and the opposite band filled. It counts every generated
    carrier as leaving the channel, which needs a body contact or fast recombination.
"""
import numpy as np
from scipy.special import expit

H2M = 3.80998212            # ħ²/2m_e, eV Å²
G0 = 7.748091729e-5         # 2q²/h, A/V (spin included)
KB = 8.617333262e-5
_trapz = getattr(np, "trapezoid", None) or np.trapz     # numpy < 2 has only trapz


def profile(x, L, lam, Ug, VDD):
    s = np.sinh(L / lam)
    return Ug + (0.0 - Ug) * np.sinh((L - x) / lam) / s + (-VDD - Ug) * np.sinh(x / lam) / s


def barrier(L, lam, top, VDD, n=1500):
    """Segment midpoints (Å), width dx (Å) and ψ with max ψ = top (lengths in Å)."""
    x = (np.arange(n) + 0.5) * L / n
    lo, hi = -2.0, 5.0
    for _ in range(80):
        Ug = 0.5 * (lo + hi)
        if profile(x, L, lam, Ug, VDD).max() > top:
            hi = Ug
        else:
            lo = Ug
    return x, L / n, profile(x, L, lam, Ug, VDD)


def t_qm(E, U, dx, m, VDD):
    """Transmission through piecewise-constant U (left lead 0, right lead −VDD), vectorised over E."""
    E = np.asarray(E, float)
    k = lambda V: np.sqrt((E - V).astype(complex) * m / H2M)
    kL, kR = k(0.0), k(-VDD)
    M = np.zeros(E.shape + (2, 2), complex)
    M[..., 0, 0] = M[..., 1, 1] = 1.0
    for Uj in U:
        kj = k(Uj)
        kj = np.where(np.abs(kj) < 1e-9, 1e-9, kj)
        c, s = np.cos(kj * dx), np.sin(kj * dx)
        P = np.empty_like(M)
        P[..., 0, 0], P[..., 0, 1], P[..., 1, 0], P[..., 1, 1] = c, s / kj, -kj * s, c
        M = P @ M
    v1 = M @ np.stack([np.ones_like(kL), 1j * kL], -1)[..., None]
    v2 = M @ np.stack([np.ones_like(kL), -1j * kL], -1)[..., None]
    v1, v2 = v1[..., 0], v2[..., 0]
    r = -(v1[..., 1] - 1j * kR * v1[..., 0]) / (v2[..., 1] - 1j * kR * v2[..., 0])
    t = v1[..., 0] + r * v2[..., 0]
    T = np.real(kR) / np.real(kL) * np.abs(t) ** 2
    return np.where(np.real(kL) > 0, T, 0.0)


def t_wkb(E, U, dx, m):
    kap = np.sqrt(np.clip(U[None, :] - np.asarray(E)[:, None], 0, None) * m / H2M)
    return np.exp(-2 * kap.sum(axis=1) * dx)


def sdt(m, lam_nm, E_bF, delta_S, VDD=0.7, L_nm=15.0, T=300.0, nE=1200):
    """OFF-state source–drain tunnelling ratios (transfer matrix, WKB)."""
    kT = KB * T
    L, lam = 10 * L_nm, 10 * lam_nm
    top = delta_S + E_bF
    x, dx, U = barrier(L, lam, top, VDD)
    E = np.linspace(1e-4, top + 25 * kT, nE)
    occ = expit(-(E - delta_S) / kT) - expit(-(E - delta_S + VDD) / kT)
    Tc = (E > U.max()).astype(float)
    Ic = _trapz(Tc * occ, E)
    Tw = np.where(E > U.max(), 1.0, t_wkb(E, U, dx, abs(m)))
    # the transfer matrix loses precision in deep tunnelling (entries ~ e^{∫κ}); there T < e^-30 and the WKB value
    # is used instead (its share of the current is < 1e-13 of the thermionic part)
    Tq = np.where(Tw < np.exp(-30.0), Tw, t_qm(E, U, dx, abs(m), VDD))
    Rq = _trapz(Tq * occ, E) / Ic
    Rw = _trapz(Tw * occ, E) / Ic
    return Rq, Rw


def gidl_bound(m_r, Eg, lam_nm, E_bF, delta_S, VDD=0.7, L_nm=15.0, T=300.0, nE=800):
    """Upper bound on the drain-side interband current (A) and the window E_bF + VDD − Eg (eV)."""
    kT = KB * T
    L, lam = 10 * L_nm, 10 * lam_nm
    top = delta_S + E_bF
    x, dx, U = barrier(L, lam, top, VDD)
    window = E_bF + VDD - Eg
    Evmax = top - Eg
    if Evmax <= -VDD:
        return 0.0, window
    E = np.linspace(-VDD, Evmax, nE)
    i0 = int(np.argmax(U))
    Ud = U[i0:]                                  # drain side of the barrier
    Ec, Ev = Ud[None, :], Ud[None, :] - Eg
    arg = (E[:, None] - Ev) * (Ec - E[:, None]) / Eg
    kap = np.sqrt(np.clip(arg, 0, None) * m_r / H2M)
    Tk = np.exp(-2 * kap.sum(axis=1) * dx)
    empty = expit((E - (delta_S - VDD)) / kT)    # 1 − f_D
    return G0 * _trapz(Tk * empty, E), window
