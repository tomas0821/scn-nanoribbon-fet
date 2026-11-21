"""Polar-optical (Fröhlich) mean free path and the ballisticity B (plan D6; settings fits/PREREG_corr.md).

Two couplings, both from the measured bulk ScN constants (Grümbel et al., PRMaterials 8, L071601, 2024):

  * "2d": Sohier, Calandra & Mauri, PRB 94, 085415 (2016), isotropic model (their Eqs. 8, 12, 13), used for the
    two-plane sheet: |g(q)|² A = K_A / (eps12 + r_eff q)², r_eff = eps_inf t / 2. K_A is fixed by the large-q limit,
    where the 2D coupling must equal the bulk one (their remark after Eq. 13):
        K_A = (e² ħω_LO / 2 ε_vac) (1/eps_inf − 1/eps_0) eps_inf² t / 4.
    The out-scattering rate for a parabolic band is their Eq. 14. With the final states on the circle |k'|,
        1/τ±(E) = (m N± K_A / π ħ³) ∫_0^π dθ / (eps12 + r_eff |k' − k|)²,  ħ²k'²/2m = E ± ħω_LO.
  * "3d": the bulk rate for a parabolic band (Lundstrom, Fundamentals of Carrier Transport):
        1/τ(E) = 2 α ω_LO √(ħω/E) [N0 asinh √(E/ħω) + (N0 + 1) asinh √(E/ħω − 1)],
        α = (1/eps_inf − 1/eps_0) √(m Ry / ħω).

Mean free path for backscattering in the 1D channel: every event is treated as randomising (conservative for the
forward-peaked Fröhlich coupling), so half reverse the velocity: λ(E) = 2 v(E) τ(E). The low-energy (kT-layer)
value is the 1D Landauer average with one mode, λ0 = ∫ λ(E) e^{−E/kT} dE / kT, and B = λ0 / (λ0 + L).

Units: eV, Å, s; λ returned in nm.
"""
import numpy as np
from scipy.integrate import quad

HBAR = 6.582119569e-16      # eV s
H2M = 3.80998212            # ħ²/2m_e, eV Å²
E2 = 14.39964548            # e²/4πε_vac, eV Å
RY = 13.605693122994        # eV
KB = 8.617333262e-5         # eV/K

# measured bulk ScN constants (Grümbel 2024, Table 2)
EPS_INF, EPS_0, HW_LO = 8.3, 29.5, 0.0849
T_SHEET = 4.50              # Å: two planes at the bulk (001) spacing a/2 = 2.25 Å


class Frohlich:
    def __init__(self, m, eps_inf=EPS_INF, eps_0=EPS_0, hw=HW_LO, t=T_SHEET, eps12=1.0, T=300.0):
        self.m, self.eps_inf, self.eps_0, self.hw, self.t, self.eps12 = abs(m), eps_inf, eps_0, hw, t, eps12
        self.kT = KB * T
        self.N0 = 1.0 / np.expm1(hw / self.kT)
        self.r = eps_inf * t / 2.0
        self.K_A = (4 * np.pi * E2 * hw / 2.0) * (1 / eps_inf - 1 / eps_0) * eps_inf ** 2 * t / 4.0
        self.alpha = (1 / eps_inf - 1 / eps_0) * np.sqrt(self.m * RY / hw)

    def k(self, E):
        return np.sqrt(E * self.m / H2M)

    def v(self, E):
        """Group velocity of the parabolic band, m/s."""
        return 2 * H2M * self.k(E) / self.m / HBAR * 1e-10

    def rate_2d(self, E):
        k, tot = self.k(E), 0.0
        for D, N in ((self.hw, self.N0), (-self.hw, self.N0 + 1)):
            if E + D <= 0:
                continue
            kp = np.sqrt(k * k + D * self.m / H2M)
            f = lambda th: 1.0 / (self.eps12 + self.r * np.sqrt(k * k + kp * kp - 2 * k * kp * np.cos(th))) ** 2
            tot += N * quad(f, 0, np.pi, limit=200)[0]
        return (self.m / (2 * H2M)) * self.K_A * tot / (np.pi * HBAR)

    def rate_3d(self, E):
        x = E / self.hw
        s = self.N0 * np.arcsinh(np.sqrt(x)) + ((self.N0 + 1) * np.arcsinh(np.sqrt(x - 1)) if x > 1 else 0.0)
        return 2 * self.alpha * (self.hw / HBAR) * np.sqrt(1 / x) * s

    def lam(self, E, model):
        rate = self.rate_2d(E) if model == "2d" else self.rate_3d(E)
        return 2 * self.v(E) / rate * 1e9

    def lam0(self, model):
        kT = self.kT
        f = lambda E: self.lam(E, model) * np.exp(-E / kT) / kT
        return quad(f, 1e-7, 40 * kT, limit=400, points=[self.hw])[0]


def ballisticity(lam0_nm, L_nm=15.0):
    return lam0_nm / (lam0_nm + L_nm)
