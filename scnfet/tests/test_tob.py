"""Analytic checks of scnfet.w90 and scnfet.tob (no external data).

Run: ~/.venvs/scnfet/bin/python -m pytest scnfet/tests -q   (from the project root)
"""
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import expit

from scnfet import tob, w90

HBAR = 1.054571817e-34
ME = 9.1093837015e-31


def chain_hr(e0=0.0, t=-2.0):
    """One-orbital chain: E(k) = e0 + 2t cos(2*pi*k)."""
    R = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
    H = np.array([[[e0]], [[t]], [[t]]], complex)
    return w90.HR(R=R, deg=np.ones(3), H=H)


def two_band_chain(gap=1.0, t=-1.5):
    """Valence and conduction cosine bands, uncoupled: Ev(k) = -gap/2 - 4|t| + 2|t|(1 - cos) mirrored."""
    R = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
    H = np.zeros((3, 2, 2), complex)
    H[0] = np.diag([-gap / 2 - 2 * abs(t), gap / 2 + 2 * abs(t)])
    H[1] = H[2] = np.diag([abs(t), -abs(t)])           # VB: max at k=0; CB: min at k=0
    return w90.HR(R=R, deg=np.ones(3), H=H)


def test_hr_roundtrip_and_bands(tmp_path):
    hr = chain_hr(0.3, -2.0)
    p = tmp_path / "x_hr.dat"
    w90.write_hr(p, hr)
    hr2 = w90.read_hr(p)
    k, E = w90.bands_1d(hr2, nk=400)
    assert np.allclose(E[:, 0], 0.3 - 4.0 * np.cos(2 * np.pi * k), atol=1e-9)


def test_modes_count():
    k, E = w90.bands_1d(two_band_chain(gap=1.0, t=-1.5), nk=2000)
    ch = tob.Channel(E=E, a=3e-10, n_val=1)
    assert ch.Eg == pytest.approx(1.0, abs=1e-6)
    grid, M = ch.modes("n")
    inside = (grid > ch.Ec + 0.01) & (grid < ch.Ec + 5.9)
    assert np.all(M[inside] == 1) and np.all(M[grid < ch.Ec - 0.001] == 0)
    gv, Mv = ch.modes("p")
    assert np.all(Mv[(gv < ch.Ev - 0.01) & (gv > ch.Ev - 1.9)] == 1)


@pytest.mark.parametrize("eta_S, eta_D", [(-5.0, -15.0), (0.0, -10.0), (3.0, 1.0)])
def test_current_matches_1d_landauer(eta_S, eta_D):
    """One mode above Ec: I = (2q/h) kT [ln(1+e^etaS) - ln(1+e^etaD)]."""
    k, E = w90.bands_1d(chain_hr(0.0, -2.0), nk=4000)       # band from -4 to +4 eV
    ch = tob.Channel(E=np.hstack([E - 100.0, E]), a=3e-10, n_val=1, window=6.0)  # dummy VB far below
    kT = tob.KB * 300
    muS, muD = ch.Ec + eta_S * kT, ch.Ec + eta_D * kT
    I = ch.current(0.0, muS, muD, kT, "n")
    I_ref = 2 * tob.Q / tob.H_PLANCK * kT * tob.Q * (np.log1p(np.exp(eta_S)) - np.log1p(np.exp(eta_D)))
    assert I == pytest.approx(I_ref, rel=1e-9)            # segment integration is exact


@pytest.mark.parametrize("eta", [-6.0, -2.0, 0.0, 2.0])
def test_density_matches_parabolic_F_minus_half(eta):
    """Near the band bottom a wide cosine band is parabolic: n = sqrt(2 m kT / pi) / hbar * F_{-1/2}(eta)."""
    t, a = 3.0, 3e-10
    k, E = w90.bands_1d(chain_hr(0.0, -t), nk=20000)
    ch = tob.Channel(E=np.hstack([E - 100.0, E]), a=a, n_val=1)
    kT = tob.KB * 300
    mu = ch.Ec + eta * kT
    n = ch.density(0.0, mu, mu, kT, "n")
    m = HBAR ** 2 / (2 * t * tob.Q * a ** 2)                 # E ≈ Ec + t a^2 k^2
    F = quad(lambda x: x ** -0.5 * expit(eta - x), 0, np.inf)[0] / np.sqrt(np.pi)
    n_ref = np.sqrt(2 * m * kT * tob.Q / np.pi) / HBAR * F
    assert n == pytest.approx(n_ref, rel=1e-2)


def make_fet(carrier="n", EOT=1e-9, W=2e-9):
    k, E = w90.bands_1d(two_band_chain(gap=1.2, t=-1.0), nk=3000)
    ch = tob.Channel(E=E, a=3.1e-10, n_val=1)
    return tob.ToBFET(ch=ch, W=W, EOT=EOT, carrier=carrier, n_atoms_per_cell=1.0)


@pytest.mark.parametrize("carrier", ["n", "p"])
def test_ideal_swing_and_ioff_pin(carrier):
    fet = make_fet(carrier)
    m = fet.metrics(VDD=0.7)
    assert m["SS_mV_dec"] == pytest.approx(1e3 * np.log(10) * tob.KB * 300, rel=2e-3)   # 59.5 mV/dec
    assert m["I_off"] == pytest.approx(10e-9 / 1e-6 * 1e-3, rel=1e-2)                 # 1e-5 mA/um
    assert m["I_on"] > 100 * m["I_off"]


def test_mirror_bands_give_equal_n_and_p():
    """The toy VB is the mirror image of the CB, so n- and p-FET metrics must coincide."""
    mn, mp = make_fet("n").metrics(), make_fet("p").metrics()
    assert mp["I_on"] == pytest.approx(mn["I_on"], rel=1e-6)
    assert mp["shift"] == pytest.approx(-mn["shift"], abs=1e-6)


def test_doping_sets_mu_and_charge_control():
    fet = make_fet("n")
    assert fet.N0 == pytest.approx(fet.n_SD, rel=1e-8)
    # electrostatic charge control: q dn / dV_G cannot exceed C_G (C_Q in series)
    _, n1, _ = fet.solve(0.8, 0.7)
    _, n2, _ = fet.solve(0.9, 0.7)
    assert 0 < tob.Q * (n2 - n1) / 0.1 <= fet.C_G * 1.0001


@pytest.mark.parametrize("carrier", ["n", "p"])
def test_ion_independent_of_doping(carrier):
    """With ideal electrostatics and a fixed C_G, the doping only translates I_D(V_GS); the I_OFF pin undoes it.

    This is why the per-atom vs per-formula-unit reading of Matić's doping cannot move I_ON (signed-off
    conventions, 2026-11-11).
    """
    k, E = w90.bands_1d(two_band_chain(gap=1.2, t=-1.0), nk=3000)
    ms = [tob.ToBFET(ch=tob.Channel(E=E, a=3.1e-10, n_val=1), W=2e-9, EOT=1e-9, carrier=carrier,
                     doping=d).metrics() for d in (0.01, 0.03)]
    assert ms[1]["I_on"] == pytest.approx(ms[0]["I_on"], rel=1e-4)
    assert ms[1]["n_on"] == pytest.approx(ms[0]["n_on"], rel=1e-4)
    assert ms[1]["muS"] != pytest.approx(ms[0]["muS"], abs=1e-3)      # the doping does move the Fermi level


def test_thinner_oxide_more_on_current():
    assert make_fet(EOT=0.6e-9).metrics()["I_on"] > make_fet(EOT=1e-9).metrics()["I_on"]


@pytest.mark.parametrize("carrier", ["n", "p"])
def test_wider_sweep_same_ion(carrier):
    """The gate-sweep span only brackets the I_OFF point: a wider span at the same step leaves I_ON unchanged."""
    f = make_fet(carrier=carrier)
    a, b = f.metrics(), f.metrics(span=4.0, n_sweep=381)
    assert b["I_on"] == pytest.approx(a["I_on"], rel=1e-3)
    assert b["I_off"] == pytest.approx(a["I_off"], rel=1e-3)
