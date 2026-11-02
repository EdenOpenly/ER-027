"""Calibration sensitivity for the ER-027 Appendix Q model (v1.2).

The main model (cpa_quarter_sim.py) is calibrated by fit_biexp(), whose peak-time
residual is read off a fixed 0.005 h grid with argmax. That residual is piecewise
constant, so the optimiser gets no gradient from it and the reported absorption and
distribution parameters are not securely identified. Four group summary statistics
(Cmax, tmax, AUC0-72, AUCinf) cannot in any case pin down a four-parameter curve.

This script:
  1. computes tmax continuously (root of dC/dt) and refits from 400 random starts;
  2. profiles the fit over fixed distribution half-lives (1-12 h), refitting the rest;
  3. reruns the key scenarios for every calibration whose squared-residual sum is no
     worse than the published fit's, and for the published fit itself.

Output: calibration_sensitivity.json. Run time a few minutes; N = 3,000 per scenario.
"""
import json
import numpy as np
from scipy.optimize import least_squares, brentq
import cpa_quarter_sim as m

SRC = m.SRC
BETA = np.log(2) / SRC["thalf"]


def C(t, A, al, B, ka):
    return A * np.exp(-al * t) + B * np.exp(-BETA * t) - (A + B) * np.exp(-ka * t)


def dC(t, A, al, B, ka):
    return -A * al * np.exp(-al * t) - B * BETA * np.exp(-BETA * t) + (A + B) * ka * np.exp(-ka * t)


def tmax(A, al, B, ka):
    g = np.linspace(1e-3, 24, 2001)
    d = dC(g, A, al, B, ka)
    i = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
    return brentq(dC, g[i[0]], g[i[0] + 1], args=(A, al, B, ka)) if len(i) else np.nan


def auc(A, al, B, ka, T=None):
    if T is None:
        return A / al + B / BETA - (A + B) / ka
    return (A / al * (1 - np.exp(-al * T)) + B / BETA * (1 - np.exp(-BETA * T))
            - (A + B) / ka * (1 - np.exp(-ka * T)))


def resid(A, al, B, ka):
    tm = tmax(A, al, B, ka)
    if not np.isfinite(tm):
        return [10.0] * 4
    return [(C(tm, A, al, B, ka) - SRC["Cmax"]) / SRC["Cmax"], (tm - SRC["tmax"]) / SRC["tmax"],
            (auc(A, al, B, ka, 72) - SRC["AUC72"]) / SRC["AUC72"],
            (auc(A, al, B, ka) - SRC["AUCinf"]) / SRC["AUCinf"]]


def ssr(A, al, B, ka):
    return float(np.sum(np.square(resid(A, al, B, ka))))


def fast_net_share(A, al, B, ka):
    """Net area of A(e^-at - e^-kat) as a share of total AUC (one convention of several)."""
    return A * (1 / al - 1 / ka) / auc(A, al, B, ka)


out = {}
pub = m.fit_biexp()
pp = (pub["A"], pub["alpha"], pub["B"], pub["ka"])
out["published_fit"] = dict(A=pp[0], alpha=pp[1], B=pp[2], ka=pp[3], ssr_continuous_tmax=ssr(*pp),
                            tmax=tmax(*pp), fast_net_share=fast_net_share(*pp),
                            fast_gross_share=(pp[0] / pp[1]) / (pp[0] / pp[1] + pp[2] / BETA))

# 1. free multistart
rng = np.random.default_rng(1)
best = None
for _ in range(400):
    x0 = np.log([rng.uniform(20, 2000), rng.uniform(0.03, 3), rng.uniform(5, 80), rng.uniform(0.1, 5)])
    try:
        r = least_squares(lambda p: resid(*np.exp(p)), x0, max_nfev=4000)
    except Exception:
        continue
    if best is None or r.cost < best.cost:
        best = r
bA, bal, bB, bka = np.exp(best.x)
out["free_best_fit"] = dict(A=bA, alpha=bal, B=bB, ka=bka, ssr=ssr(bA, bal, bB, bka),
                            note="absorption and distribution rates converge (alpha ~ ka): not identified")

# 2. profile over distribution half-life
prof = []
for th in (1.0, 1.55, 2.0, 2.9, 4.0, 6.0, 8.0, 12.0):
    al = np.log(2) / th
    pb = None
    for A0 in (50, 200, 800, 5000):
        for k0 in (0.2, 0.5, 1, 2, 4):
            try:
                r = least_squares(lambda p: resid(np.exp(p[0]), al, np.exp(p[1]), np.exp(p[2])),
                                  np.log([A0, 40, k0]), max_nfev=3000)
            except Exception:
                continue
            if pb is None or r.cost < pb.cost:
                pb = r
    A, B, ka = np.exp(pb.x)
    prof.append(dict(dist_thalf_h=th, A=A, alpha=al, B=B, ka=ka, abs_thalf_h=np.log(2) / ka,
                     ssr=ssr(A, al, B, ka), tmax=tmax(A, al, B, ka), fast_net_share=fast_net_share(A, al, B, ka)))
out["profile"] = prof

# 3. key scenarios across plausible calibrations
cals = [("published fit (main model)", dict(A=pp[0], alpha=pp[1], B=pp[2], beta=BETA, ka=pp[3]))]
for d in prof:
    if d["ssr"] <= out["published_fit"]["ssr_continuous_tmax"] * 1.25:
        cals.append((f"distribution t1/2 {d['dist_thalf_h']} h", dict(A=d["A"], alpha=d["alpha"], B=d["B"],
                                                                         beta=BETA, ka=d["ka"])))
N = 3000
sens = {}
for name, p in cals:
    r = np.random.default_rng(20261101)
    H = m.kernels(p)
    Cref = m.concentrations(m.even(1), H)
    row = {}
    for lab, doses in (("home CV15", m.quarters_home(N, 0.15, 0.0, r)),
                       ("pooled CV30", m.quarters_pooled(N, 0.30, 0.0, r)),
                       ("missed 1 per 28", m.missed_one_per_28(N, r)),
                       ("25 mg alternate days", m.alt_day_25(1))):
        s = m.summarise(m.metrics(m.concentrations(doses, H), Cref))
        tk = [k for k in s if k.startswith("daily trough")][0]
        pk = [k for k in s if k.startswith("daily peak")][0]
        row[lab] = dict(trough_95=s[tk], peak_95=s[pk], pct_days_trough_lt80=s["% of days trough <80% of ref"])
    sens[name] = row
out["scenario_sensitivity"] = sens

with open("calibration_sensitivity.json", "w") as f:
    json.dump(out, f, indent=1, default=float)
print(json.dumps({k: sens[k] for k in sens}, indent=1, default=float))
