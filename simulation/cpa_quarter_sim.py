"""
Monte Carlo simulation: how much does uneven quartering of a cyproterone acetate
(CPA) 50 mg tablet change steady-state plasma exposure?

Eden Openly — ER-027 §12 appendix. Version 1.0, 1 November 2026 (UTC).

Approach
- Linear pharmacokinetics, so results are expressed relative to an even 12.5 mg daily
  reference and no absolute concentration is needed.
- Disposition: oral bi-exponential (two-compartment with first-order absorption),
  reconstructed from published single-dose summary statistics for Androcur 50 mg
  (AA Pharma 2010 monograph, comparative bioavailability table, Androcur arm):
  Cmax 134.1 ng/mL, tmax 2.50 h, AUC0-72 2148.5 ng·h/mL, AUCinf 2806.8 ng·h/mL,
  t1/2 38.21 h. A one-compartment model using only the terminal half-life is run as
  a sensitivity check.
- Quarter content: each tablet (content ~ N(50 mg, CV 2%)) is reduced by a mass-loss
  fraction, then divided into four shares by a symmetric Dirichlet distribution tuned
  to a target quarter-to-quarter CV.
- Dosing scenarios: (A) home splitting, the four quarters of one tablet on four
  consecutive days in random order; (B) pharmacy pooled, every daily quarter from a
  different tablet; (C) adverse ordering, largest quarter first, smallest last.
- Benchmarks: even 12.5 mg daily with one missed dose per 28 days; 25 mg every other day.
"""
import json
import numpy as np
from scipy.optimize import least_squares

RNG = np.random.default_rng(20261101)

# ---------------------------------------------------------------- PK model
SRC = dict(Cmax=134.074, tmax=2.50, AUC72=2148.491, AUCinf=2806.800, thalf=38.21)  # Androcur arm


def biexp(t, A, alpha, B, beta, ka):
    t = np.asarray(t, float)
    return A * np.exp(-alpha * t) + B * np.exp(-beta * t) - (A + B) * np.exp(-ka * t)


def fit_biexp(src=SRC):
    beta = np.log(2) / src["thalf"]

    def auc(A, alpha, B, ka, T=None):
        if T is None:
            return A / alpha + B / beta - (A + B) / ka
        return (A / alpha * (1 - np.exp(-alpha * T)) + B / beta * (1 - np.exp(-beta * T))
                - (A + B) / ka * (1 - np.exp(-ka * T)))

    def resid(p):
        A, alpha, B, ka = np.exp(p)
        tg = np.linspace(0, 24, 4801)
        c = biexp(tg, A, alpha, B, beta, ka)
        i = np.argmax(c)
        return [
            (c[i] - src["Cmax"]) / src["Cmax"],
            (tg[i] - src["tmax"]) / src["tmax"],
            (auc(A, alpha, B, ka, 72) - src["AUC72"]) / src["AUC72"],
            (auc(A, alpha, B, ka) - src["AUCinf"]) / src["AUCinf"],
        ]

    best = None
    for a0 in (0.1, 0.2, 0.4, 0.8):
        for k0 in (0.5, 1.0, 2.0):
            r = least_squares(resid, np.log([200, a0, 45, k0]))
            if best is None or r.cost < best.cost:
                best = r
    A, alpha, B, ka = np.exp(best.x)
    return dict(A=A, alpha=alpha, B=B, beta=beta, ka=ka, fit_resid=list(map(float, resid(best.x))))


def onecomp_params(thalf=SRC["thalf"], tmax=SRC["tmax"]):
    ke = np.log(2) / thalf
    # solve ka from tmax = ln(ka/ke)/(ka-ke)
    from scipy.optimize import brentq
    ka = brentq(lambda k: np.log(k / ke) / (k - ke) - tmax, ke * 1.01, 50)
    return dict(A=0.0, alpha=1.0, B=1.0, beta=ke, ka=ka)


# ---------------------------------------------------------------- simulation core
DT = 0.5                      # h, within-day grid
OFFS = np.arange(0, 24, DT)   # 48 offsets per day
NDAYS = 120
BURN = 30                     # discard days 1-30


def kernels(p):
    """H[o, k] = concentration per mg at time offset o, k days after a dose."""
    k = np.arange(NDAYS)
    t = OFFS[:, None] + 24 * k[None, :]
    return biexp(t, p["A"], p["alpha"], p["B"], p["beta"], p["ka"])


def concentrations(doses, H):
    """doses (N, NDAYS) mg given at start of each day.
    Returns dict: trough (N, NDAYS) = pre-dose value at start of each day,
    peak (N, NDAYS) = max within each day, mean (N,) = mean over days after burn-in."""
    N = doses.shape[0]
    idx = np.arange(NDAYS)
    lag = idx[:, None] - idx[None, :]
    peak = np.full((N, NDAYS), -np.inf)
    tot = np.zeros(N)
    trough = None
    for o in range(len(OFFS)):
        T = np.where(lag >= 0, H[o][np.clip(lag, 0, None)], 0.0)
        c = doses @ T.T
        if o == 0:
            trough = c
        peak = np.maximum(peak, c)
        tot += c[:, BURN:].sum(axis=1)
    return dict(trough=trough, peak=peak, mean=tot / ((NDAYS - BURN) * len(OFFS)))


def metrics(C, Cref):
    """Trough of day j = pre-dose value at start of day j+1; peak = within-day max.
    All relative to the even-split reference."""
    tr = C["trough"][:, BURN + 1:NDAYS] / Cref["trough"][0, BURN + 1:NDAYS]
    pk = C["peak"][:, BURN:NDAYS - 1] / Cref["peak"][0, BURN:NDAYS - 1]
    return dict(cavg_rel=C["mean"] / Cref["mean"][0], trough_rel=tr, peak_rel=pk)


def summarise(m):
    tr, pk, ca = m["trough_rel"], m["peak_rel"], m["cavg_rel"]
    q = lambda x, p: float(np.percentile(x, p))
    worst_tr = tr.min(axis=1)
    best_tr = tr.max(axis=1)
    return {
        "Cavg % of ref, median [2.5-97.5 pct of patients]":
            [round(100 * q(ca, 50), 1), round(100 * q(ca, 2.5), 1), round(100 * q(ca, 97.5), 1)],
        "daily trough % of ref, 2.5-97.5 pct of all days":
            [round(100 * q(tr, 2.5), 1), round(100 * q(tr, 97.5), 1)],
        "daily peak % of ref, 2.5-97.5 pct of all days":
            [round(100 * q(pk, 2.5), 1), round(100 * q(pk, 97.5), 1)],
        "% of days trough <90% of ref": round(100 * float((tr < 0.9).mean()), 2),
        "% of days trough <80% of ref": round(100 * float((tr < 0.8).mean()), 2),
        "% of days peak outside 80-125% of ref": round(100 * float(((pk < 0.8) | (pk > 1.25)).mean()), 2),
        "patient's lowest trough over 90 days, median % of ref": round(100 * q(worst_tr, 50), 1),
        "patient's lowest trough, 2.5th pct of patients": round(100 * q(worst_tr, 2.5), 1),
        "patient's highest trough, 97.5th pct of patients": round(100 * q(best_tr, 97.5), 1),
    }


# ---------------------------------------------------------------- dose generators
def dirichlet_alpha(cv):
    # Dirichlet(a,a,a,a): CV of each share = sqrt(3/(4a+1))
    return (3 / cv ** 2 - 1) / 4


def tablets(n, rng):
    return rng.normal(50.0, 1.0, n)  # CV 2%


def quarters_home(N, cv, loss, rng, order="random"):
    ntab = NDAYS // 4
    tab = tablets(N * ntab, rng) * (1 - loss)
    a = dirichlet_alpha(cv)
    sh = rng.dirichlet([a] * 4, size=N * ntab)
    if order == "largest_first":
        sh = -np.sort(-sh, axis=1)
    q = tab[:, None] * sh
    return q.reshape(N, ntab * 4)


def quarters_pooled(N, cv, loss, rng):
    a = dirichlet_alpha(cv)
    tab = tablets(N * NDAYS, rng) * (1 - loss)
    share = rng.beta(a, 3 * a, size=N * NDAYS)   # marginal of Dirichlet share
    return (tab * share).reshape(N, NDAYS)


def even(N):
    return np.full((N, NDAYS), 12.5)


def missed_one_per_28(N, rng):
    """One missed dose per 28 days, at a uniform rate throughout.

    v1.1: the last block (days 112-119) is only 8 days long. v1.0 gave it a missed
    dose in every patient, which raised the rate in the analysed window above one
    per 28 days. It now gets one with probability 8/28. The partial-block draw uses
    its own generator, so the main RNG stream, and every other output, is unchanged.
    """
    d = even(N)
    rng_partial = np.random.default_rng(2026110128)
    rows_all = np.arange(N)
    for start in range(0, NDAYS, 28):
        end = min(start + 28, NDAYS)
        days = rng.integers(start, end, N)
        rows = rows_all
        if end - start < 28:
            keep = rng_partial.random(N) < (end - start) / 28
            rows, days = rows[keep], days[keep]
        d[rows, days] = 0.0
    return d


def alt_day_25(N):
    d = np.zeros((N, NDAYS))
    d[:, ::2] = 25.0
    return d


# ---------------------------------------------------------------- run
def run(p, label, N=10000):
    H = kernels(p)
    Cref = concentrations(even(1), H)
    out = {"model": label, "params": {k: float(v) for k, v in p.items() if k != "fit_resid"}}
    out["ref_steady_state_trough_over_peak"] = round(float(Cref["trough"][0, -1] / Cref["peak"][0, -2]), 3)
    res = {}
    for cv in (0.10, 0.15, 0.20, 0.30):
        for loss in (0.0, 0.03, 0.06, 0.10):
            key = f"CV{int(cv*100)}_loss{int(loss*100)}"
            res[f"A_home_{key}"] = summarise(metrics(concentrations(quarters_home(N, cv, loss, RNG), H), Cref))
            res[f"B_pooled_{key}"] = summarise(metrics(concentrations(quarters_pooled(N, cv, loss, RNG), H), Cref))
            res[f"C_largestfirst_{key}"] = summarise(
                metrics(concentrations(quarters_home(N, cv, loss, RNG, "largest_first"), H), Cref))
    res["BENCH_missed_1_per_28"] = summarise(metrics(concentrations(missed_one_per_28(N, RNG), H), Cref))
    res["BENCH_25mg_alt_day"] = summarise(metrics(concentrations(alt_day_25(1), H), Cref))
    out["results"] = res
    return out, Cref


if __name__ == "__main__":
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
    p2 = fit_biexp()
    print("bi-exponential fit:", json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in p2.items()}))
    allout = {"fit": p2}
    allout["main"], _ = run(p2, "two-compartment (fitted), t1/2 38.2 h", N)
    # sensitivity: keep fitted absorption and fast phase, vary terminal half-life
    for th in (30.0, 44.0, 56.0):
        pth = dict(p2, beta=np.log(2) / th)
        allout[f"sens_thalf_{int(th)}"], _ = run(pth, f"two-compartment, terminal t1/2 {th} h", N // 2)
    allout["sens_onecomp"], _ = run(onecomp_params(), "one-compartment, t1/2 38.2 h", N // 2)
    with open("results.json", "w") as f:
        json.dump(allout, f, indent=1, default=float)
    print("done")
