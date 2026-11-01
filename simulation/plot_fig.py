import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cpa_quarter_sim as m
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#1f1f1e", "#6b6a65", "#e4e3de"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED})
rng = np.random.default_rng(7)
p = m.fit_biexp(); H = m.kernels(p); Cref = m.concentrations(m.even(1), H)
def troughs(d):
    C = m.concentrations(d, H)
    return 100 * C["trough"][0, m.BURN+1:] / Cref["trough"][0, m.BURN+1:]
series = [("Home split, CV 15%", troughs(m.quarters_home(1, .15, 0, rng)), BLUE),
          ("Pharmacy pooled, CV 30%", troughs(m.quarters_pooled(1, .30, 0, rng)), ORANGE),
          ("Even dose, 1 missed per 28 days", troughs(m.missed_one_per_28(1, rng)), AQUA)]
d = json.load(open("results.json"))["main"]["results"]
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 7.4), gridspec_kw={"height_ratios": [1, 1.25]})
days = np.arange(m.BURN + 2, m.NDAYS + 1)[:60]
for lab, y, c in series:
    a1.plot(days, y[:60], color=c, lw=2, label=lab)
a1.axhline(100, color=MUTED, lw=1, ls=(0, (3, 3)))
a1.set_ylabel("Daily trough, % of even-split level"); a1.set_xlabel("Day of treatment")
a1.set_title("A. One simulated patient (two-compartment model, no mass loss)", loc="left", color=INK, fontsize=10)
a1.legend(frameon=False, loc="lower left", ncol=1, fontsize=8)
a1.grid(axis="y", color=GRID, lw=0.8); a1.set_axisbelow(True)
for s in ("top", "right"): a1.spines[s].set_visible(False)
# Panel B
rows = []
for cv in (10, 15, 20, 30):
    rows.append((f"Home split, CV {cv}%", d[f"A_home_CV{cv}_loss0"], BLUE))
    rows.append((f"Pharmacy pooled, CV {cv}%", d[f"B_pooled_CV{cv}_loss0"], ORANGE))
rows.append(("Even dose, 1 missed per 28 days", d["BENCH_missed_1_per_28"], AQUA))
rows.append(("25 mg every other day", d["BENCH_25mg_alt_day"], AQUA))
key = "daily trough % of ref, 2.5-97.5 pct of all days"
a2.axvspan(62.9, 159.0, color="#f0efec", zorder=0)
a2.text(158, -0.9, "Between-person exposure range\n(95%, AUC CV 24%)", ha="right", va="center", color=MUTED, fontsize=7.5)
for i, (lab, r, c) in enumerate(rows):
    lo, hi = r[key]; y = len(rows) - 1 - i
    a2.plot([lo, hi], [y, y], color=c, lw=3, solid_capstyle="round")
    a2.text(hi + 2, y, f"{lo:.0f}–{hi:.0f}%", va="center", color=INK, fontsize=7.5)
a2.set_yticks(range(len(rows))); a2.set_yticklabels([r[0] for r in rows[::-1]], fontsize=8)
a2.axvline(100, color=MUTED, lw=1, ls=(0, (3, 3)))
a2.set_xlim(50, 175); a2.set_ylim(-1.6, len(rows) - 0.4)
a2.set_xlabel("Daily trough, % of even-split level (middle 95% of days)")
a2.set_title("B. Spread of daily troughs, 10,000 simulated patients, no mass loss", loc="left", color=INK, fontsize=10)
a2.grid(axis="x", color=GRID, lw=0.8); a2.set_axisbelow(True)
for s in ("top", "right"): a2.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("ER-027-appendix-quartering-simulation-fig1.png", dpi=200)
