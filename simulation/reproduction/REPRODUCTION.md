# Reproduction record

## 1. Exact re-run of the published code — 1 November 2026 (UTC)

- **Environment:** Python 3.13, numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2, on Linux.
- **Command:** `python cpa_quarter_sim.py 10000`
- **`results.json`:** all 3,294 values identical to the published file.
- **`plot_fig.py`:** figure identical pixel for pixel (1440 × 1480 px). The file checksum differs only because PNG metadata varies between runs.

## 2. Independent re-implementation — 1 November 2026 (UTC)

`check_independent.py` was written separately, from the parameters stated in the appendix rather than from the published code. It uses 2,000 simulated patients per condition and a different random seed, so agreement to about one percentage point, not exact equality, is the expected result.

Middle 95% of daily troughs, as % of the even-split level, and % of days below 80%:

| Scenario | Independent | Published |
|---|---|---|
| Home split, CV 15% | 90–111%; 0.01% | 90–110%; 0.01% |
| Pharmacy pooled, CV 15% | 87–114%; 0.1% | 87–114%; 0.1% |
| Home split, CV 30% | 81–121%; 2.1% | 81–121%; 2.1% |
| Pharmacy pooled, CV 30% | 74–129%; 7.1% | 74–129%; 7.0% |
| One dose missed per 28 days | 64–100%; 8.5% | 64–100%; 8.5% |
| 25 mg every other day | 78–122%; 49.4% | 78–122%; 49.4% |
| Accumulation factor | 2.83 | 2.83 |

Both reproductions were carried out with AI assistance (Claude, Anthropic). They have not yet been repeated by a human.
