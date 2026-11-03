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
| Terminal-phase accumulation factor, 1/(1 − e^(−24β)) | 2.83 | 2.83 |

## 3c. Figure label rounding (v1.3) — 3 November 2026 (UTC)

Panel B of the figure labelled the home-split, CV 15% interval "90–110%". `plot_fig.py` formatted every label with `:.0f`, which rounds a half to the nearest even number, so 110.5 printed as 110. ER-027's Appendix Q gives this interval to one decimal place, 89.9–110.5%, because the rounded "90–110%" had been read as "within 10% on 95% of days", which is narrower than the result. `plot_fig.py` now labels that row to one decimal place and every other row to whole per cent, as the appendix does. No other label is affected: no other interval endpoint ends in exactly .5. No model output changes; `results.json` and `calibration_sensitivity.json` are unchanged. Re-drawn figure checked against v1.2: only the pixels of that one label differ.

## 3b. Calibration sensitivity and a label correction (v1.2) — 1 November 2026 (UTC)

An external reviewer found that `fit_biexp()` takes the peak time from a fixed grid with `argmax`, so the optimiser gets no gradient from that residual. `calibration_sensitivity.py` refits with a continuous peak time from 400 random starts, profiles the fit over fixed distribution half-lives, and reruns the key scenarios (N = 3,000) for selected profile calibrations with a squared-residual sum up to 25% above the main one's (an exploratory rule; the 4 h profile is about 21% worse, the others better). The reviewer's independent fit (squared-residual sum about 0.000335, against 0.00163–0.00169 for the published fit) is reproduced: 0.00034, with absorption and distribution rates converging. Results across the selected calibrations, distribution half-life 1.55–4 h:

| Scenario | Range across calibrations |
|---|---|
| Home splitting, CV 15%: daily trough, middle 95% | 90.1–110.3% to 89.5–111.0% |
| Pooled, CV 30%: days with trough <80% | 6.85–7.51% |
| One missed per 28 days: daily trough, middle 95% | 61.7–100% to 64.6–100% |
| Home splitting, CV 15%: daily peak, middle 95% | 81.6–120.3% to 82.5–119.4% |

The main model is unchanged. The `results.json` key "lowest trough over 90 days" is renamed "over days 32–120 (89 troughs)", which is what it measures; every value is identical to v1.1.

## 3a. Note on the accumulation factor — 1 November 2026 (UTC)

The 2.83 in section 2 is the accumulation factor of the terminal exponential alone. An external reviewer pointed out that it is not the accumulation of the whole fitted curve. Summing the fitted single-dose curve at 24-hour intervals gives: steady-state peak 199.7 against single-dose peak 133.6 ng/mL (1.49-fold); 24-hour exposure 1.94-fold; trough 2.77-fold. The agreement of 2.83 with the label's "about 3" is a plausibility check on the terminal phase, not a validation of the model. No output changes.

## 3. Correction to the missed-dose benchmark (v1.1) — 1 November 2026 (UTC)

An external reviewer found that the benchmark "one dose missed per 28 days" gave every patient a missed dose in the last block, days 112–119, which is only 8 days long. That raised the missed-dose rate in the analysed window (days 31–120) above one per 28 days. The independent re-implementation had the same flaw, so its agreement in section 2 did not test this point.

In v1.1, both scripts give the short last block a missed dose with probability 8/28, so the rate is one per 28 days throughout. In `cpa_quarter_sim.py` that draw uses its own generator, so the main random stream and every other output are unchanged: only the `BENCH_missed_1_per_28` entries of `results.json` differ, in the main run and the four sensitivity runs. The figure is unchanged pixel for pixel.

| One dose missed per 28 days | v1.0 | v1.1 | Independent, v1.1 |
|---|---|---|---|
| Daily trough, middle 95% of days | 64–100% | 64–100% | 64–100% |
| Days with trough <90% | 12.7% | 10.8% | — |
| Days with trough <80% | 8.5% | 7.2% | 7.1% |
| Lowest trough over days 32–120 (89 troughs), median | 64% | 64% | — |

Both reproductions were carried out with AI assistance (Claude, Anthropic). They have not yet been repeated by a human.
