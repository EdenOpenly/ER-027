# ER-027 — companion repository

Code and data for **ER-027**, an Eden Openly Evidence Review on cyproterone acetate dosing in Korea: *Cyproterone acetate in Korea: from the 50 mg tablet to the guideline dose*, version 1.0, doi:[10.5281/zenodo.23909671](https://doi.org/10.5281/zenodo.23909671).

The tag `v1.0` marks the state of this repository that ER-027 version 1.0 reports. The code inside it is at release 1.3 (the simulation script itself is unchanged since 1.2).

## What is here

`simulation/` holds the pharmacokinetic simulation reported in ER-027's Appendix Q, which supports its practice section (§12). Korea sells only a 50 mg cyproterone acetate tablet, so a 12.5 mg daily dose means taking a quarter of a tablet. The simulation asks how much uneven quartering, and tablet lost as powder when splitting, changes the drug concentrations that reach plasma.

- `cpa_quarter_sim.py` — the simulation (version 1.2). It writes `results.json`.
- `calibration_sensitivity.py` — refits the model with a continuous peak time, profiles the fit, and reruns the key scenarios across selected profile calibrations (an exploratory analysis, not an exhaustive one). It writes `calibration_sensitivity.json`.
- `plot_fig.py` — draws the figure (Figure Q1) from the model and `results.json` (version 1.3).
- `results.json` — the full published outputs.
- `ER-027-appendix-quartering-simulation-fig1.png` — the published figure.
- `requirements.txt` — the package versions the published outputs were produced with (Python 3.13).
- `INPUTS.md` — every model input and its source.
- `SHA256SUMS` — checksums of the published files.
- `reproduction/` — an independent re-implementation and the record of both reproductions.

## How to reproduce

```
cd simulation
pip install -r requirements.txt
python cpa_quarter_sim.py 10000   # about 2 minutes; overwrites results.json
python plot_fig.py                # overwrites the figure
python calibration_sensitivity.py # a few minutes; writes calibration_sensitivity.json
sha256sum -c SHA256SUMS
```

The random seed is fixed, so with the pinned package versions `results.json` is reproduced exactly. The figure is reproduced pixel for pixel, but its file checksum can differ because PNG metadata varies between runs (see `reproduction/REPRODUCTION.md`).

## What the simulation is, and is not

It is a model, not a measurement: an illustrative model calibrated to the published group summary statistics of one single-dose bioavailability study of Androcur 50 mg in 18 healthy men. Those four statistics do not uniquely fix the concentration–time curve or its absorption and distribution parameters, so predictions are conditional on the calibration chosen; `calibration_sensitivity.py` shows how much they move across a selection of calibrations that fit about as well. It assumes linear kinetics at 12.5 mg and does not model testosterone. The splitting inputs (quarter-to-quarter variation, mass loss) are borrowed from studies of other drugs' tablets. No study identified in ER-027 has measured them for Androcur. ER-027 reports the simulation in Appendix Q as context for its §12 and gives it no certainty grade.

## Licence

Code: MIT (see `LICENSE`). Text, data and figure: Creative Commons Attribution 4.0 International (CC BY 4.0).

## Contact

Eden Openly — editor@edenopenly.com · https://edenopenly.com. Errors are corrected in public: https://edenopenly.com/editorial/corrections/
