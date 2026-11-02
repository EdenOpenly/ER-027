# Model inputs and their sources

## Pharmacokinetics (fitted in `fit_biexp()`)

Fitted to the Androcur 50 mg arm of a comparative bioavailability study (18 healthy men, single dose, arithmetic means), as reported in: AA Pharma Inc., *Cyproterone (cyproterone acetate tablets BP) 50 mg — Product Monograph*, 2010.

| Input | Value | Fitted |
|---|---|---|
| Cmax | 134.074 ng/mL | 133.6 ng/mL |
| tmax | 2.50 h | 2.43 h |
| AUC0–72 | 2148.491 ng·h/mL | 2195 ng·h/mL |
| AUCinf | 2806.800 ng·h/mL | 2753 ng·h/mL |
| Terminal half-life | 38.21 h | fixed |

Fitted parameters: A 252.47, α 0.23667 h⁻¹, B 37.392, β 0.01814 h⁻¹ (fixed), kₐ 0.77359 h⁻¹. That is an absorption half-life of 0.90 h and a distribution-phase half-life of 2.9 h, with the distribution phase carrying 34% of disposition AUC.

**Calibration, not recovery (v1.2 note).** `fit_biexp()` reads the peak time off a fixed grid with `argmax`, so that residual gives the optimiser no gradient, and the four group summary statistics do not identify four parameters in any case. Refitted with a continuous peak time (`calibration_sensitivity.py`), the best free fit drives the absorption and distribution rates together (both half-lives about 1.55 h) and fits the four statistics more closely (squared-residual sum 0.00034 against 0.00169). Fixing the distribution half-life anywhere from about 1.5 to 4 h fits as well as the main calibration. Across those calibrations, the main trough results move by about one percentage point (home splitting, CV 15%: 89.5–111.0% at the widest; pooled, CV 30%: 6.9–7.5% of days below 80%). The main model keeps the v1.0 calibration, which lies inside that family; its absorption and distribution half-lives are calibration values, not measured ones. The share of exposure in the fast phase depends on convention: 34% counting the disposition terms' gross areas, 27% net of the absorption term.

### Cross-checks against product labels

- **Androcur 50 SmPC, Netherlands** (CBG h06540, version 08-2022): terminal half-life 43.9 ± 12.8 h.
- **Androcur 10 mg labels, Germany** (Advanz Pharma, Fachinformation 09.2024) **and Israel** (Bayer Israel, May 2022):
  - biphasic decline, half-lives about 0.8 h and 2.3 days;
  - two- to two-and-a-half-fold accumulation with daily dosing.
- **Androcur 50/100 mg package insert, South Africa** (Bayer, SAHPRA-approved, 2009): accumulation "about 3".

The model's terminal-phase accumulation factor, 1/(1 − e^(−24β)), is 2.83, consistent with the labels' figures, which are also based on the terminal half-life. That is a plausibility check on the terminal phase, not a validation of the whole curve: summed over daily doses, the fitted model's steady-state peak is 1.49 times the single-dose peak, 24-hour exposure is 1.94 times the first day's, and the trough rises 2.77-fold. Its 38.21 h terminal half-life is the shortest of these sources, so it shows the most day-to-day variation. The sensitivity runs use 30, 44 and 56 h.

## Splitting inputs

Quarter-to-quarter coefficients of variation of 10, 15, 20 and 30% were simulated, and mass losses of 0, 3, 6 and 10%. They are anchored to quartering studies of other drugs:

- **Madathilethu J, et al.** *BMJ Paediatr Open* 2018;2(1):e000198. doi:10.1136/bmjpo-2017-000198. Hydrocortisone 10 mg, four operators, commercial splitter: content CV 15.2%, about 2% of mass lost.
- **Endocrine Abstracts** 2009;23:OC3.1 (BSPED 2009). Hydrocortone 10 mg quarters: content CV 29.7%, weight CV 13.4%.
- **Dosti B, Malaj L, Xhafaj D.** *Int J Pharm Pharm Sci* 2016;8(6):222–224. Enalapril 20 mg, scored on one side, quarters: weight RSD 15.2%.

Tablet content: normal, mean 50 mg, CV 2% (assumed).
