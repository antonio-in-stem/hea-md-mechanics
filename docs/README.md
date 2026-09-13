# Technical Manual and Documentation

[Article (PDF)](main.pdf) · [Supplement (PDF)](supplement.pdf)

This document provides complete technical specifications for the **hea-md-mechanics** study (Hf-Nb-Ta-Ti-Zr atomistic compression), covering mathematical methods, data dictionary, computational reproducibility, and parameter sources.

---

## 1. Analysis Methods

### Loading Records & Conventions
Engineering compressive strain $\varepsilon_c$ and stress $\sigma_{c,\mathrm{eng}}$ are positive in compression:

$$
\varepsilon_c = 1 - \frac{L_z}{L_{z0}}, \qquad \sigma_{c,\mathrm{eng}} \; (\mathrm{GPa}) = 10^{-4} P_{zz} \; (\mathrm{bar}) \frac{L_x L_y}{L_{x0} L_{y0}}
$$

LAMMPS pressure includes kinetic and virial terms. Strain rate is $\dot{\varepsilon}_c = 10^9 \; \mathrm{s}^{-1}$ (`erate = -0.001 ps^-1`). The primary comparison interval is $0 \le \varepsilon_c \le 0.10$ across all nine records (two end at 10% strain; seven extend to 20%).

### Composition Path Substitution Vector
Composition variation follows equal-fraction complementary substitution:

$$
c_{\mathrm{target}} = s, \qquad c_{i \neq \mathrm{target}} = \frac{1-s}{4}, \qquad \sum_i c_i = 1
$$

Adjusting a target element scales all four complementary species fractions simultaneously through Hf-Nb-Ta-Ti-Zr composition space.

### Mechanical Descriptors
1. **Initial Fitted Stiffness ($E_{0-2\%}$)**: Ordinary least-squares slope evaluated over $0 \le \varepsilon_c \le 0.02$ retaining the intercept $b$:

$$
\sigma_i \approx E_{0-2\%} \varepsilon_i + b
$$

Retaining intercept $b$ accounts for initial thermal and virial fluctuations near the origin.

2. **Operational 0.2% Offset Stress ($\sigma_{0.2}^{\mathrm{op}}$)**: First downward intersection of the 21-sample centered mean stress with the slope line shifted by $0.002$ strain, requiring 3-sample persistence.
3. **Common-Window Peak Stress ($\sigma_{\mathrm{p}}^{0-10\%}$)**: Maximum raw stress within $0 \le \varepsilon_c \le 0.10$.
4. **Complete-Record Maxima**: Evaluated over $0 \le \varepsilon_c \le 0.20$ for the seven extended loading records (secondary stress maxima occur between 18.41% and 19.29% strain).

### Parameter Sensitivity Design
Evaluates 60 processing settings per record (5 fit limits $\times$ 4 smoothing lengths $\times$ 3 persistence counts), producing 540 total evaluations. Sensitivity ranges measure descriptor extraction bounds across matched parameter combinations.

---

## 2. Data Dictionary & Outputs

### Primary Dataset
- 16,009 paired strain-stress observations across 9 nominal composition paths.
- Sampling spacing is $\approx 0.0001$ strain (100 steps). Every common interval contains 1,001 samples and primary fits contain 201 samples.

### Numerical Output Files (`results/`)
| Export File | Contents |
|---|---|
| `summary.csv` | Primary mechanical descriptors and relative percentage changes |
| `analysis.json` | Full campaign metadata, fit diagnostics, crossings, and 540 evaluations |
| `sensitivity.csv` | Full table of 540 curve–setting evaluations |
| `method-sensitivity.csv` | Per-case processing ranges and valid descriptor counts |
| `paired-sensitivity.csv` | Relative percentage changes under matched sample/reference extraction settings |

### Field Suffix Conventions
- `_GPa`: Values in GPa.
- `_percent`: Relative percentage changes ($100 \times \Delta$).
- Strain fields (including `fit_max`): Dimensionless fractions.

---

## 3. Reproduction & Verification

### Python Environment
Requires Python 3.10+. Install package dependencies from the repository root:
```bash
python -m pip install -e ".[test,plot]"
```

### Verification Commands
```bash
python -m hea_md verify
python -m pytest -q
```
- `hea_md verify`: Recomputes the campaign and validates JSON and CSV exports.
- `pytest`: Runs the full 179-test regression suite.

---

## 4. Model References and Sources

| Reference / Parameter Source | Role in Study |
|---|---|
| Huang et al. (2021), *Materials & Design* 202, 109560, [DOI](https://doi.org/10.1016/j.matdes.2021.109560) | MEAM model parameterization and short-range order study |
| [NIST Hf–Nb–Ta–Ti–Zr repository](https://www.ctcms.nist.gov/potentials/entry/2021--Huang-X-Liu-L-Duan-X-et-al--Hf-Nb-Ta-Ti-Zr/) | Five-element MEAM parameterization files |
| [LAMMPS documentation](https://docs.lammps.org/) | Command specifications for `units`, `fix deform`, `compute pressure`, `fix nh`, `set`, and `pair_style meam` |
