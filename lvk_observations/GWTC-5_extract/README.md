# GWTC-5.0 BBH population results — small extract for GROWL

`gwtc5_bbh_rate_and_mass.h5` (6.7 MB) holds the two LVK population results that the GROWL
onboarding tutorials compare simulations against:

| quantity | group | what it is |
|---|---|---|
| BBH merger rate density $\mathcal R(z)$ | `default_bbh/rate_vs_redshift`, `pixelpop/rate_vs_redshift` | Gpc⁻³ yr⁻¹ on a redshift grid |
| primary-mass distribution at $z = 0.2$ | `default_bbh/dR_dm1_z0.2`, `pixelpop/dR_dm1_z0.2` | $\mathrm d\mathcal R/\mathrm dm_1$ in Gpc⁻³ yr⁻¹ M☉⁻¹ on an $m_1$ grid |

Each group has `positions` (the grid), `quantiles` (rows = 5 %, 50 %, 95 % over the posterior
draws) and `draws` (500 thinned posterior draws, float32). `default_bbh/lamb` are the posterior
samples of the redshift-evolution index $\kappa_z$ ($\mathcal R \propto (1+z)^{\kappa_z}$).
File-level attributes record the 1st/99th percentiles of the events' posterior samples in $z$
and $m_1$ (the region where the data actually constrain the fits; the paper hatches outside it).

"Default BBH" is the LVK parametric model (Broken Power Law + Two Peaks in mass, power-law in
redshift); "PixelPop" is their non-parametric model. Both are as plotted in Figs. 2 and 5 of the
GWTC-5.0 population paper.

## Provenance

Everything here is derived from the **LVK GWTC-5.0 population data release**,
https://doi.org/10.5281/zenodo.20292639 — please cite the GWTC-5.0 population paper and that DOI,
not this repository, when you use these numbers. The full release is 25 GB compressed (42 GB
unpacked); the extract was made with `extract_gwtc5_bbh_rate_and_mass.py` (in this folder) from
these files of the release:

| file (inside `popsummary_files.tar.gz`) | size | used for |
|---|---|---|
| `gwtc5_updated_default_mmax_mass_TwoPeakBrokenPowerLawSmoothedMassDistribution_redshift_PowerLawRedshift_magnitude_iid_spin_magnitude_gaussian_tilt_iid_spin_orientation_popsummary_result.h5` | 727 MB | Default BBH $\mathcal R(z)$ and $m_1$ |
| `may26_datarelease/z_varcut1/z_varcut1_popsummary.h5` | 151 MB | PixelPop $\mathcal R(z)$ |
| `may26_datarelease/m1m2_varcut1/m1m2_varcut1_popsummary.h5` | 295 MB | PixelPop $m_1$ |
| `figure_scripts/one_percent_pe_limits.json` (in `scripts.tar.gz`) | 50 kB | PE percentile limits |

The recipe for turning the release's `rates_on_grids` into these quantities (in particular the
$(1+z)^{\kappa_z}$ factor that moves the $z=0$ mass distribution to $z=0.2$, and the
$\mathrm dR/\mathrm d\ln m_1 \to \mathrm dR/\mathrm dm_1$ conversion for PixelPop) is copied from
the release's own `figure_scripts/make_fig_2.ipynb` and `make_fig_5.ipynb`.

To regenerate, download and unpack the release, then:

```bash
python extract_gwtc5_bbh_rate_and_mass.py /path/to/unpacked/GWTC-5
```

Reading the extract needs only `h5py`:

```python
import h5py
with h5py.File("gwtc5_bbh_rate_and_mass.h5") as f:
    z = f["default_bbh/rate_vs_redshift/positions"][()]
    R_lo, R_med, R_hi = f["default_bbh/rate_vs_redshift/quantiles"][()]
```
