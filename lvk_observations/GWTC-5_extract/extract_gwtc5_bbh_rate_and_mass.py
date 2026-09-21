"""
Extract the two GWTC-5.0 population results used in GROWL onboarding tutorial 3 -- the BBH merger
rate density R(z) and the primary-mass distribution dR/dm1 at z = 0.2 -- from the full LVK
data release into one small HDF5 file that can live in git.

Source: "GWTC-5.0: Population Properties of Merging Compact Binaries", data release
        https://doi.org/10.5281/zenodo.20292639  (popsummary_files.tar.gz, 25 GB compressed)

Files read (paths inside popsummary_files/):
  1. gwtc5_updated_default_mmax_mass_TwoPeakBrokenPowerLawSmoothedMassDistribution_redshift_PowerLawRedshift
       _magnitude_iid_spin_magnitude_gaussian_tilt_iid_spin_orientation_popsummary_result.h5
       -> the "Default BBH" parametric model (Broken Power Law + Two Peaks masses, power-law redshift).
          Used by the paper's Fig. 2 (m1) and Fig. 5 (R(z)).   [727 MB]
  2. may26_datarelease/z_varcut1/z_varcut1_popsummary.h5
       -> PixelPop (non-parametric) R(z).  Used by Fig. 5.        [151 MB]
  3. may26_datarelease/m1m2_varcut1/m1m2_varcut1_popsummary.h5
       -> PixelPop (non-parametric) m1 distribution. Used by Fig. 2.  [295 MB]
  4. figure_scripts/one_percent_pe_limits.json
       -> 1st/99th percentiles of the events' posterior samples (where the data constrain the fit).

How the quantities are built (copied from the paper's figure_scripts/make_fig_2.ipynb and make_fig_5.ipynb):
  * R(z):           posterior/rates_on_grids/redshift/rates, one row per posterior draw   [Gpc^-3 yr^-1]
  * dR/dm1 at z=0.2: posterior/rates_on_grids/mass_1/rates (evaluated at z = 0) x (1 + 0.2)^lamb, with
                     lamb the redshift-evolution hyperparameter sample of the same draw   [Gpc^-3 yr^-1 Msun^-1]
                     (PixelPop tabulates dR/dln(m1) on a log grid; we divide by m1 to get dR/dm1.)

Stored for each: the grid ('positions'), the 5/50/95 % quantiles over posterior draws ('quantiles',
rows in that order), and a thinned set of posterior draws ('draws', float32) for anyone who wants
other quantiles.  Run from this directory:

    python extract_gwtc5_bbh_rate_and_mass.py  /path/to/GWTC-5   [output.h5]
"""
import json
import os
import sys

import h5py
import numpy as np

SRC = sys.argv[1] if len(sys.argv) > 1 else "../GWtC-5"
OUT = sys.argv[2] if len(sys.argv) > 2 else "gwtc5_bbh_rate_and_mass.h5"
Z_EVAL = 0.2
N_DRAWS_KEEP = 500          # thinned posterior draws to store (the quantiles use all draws)
Q = [0.05, 0.50, 0.95]

PS = os.path.join(SRC, "popsummary_files")
DEFAULT_BBH = os.path.join(PS, "gwtc5_updated_default_mmax_mass_TwoPeakBrokenPowerLawSmoothedMassDistribution_"
                               "redshift_PowerLawRedshift_magnitude_iid_spin_magnitude_gaussian_tilt_iid_spin_"
                               "orientation_popsummary_result.h5")
PIXELPOP_Z = os.path.join(PS, "may26_datarelease/z_varcut1/z_varcut1_popsummary.h5")
PIXELPOP_M1 = os.path.join(PS, "may26_datarelease/m1m2_varcut1/m1m2_varcut1_popsummary.h5")
PE_LIMITS = os.path.join(SRC, "figure_scripts/one_percent_pe_limits.json")


def hyper(f, name):
    """One hyperparameter's posterior samples, looked up by name in the popsummary attrs."""
    names = [s.strip("'") for s in str(f.attrs["hyperparameters"]).strip("[]").split()]
    return f["posterior/hyperparameter_samples"][:, names.index(name)]


def store(grp, positions, draws, description, unit):
    rng = np.random.default_rng(0)
    keep = np.sort(rng.choice(len(draws), size=min(N_DRAWS_KEEP, len(draws)), replace=False))
    grp.create_dataset("positions", data=np.asarray(positions, dtype=np.float64))
    grp.create_dataset("quantiles", data=np.quantile(draws, Q, axis=0))
    grp.create_dataset("draws", data=np.asarray(draws[keep], dtype=np.float32), compression="gzip")
    grp.attrs["description"] = description
    grp.attrs["unit"] = unit
    grp.attrs["quantile_levels"] = Q
    grp.attrs["n_draws_total"] = len(draws)


with h5py.File(OUT, "w") as out:
    out.attrs["source"] = "GWTC-5.0: Population Properties of Merging Compact Binaries, https://doi.org/10.5281/zenodo.20292639"
    out.attrs["extracted_by"] = os.path.basename(__file__)
    out.attrs["z_eval_for_mass_distribution"] = Z_EVAL
    with open(PE_LIMITS) as fh:
        lim = json.load(fh)
    out.attrs["pe_redshift_1st_99th_percentile"] = [lim["redshift"]["1st"], lim["redshift"]["99th"]]
    out.attrs["pe_mass_1_source_1st_99th_percentile"] = [lim["mass_1_source"]["1st"], lim["mass_1_source"]["99th"]]

    # ---- Default BBH (parametric) --------------------------------------------------------------
    with h5py.File(DEFAULT_BBH) as f:
        g = out.create_group("default_bbh")
        g.attrs["source_file"] = os.path.relpath(DEFAULT_BBH, SRC)
        g.attrs["model"] = str(f.attrs["model_names"])
        z = f["posterior/rates_on_grids/redshift/positions"][0]
        Rz = f["posterior/rates_on_grids/redshift/rates"][()]
        store(g.create_group("rate_vs_redshift"), z, Rz, "BBH merger rate density R(z)", "Gpc^-3 yr^-1")
        m1 = f["posterior/rates_on_grids/mass_1/positions"][0]
        lamb = hyper(f, "lamb")
        dRdm1 = f["posterior/rates_on_grids/mass_1/rates"][()] * (1 + Z_EVAL) ** lamb[:, None]
        store(g.create_group("dR_dm1_z0.2"), m1, dRdm1, f"dR/dm1 at z = {Z_EVAL} (z=0 rate x (1+z)^lamb)", "Gpc^-3 yr^-1 Msun^-1")
        g.create_dataset("lamb", data=lamb.astype(np.float32))
        g["lamb"].attrs["description"] = "redshift-evolution hyperparameter: R(z) ~ (1+z)^lamb"

    # ---- PixelPop (non-parametric) -------------------------------------------------------------
    with h5py.File(PIXELPOP_Z) as f:
        g = out.create_group("pixelpop")
        g.attrs["source_file_redshift"] = os.path.relpath(PIXELPOP_Z, SRC)
        z = f["posterior/rates_on_grids/redshift/positions"][0]
        Rz = f["posterior/rates_on_grids/redshift/rates"][()].astype(np.float64)
        store(g.create_group("rate_vs_redshift"), z, Rz, "BBH merger rate density R(z), PixelPop", "Gpc^-3 yr^-1")
    with h5py.File(PIXELPOP_M1) as f:
        g.attrs["source_file_mass"] = os.path.relpath(PIXELPOP_M1, SRC)
        lnm1 = f["posterior/rates_on_grids/log_mass_1/positions"][0].astype(np.float64)
        m1 = np.exp(lnm1)
        lamb = hyper(f, "lamb").astype(np.float64)
        dRdlnm1 = f["posterior/rates_on_grids/log_mass_1/rates"][()].astype(np.float64)
        dRdm1 = dRdlnm1 / m1[None, :] * (1 + Z_EVAL) ** lamb[:, None]
        store(g.create_group("dR_dm1_z0.2"), m1, dRdm1, f"dR/dm1 at z = {Z_EVAL}, PixelPop (dR/dln m1 / m1 x (1+z)^lamb)", "Gpc^-3 yr^-1 Msun^-1")

print(f"wrote {OUT}: {os.path.getsize(OUT)/1e6:.1f} MB")
with h5py.File(OUT) as f:
    for name in ["default_bbh/rate_vs_redshift", "pixelpop/rate_vs_redshift", "default_bbh/dR_dm1_z0.2", "pixelpop/dR_dm1_z0.2"]:
        pos, q = f[name]["positions"][()], f[name]["quantiles"][()]
        i = np.argmin(np.abs(pos - (Z_EVAL if "redshift" in name else 35.0)))
        print(f"  {name:32s} grid {pos.min():.3g}..{pos.max():.3g} ({len(pos)} pts);  at {pos[i]:.3g}: "
              f"5/50/95% = {q[0,i]:.3g} / {q[1,i]:.3g} / {q[2,i]:.3g}")
