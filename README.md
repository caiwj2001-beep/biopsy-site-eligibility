# biopsy-site-eligibility

Analysis code for the manuscript:

> **Where the biopsy is taken changes treatment eligibility: within-patient discordance between primary tumours and metastases** (under review, *Annals of Oncology*, 2026)

The study asks a question that guidelines leave open: when several lesions can be biopsied, does the choice change whether a patient is eligible for a targeted therapy? It quantifies that discordance in 2,278 patients with paired primary and metastatic sequencing, and sets it against a **within-organ noise floor** — the discordance observed between two metastases sampled from the same organ in the same patient.

## Data sources (public, de-identified)

| Dataset | Access | Use |
|---|---|---|
| MSK-IMPACT clinical sequencing cohort (`msk_impact_50k_2026`) | cBioPortal, https://www.cbioportal.org/study/summary?id=msk_impact_50k_2026 | main analysis (54,331 samples, 48,179 patients) |
| AURORA US Network (`brca_aurora_2023`) | cBioPortal | independent directional replication (48 pairs) |

No individual-level identifiable data are redistributed here. Download both studies from cBioPortal into a local directory before running.

## Before running

Each script declares its input directory at the top. Set it to your local copy of `msk_impact_50k_2026`:

```python
BASE = r"<path to>/msk_impact_50k_2026"
```

Required files from that directory: `data_clinical_sample.csv`, `data_mutations.csv`, `data_cna.txt`. `replication_aurora.py` additionally reads the AURORA study files.

## What each script does

| Script | Purpose | Output |
|---|---|---|
| `analysis_main.py` | Eligibility classification (30 genes, ESCAT-based list), primary–metastasis discordance, within-organ noise floor, organ-level excess, decision impact, TMB/MSI, Cox model | `results/p1_*.csv`, `p2_*.csv`, `p4_*.csv`, `os_cox.csv`, `tmb_discordance.csv`, `msi_discordance.csv` |
| `standardisation.py` | Primary endpoint after direct standardisation to the tumour-type distribution | `results/p1_standardised_excess.csv` |
| `layer_controls.py` | Three-layer comparison (primary vs metastasis; metastasis vs metastasis across organs; within organ) and amplification decomposition | `results/layer_controls.csv`, `amplification_decomposition.csv` |
| `layer_controls_restricted.py` | Same three layers, restricted to patients contributing to the paired cohort | `results/layer_controls_restricted.csv` |
| `power_boundary.py` | Precision required for organ-specific testing rules | `results/p4_power_boundary.csv` |
| `replication_aurora.py` | AURORA replication of the cross-lesion rate | `results/replication_aurora.csv` |
| `export_figure_data.py` | Export of the exact numbers plotted in Figures 1–5 | `results/os_km_input.csv` and figure inputs |

## Key numbers reproduced by this code

| Quantity | Value |
|---|---|
| Primary–metastasis discordance | 19.1% (435/2,278) |
| Same-organ replicate discordance (noise floor) | 11.2% (62/555) |
| Excess, unstandardised | 7.9 percentage points (95% CI 4.8–11.0) |
| Excess, standardised to tumour type | 6.9 percentage points (95% CI 3.3–10.1) |
| Metastasis vs metastasis, different organs | 15.3% (152/992) |
| Loss if only the primary were tested | 12.5% (284/2,278), 99 of them amplification-only |
| Discordance vs overall survival | adjusted HR 0.97 (0.82–1.14) |

## Environment

Python 3.12 with `pandas`, `numpy`, `statsmodels`, `lifelines`; R 4.6.1 with `survival` and `ggplot2` for figures (figure scripts are not included in this repository).

## Notes and limitations

- Eligibility is defined at **gene level**: any non-synonymous variant or high-level amplification (discretised copy number >= 2) in a 30-gene panel. Variant-level pathogenicity annotation was not applied, so the rates are upper bounds.
- The public release contains **no treatment records and no specimen dates**; the survival analysis is therefore associative and is reported as a falsification test only.
- Structural-variant (fusion) calls were not available from the release used here; ALK, ROS1, RET and NTRK fusions are under-ascertained.

## Citation

If you use this code, please cite the manuscript above. A DOI for the code release will be added once the paper is accepted.

## Contact

Wen-Jie Cai, Department of Radiation Oncology, Quanzhou First Hospital Affiliated to Fujian Medical University — cwj@qzdyyy.com

## License

MIT (see `LICENSE`).
