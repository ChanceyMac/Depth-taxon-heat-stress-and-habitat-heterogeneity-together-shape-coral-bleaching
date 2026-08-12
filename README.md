# Depth, taxon, heat stress, and habitat heterogeneity together shape coral bleaching

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17476109.svg)](https://doi.org/10.5281/zenodo.17476109)

Code and data supporting:

> **MacDonald C., Djakiman C., Masdar H., Bachman S., Green R., Fox H., Limmon G., Beger M.** (in review). *Depth, taxon, heat stress, and habitat heterogeneity together shape coral bleaching.* **Ecology and Evolution**.

**Corresponding author:** Chancey MacDonald — [chancey.macdonald@newcastle.ac.uk](mailto:chancey.macdonald@newcastle.ac.uk)

**Author affiliations**

1. School of Natural and Environmental Sciences, Newcastle University, Newcastle Upon Tyne, United Kingdom
2. School of Biology, University of Leeds, Leeds, United Kingdom
3. Maritime and Marine Science Center of Excellence, Pattimura University, Ambon, Indonesia
4. Faculty of Marine Sciences and Fisheries, Hasanuddin University, Makassar, Sulawesi Selatan, Indonesia
5. [C]Worthy, LLC, Boulder, Colorado, United States of America
6. Coral Reef Alliance, San Francisco, California, United States of America
7. Centre for Biodiversity and Conservation Science, School of the Environment, University of Queensland, Brisbane, Queensland, Australia

**Keywords:** Coral Bleaching · Depth Refuge · Habitat Heterogeneity · Heat Stress

---

## Abstract

The depth refugia hypothesis proposes that corals find refuge from heat-stress induced coral bleaching at depth. In addition, geomorphologically complex reefs may provide additional micro-refugia. However, evidence for pervasive community-level depth and micro refugia has been equivocal, necessitating a nuanced exploration of local coral bleaching responses. We assess the effects of depth and morphological habitat heterogeneity on coral bleaching. We include over 6000 coral colonies between two and eighteen meters of depth over a heat stress gradient from 1 to 6 degree heating weeks, and with varying site-level habitat heterogeneity. Our results demonstrate strong interactions in bleaching susceptibility along gradients of depth, heat stress, and habitat heterogeneity that were taxa-specific and varied with colony level bleaching severity. Community-level bleaching susceptibility decreased with depth, but severe bleaching of colonies (> 50 % of coral surface) increased with depth. The coral community as a whole was more susceptible to bleaching as cumulative heat stress exposure increased, but did not vary uniformly with large-scale habitat heterogeneity. Genus-level bleaching susceptibility mostly followed expectations, with *Seriatopora* and *Pocillopora* being the most susceptible genera overall, and *Porites* the least susceptible. However, *Acropora* corals bleached less than expected overall, and in unexpected patterns along interacting depth and heat stress gradients. Overall, our results demonstrate that bleaching susceptibility varies substantially with multiple interacting biological and environmental factors, and that the influence of these can change with the level of bleaching assessed. Consequently, depth effects on coral bleaching can be highly nuanced and care should be taken when investigating and reporting on the 'depth refugia' potential for coral bleaching and when interpreting community-level bleaching outcomes that are naïve to taxa and water depth.

---

## Overview

This repository contains the cleaned working data and analysis code used to model the probability of coral bleaching across depth, taxa, satellite-derived heat stress (Degree Heating Weeks), in-situ temperature, and satellite-derived habitat heterogeneity during the 2024 bleaching event in the Banda Sea and Southeast Sulawesi region, Indonesia.

All models are Bayesian generalised linear mixed models fit with [`brms`](https://paul-buerkner.github.io/brms/); posterior summaries, contrasts, and figures are produced from posterior draws using [`tidybayes`](https://mjskay.github.io/tidybayes/).

## Repository structure

```
.
├── README.md                                       # this file
├── LICENSE                                         # code + data licence(s)
├── .gitignore
├── R/
│   ├── SES_bleaching_library.R                     # package loader
│   └── SES_bleaching_analysis_annotated.R          # analysis script
├── Data/
│   └── bleaching_data_cleaned.csv                  # cleaned, joined working data
└── Outputs/                                        # figures + contrast tables written here
    └── .gitkeep
```

## Data

### `Data/bleaching_data_cleaned.csv`

The cleaned, colony-level working dataset used for all analyses in the manuscript. Colony-level in-water bleaching scores were joined with site-level benthic cover and remote-sensed habitat and heat-stress metrics.

Key columns:

| Column | Description |
|---|---|
| `Site_code` / `Site` | Site code and full location name |
| `Date` | Survey date |
| `Depth_m` | Sampling depth (m), numeric (2, 6, 12, 18) |
| `Rep` | Quadrat identifier within site × depth |
| `Taxa` | Coral genus / taxonomic grouping |
| `Bleach_cat` | Siebeck-style 6-point severity score (1 = unbleached, 6 = fully bleached/dead) |
| `Bleached` | Binary indicator of any bleaching (categories 2–6) |
| `Bleached_sev` | Binary indicator of severe bleaching (categories 4–6, > 50 % of colony) |
| `Temp` | In-situ water temperature at sampling (°C) |
| `CRW_DHW` | NOAA Coral Reef Watch cumulative Degree Heating Weeks at time of sampling |
| `HabHet_mean_150m`, `HabHet_max_150m`, `HabHet_mean_1000m`, `HabHet_max_1000m` | Satellite-derived habitat heterogeneity metrics (mean / max at 150 m and 1000 m radii) |
| `Coral_cover`, `Dead_coral`, `Dead_coral_recent`, `Algae`, `Soft_coral`, `Sponge`, `Rock`, `Rubble`, `Sand`, `Silt`, `TWS`, `Other` | Site-level mean benthic cover (proportions), derived from point-intercept transect data |

*Note:* CSV format does not preserve factor level ordering. The analysis script coerces `Taxa`, `Site`, and `Rep` to factors on read; downstream regressions therefore use the default (alphabetical) reference levels. All fitted probabilities, contrasts, and figures in the manuscript are invariant to this choice.

## Code

### `R/SES_bleaching_analysis_annotated.R`

The analysis script. Reads `Data/bleaching_data_cleaned.csv`, fits every model reported in the manuscript, produces Figures 2–6, and generates the supplementary contrast tables. Minimum working examples of the monotonic-depth sensitivity models are included.

### `R/SES_bleaching_library.R`

Loads the R packages required by the analysis script. Installs anything missing.

## Reproducing the analysis

### Requirements

- **R** ≥ 4.4.0
- **cmdstanr** with CmdStan installed (see [cmdstanr installation guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html))
- Packages loaded by `R/SES_bleaching_library.R`: `brms`, `tidybayes`, `tidyverse`, `ggplot2`, `ggridges`, `patchwork`, `HDInterval`, `parameters`, `modelr`, `paletteer`, and dependencies.

### Runtime

Each Bayesian GLMM takes 5 – 40 minutes on a modern multi-core machine depending on the model. The full script fits ~20 models; expect 4 – 8 hours end-to-end depending on hardware. Individual models can be re-run in isolation without refitting others.

### Full reproduction

```r
setwd("path/to/repo")
source("R/SES_bleaching_library.R")
source("R/SES_bleaching_analysis_annotated.R")
```

Figures and contrast tables are written to `Outputs/`.

## Session info

The analysis was run under:

```
R version 4.4.3 (2025-02-28)
Platform: aarch64-apple-darwin20 (macOS 26.3.1)

brms 2.21.0 · cmdstanr 0.9.0 · rstan 2.32.6 · StanHeaders 2.32.9
tidybayes 3.0.7 · bayestestR 0.15.2 · parameters 0.22.1 · modelr 0.1.11
tidyverse 2.0.0 · patchwork 1.2.0 · ggridges 0.5.6 · paletteer 1.6.0
```

A full `sessionInfo()` capture is preserved with the archived Zenodo release.

## Citation

If you use this code or data, please cite the manuscript:

> MacDonald, C., Djakiman, C., Masdar, H., Bachman, S., Green, R., Fox, H., Limmon, G., & Beger, M. (in review). Depth, taxon, heat stress, and habitat heterogeneity together shape coral bleaching. *Ecology and Evolution*.

and the archived repository:

> MacDonald, C., Djakiman, C., Masdar, H., Bachman, S., Green, R., Fox, H., Limmon, G., & Beger, M. (2026). *Data and code for 'Depth, taxon, heat stress, and habitat heterogeneity together shape coral bleaching'* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.17476109

## Licence

- **Code** (everything under `R/`): [MIT recommended]
- **Data** (everything under `Data/`): [CC-BY-4.0 recommended]

See [LICENSE](./LICENSE) for the chosen licence terms.

## Contributions

CM conceived the research, collected field data, led analysis and wrote the first draft of the manuscript. CD assisted with analysis, provided intellectual contributions and edited the manuscript. HM ran logistics for field work and edited the manuscript. SB developed the analysis pipeline for habitat heterogeneity measures, calculated the measures and edited the manuscript. RG calculated habitat heterogeneity measures and edited the manuscript. HF funded the research and edited the manuscript. GL funded the research and edited the manuscript. MB conceived and funded the research, provided intellectual contributions and edited the manuscript.

## Funding

This work was funded by UK aid from the UK government and by the International Development Research Centre, Ottawa, Canada as part of the Climate Adaptation and Resilience (CLARE) research programme, through the Climate REEFS project (MB, HF, GL). The views expressed herein do not necessarily represent those of the UK government, IDRC, or its Board of Governors. CM receives funding from a NERC Independent Research Fellowship. HF and RG were supported in part by grants from the Paul M Angell Family Foundation and the Builders Initiative.

## Acknowledgements

We thank Prandito Simanjuntak, Fajrin Rahayaan, and the crew of MV *Cakrawala Maritim* for logistical support for the field work; and Ben Charo for support with grant preparation.

## Contact & issues

For questions about the code or data, please open an issue on this repository. For scientific correspondence, contact [chancey.macdonald@newcastle.ac.uk](mailto:chancey.macdonald@newcastle.ac.uk).
