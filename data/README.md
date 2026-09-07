# Data directory

Only the following supplied input files are required:

```text
metadata/sample_map
metadata/sample_info.xlsx
metadata/fungi.profile
metadata/fungi.profile.genus
metadata/bacteria_species.profile
```

All other files under `data/metadata/`, `data/processed/`, and `data/random_forest/` are generated files. They should not be supplied manually for a clean run. Generate them with:

```bash
Rscript code/01_prepare_inputs.R
```

The preparation script standardizes follow-up sample IDs, creates baseline and all-sample metadata tables, extracts clinical metadata, copies the supplied abundance profiles into standardized processed filenames, and generates random-forest feature matrices.
