# Details

## Naming requirements

### General workflow

It is advised not to mix ChIP-seq experiments using different anybodies in the same `snakemake` workflow. In other words, every ChIP-seq experiment with a specific antibody should get its own `snakemake` workflow.

### Samples

`samples` must follow the following naming convention: `{source}[_{condition}]_{replicate}`, where

- `{source}` is usually either `"IP"` or `"input"` (or `"IgG"`)
- `{condition}` can refer to a treatment condition, or a particular type of cell perturbation
- `{replicate}` denotes the biological replicate, and must follow the naming convention `"rep1"`, `"rep2"` and so on.

Valid names are for example: `"IP_doxorubicin_rep1"`, `"IP_vehicle_rep2"`, `"input_rep1"`.
