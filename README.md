# deeptools-nextflow
A small collection of Nextflow workflows for reproducible deepTools-based signal processing and visualization.

## Workflows

### `signal_matrix_profile`

Runs a small deepTools workflow for generating normalized signal tracks, matrices, heatmaps, and profile plots.

The workflow can start from BAM files, existing BigWig files, or existing computeMatrix output files, depending on which steps are enabled in the config.

Supported steps:

```text
BAM -> bamCoverage -> BigWig
BigWig + BED -> computeMatrix
matrix -> plotHeatmap
matrix -> plotProfile
```

Example configuration

The included config file is an example and uses deepTools v3.5.6 through Apptainer.

As provided, it starts from a BAM file, creates a CPM-normalized BigWig with `bamCoverage`, creates a TSS-centered matrix with `computeMatrix`, and plots the result with `plotHeatmap`. `plotProfile` is available but turned off by default.

`computeMatrix` is set to `reference-point` mode at gene TSS positions, using 2 kb upstream and 2 kb downstream of each TSS.

Input files, output folders, deepTools options, Apptainer settings, and Slurm resources are all set in the config file.
