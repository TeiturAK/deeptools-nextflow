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

