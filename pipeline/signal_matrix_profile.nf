nextflow.enable.dsl=2

/*
  deepTools pipeline for signal visualization 
  over specified regions.

  BAM -> bamCoverage -> BigWig
  BigWig + BED -> computeMatrix
  matrix -> plotHeatmap / plotProfile
*/

def cleanName(x) {
  x.toString().replaceAll(/[^A-Za-z0-9_.-]/, '_')
}

def computeMatrixArgs() {
  def mode = params.compute_matrix_mode.toString()

  if (mode == 'scale-regions') {
    return "${params.computematrix_common_args} ${params.computematrix_scale_regions_args}"
  }

  if (mode == 'reference-point') {
    return "${params.computematrix_common_args} ${params.computematrix_reference_point_args}"
  }

  error "params.compute_matrix_mode must be 'scale-regions' or 'reference-point'"
}

process BAMCOVERAGE {
  tag { sample }

  publishDir "${params.outdir}/bamCoverage", mode: 'copy', pattern: '*.bw', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), path(bam), path(index)

  output:
  tuple val(sample), path("${sample}.bw"), emit: bw
  path "${sample}.bamCoverage.err", emit: err

  script:
  """
  set -euo pipefail

  {
    bamCoverage \\
      -b ${bam} \\
      -o ${sample}.bw \\
      --outFileFormat bigwig \\
      --numberOfProcessors ${task.cpus} \\
      ${params.bamcoverage_args}
  } 2> ${sample}.bamCoverage.err
  """
}

process COMPUTEMATRIX {
  tag { "${sample}.${region_set}" }

  publishDir { "${params.outdir}/computeMatrix/${region_set}" }, mode: 'copy', pattern: '*.matrix.mat.gz', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), path(bigwig), val(region_set), path(regions)

  output:
  tuple val(sample), val(region_set), path("${sample}.${region_set}.matrix.mat.gz"), emit: matrix
  path "${sample}.${region_set}.computeMatrix.err", emit: err

  script:
  """
  set -euo pipefail

  {
    computeMatrix \\
      ${params.compute_matrix_mode} \\
      -S ${bigwig} \\
      -R ${regions} \\
      -o ${sample}.${region_set}.matrix.mat.gz \\
      ${computeMatrixArgs()}
  } 2> ${sample}.${region_set}.computeMatrix.err
  """
}

process PLOTHEATMAP {
  tag { "${sample}.${region_set}" }

  publishDir { "${params.outdir}/plotHeatmap/${region_set}" }, mode: 'copy', pattern: '*.png', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(region_set), path(matrix)

  output:
  tuple val(sample), val(region_set), path("${sample}.${region_set}.heatmap.png"), emit: heatmap
  path "${sample}.${region_set}.plotHeatmap.err", emit: err

  script:
  """
  set -euo pipefail

  {
    plotHeatmap \\
      -m ${matrix} \\
      -out ${sample}.${region_set}.heatmap.png \\
      ${params.plotheatmap_args}
  } 2> ${sample}.${region_set}.plotHeatmap.err
  """
}

process PLOTPROFILE {
  tag { "${sample}.${region_set}" }

  publishDir { "${params.outdir}/plotProfile/${region_set}" }, mode: 'copy', pattern: '*.png', overwrite: true
  publishDir "${params.outdir}/pipeline_info/logs", mode: 'copy', pattern: '*.err', overwrite: true

  input:
  tuple val(sample), val(region_set), path(matrix)

  output:
  tuple val(sample), val(region_set), path("${sample}.${region_set}.profile.png"), emit: profile
  path "${sample}.${region_set}.plotProfile.err", emit: err

  script:
  """
  set -euo pipefail

  {
    plotProfile \\
      -m ${matrix} \\
      -out ${sample}.${region_set}.profile.png \\
      ${params.plotprofile_args}
  } 2> ${sample}.${region_set}.plotProfile.err
  """
}

workflow {
  if (!(params.compute_matrix_mode.toString() in ['scale-regions', 'reference-point'])) {
    error "params.compute_matrix_mode must be 'scale-regions' or 'reference-point'"
  }

  if (params.run_bamcoverage) {
    bam_files = Channel
      .fromList(params.bam_files)
      .map { row ->
        tuple(
          cleanName(row.sample),
          file(row.bam),
          file(row.index)
        )
      }

    BAMCOVERAGE(bam_files)
    signal_files = BAMCOVERAGE.out.bw
  } else {
    signal_files = Channel
      .fromList(params.bigwig_files)
      .map { row ->
        tuple(
          cleanName(row.sample),
          file(row.bigwig)
        )
      }
  }

  if (params.run_computematrix) {
    region_files = Channel
      .fromList(params.region_files)
      .map { row ->
        tuple(
          cleanName(row.name),
          file(row.bed)
        )
      }

    compute_input = signal_files
      .combine(region_files)
      .map { sample, bigwig, region_set, regions ->
        tuple(sample, bigwig, region_set, regions)
      }

    COMPUTEMATRIX(compute_input)
    matrices = COMPUTEMATRIX.out.matrix
  } else {
    matrices = Channel
      .fromList(params.matrix_files)
      .map { row ->
        tuple(
          cleanName(row.sample),
          cleanName(row.region_set),
          file(row.matrix)
        )
      }
  }

  if (params.run_plotheatmap) {
    PLOTHEATMAP(matrices)
  }

  if (params.run_plotprofile) {
    PLOTPROFILE(matrices)
  }
}
