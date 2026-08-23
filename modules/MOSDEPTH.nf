process MOSDEPTH {
    label 'process_medium'
    tag "${meta.sample_id}:coverage"
    container 'quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_1'
    publishDir "${params.outdir ?: 'results'}/coverage/native", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam), path(bai), path(reference_fai), val(calling_interval), val(thresholds)

    output:
    tuple val(meta), path("${meta.sample_id}.mosdepth.summary.txt"), path("${meta.sample_id}.mosdepth.global.dist.txt"), path("${meta.sample_id}.mosdepth.region.dist.txt"), path("${meta.sample_id}.regions.bed.gz"), path("${meta.sample_id}.regions.bed.gz.csi"), path("${meta.sample_id}.thresholds.bed.gz"), path('coverage_regions.bed'), emit: coverage

    script:
    """
    set -euo pipefail
    if [[ -n '${calling_interval}' ]]; then
        region='${calling_interval}'
        contig="\${region%%:*}"
        range="\${region#*:}"
        start="\${range%-*}"
        end="\${range#*-}"
        printf '%s\t%s\t%s\tconfigured_calling_interval\n' "\$contig" "\$((start - 1))" "\$end" > coverage_regions.bed
    else
        awk 'BEGIN{OFS="\t"} {print \$1,0,\$2,\$1}' '${reference_fai}' > coverage_regions.bed
    fi
    mosdepth --threads ${task.cpus} --no-per-base --by coverage_regions.bed --thresholds '${thresholds}' '${meta.sample_id}' '${bam}'
    test -s '${meta.sample_id}.mosdepth.summary.txt'
    test -s '${meta.sample_id}.regions.bed.gz'
    test -s '${meta.sample_id}.thresholds.bed.gz'
    """
}

process COVERAGE_SUMMARY {
    label 'process_low'
    tag "${meta.sample_id}:coverage summary"
    container 'python:3.11@sha256:d0199e2a90bf7a206a485b115323a75bc946f30b463d704c5435a454aca084dd'
    publishDir "${params.outdir ?: 'results'}/coverage/summaries", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(summary), path(global_dist), path(region_dist), path(regions), path(regions_index), path(threshold_bed), path(coverage_regions), val(analysis_type), val(assembly)

    output:
    tuple val(meta), path('*.coverage_summary.tsv'), path('*.coverage_threshold_summary.tsv'), path('*.coverage_chromosome_summary.tsv'), emit: summaries

    script:
    """
    set -euo pipefail
    summarize_mosdepth.py --summary '${summary}' --regions '${regions}' --thresholds '${threshold_bed}' --sample '${meta.sample_id}' --patient '${meta.patient_id}' --role '${meta.role}' --analysis '${analysis_type}' --assembly '${assembly}' --prefix '${meta.sample_id}'
    """
}
