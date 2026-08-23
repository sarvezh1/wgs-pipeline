process GATK_BQSR {
    label 'process_medium'
    tag "${meta.sample_id}:${meta.role}"

    publishDir "${params.outdir ?: 'results'}/bqsr", mode: 'copy', overwrite: true
    container 'broadinstitute/gatk:4.6.2.0'

    input:
    tuple val(meta), path(markdup_bam), path(markdup_bai), path(markdup_stats), path(reference_fa), path(reference_fai), path(reference_dict), path(known_sites_resources), val(interval)

    output:
    tuple val(meta), path("${meta.sample_id}.recal.bam"), path("${meta.sample_id}.recal.table"), emit: recalibrated

    script:
    def interval_arg = interval ? "--intervals '${interval}'" : ''
    def known_sites_args = known_sites_resources
        .findAll { resource -> !resource.name.endsWith('.tbi') }
        .collect { resource -> "--known-sites '${resource}'" }
        .join(' ')
    """
    set -euo pipefail

    gatk --java-options '-Xmx768m' BaseRecalibrator \
        --reference '${reference_fa}' \
        --input '${markdup_bam}' \
        ${known_sites_args} \
        ${interval_arg} \
        --output '${meta.sample_id}.recal.table'

    gatk --java-options '-Xmx768m' ApplyBQSR \
        --reference '${reference_fa}' \
        --input '${markdup_bam}' \
        --bqsr-recal-file '${meta.sample_id}.recal.table' \
        --output '${meta.sample_id}.recal.bam'
    """
}
