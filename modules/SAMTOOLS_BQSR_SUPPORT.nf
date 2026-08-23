process SAMTOOLS_BQSR_INDEX {
    label 'process_low'
    tag "${meta.sample_id}:${meta.role}"

    publishDir "${params.outdir ?: 'results'}/bqsr", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(recal_bam), path(recal_table)

    output:
    tuple val(meta), path("${meta.sample_id}.recal.bam"), path("${meta.sample_id}.recal.bam.bai"), path("${meta.sample_id}.recal.table"), emit: caller_ready

    script:
    """
    set -euo pipefail
    samtools index -@ ${task.cpus} '${recal_bam}' '${meta.sample_id}.recal.bam.bai'
    samtools quickcheck '${recal_bam}'
    samtools view -H '${recal_bam}' | grep -q '^@SQ'
    """
}
