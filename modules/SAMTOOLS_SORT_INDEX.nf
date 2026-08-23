process SAMTOOLS_SORT_INDEX {
    label 'process_medium'
    tag "${meta.sample_id}:${meta.role}"

    publishDir "${params.outdir ?: 'results'}/alignment", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(sam), path(alignment_log)

    output:
    tuple val(meta), path("${meta.sample_id}.sorted.bam"), path("${meta.sample_id}.sorted.bam.bai"), emit: aligned
    tuple val(meta), path("${meta.sample_id}.alignment.log"), emit: logs

    script:
    """
    set -euo pipefail
    samtools sort \
        -@ ${task.cpus} \
        -m 256M \
        -o '${meta.sample_id}.sorted.bam' \
        '${sam}'
    samtools index \
        -@ ${task.cpus} \
        '${meta.sample_id}.sorted.bam' \
        '${meta.sample_id}.sorted.bam.bai'
    samtools flagstat '${meta.sample_id}.sorted.bam' >> '${meta.sample_id}.alignment.log'
    samtools stats '${meta.sample_id}.sorted.bam' >> '${meta.sample_id}.alignment.log'
    """
}
