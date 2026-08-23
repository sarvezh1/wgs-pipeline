process SAMTOOLS_MARKDUP_QC {
    label 'process_medium'
    tag "${meta.sample_id}:${meta.role}"

    publishDir "${params.outdir ?: 'results'}/bam_qc", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(sorted_bam), path(sorted_bai)

    output:
    tuple val(meta), path("${meta.sample_id}.markdup.bam"), path("${meta.sample_id}.markdup.bam.bai"), path("${meta.sample_id}.markdup.stats.txt"), emit: marked
    tuple val(meta), path("${meta.sample_id}.flagstat.txt"), path("${meta.sample_id}.samtools_stats.txt"), emit: qc

    script:
    """
    set -euo pipefail

    samtools quickcheck '${sorted_bam}'
    samtools collate \
        -@ ${task.cpus} \
        -u \
        -o '${meta.sample_id}.namecollated.bam' \
        '${sorted_bam}'
    samtools fixmate \
        -@ ${task.cpus} \
        -m \
        '${meta.sample_id}.namecollated.bam' \
        '${meta.sample_id}.fixmate.bam'
    samtools sort \
        -@ ${task.cpus} \
        -m 256M \
        -o '${meta.sample_id}.fixmate.sorted.bam' \
        '${meta.sample_id}.fixmate.bam'
    samtools markdup \
        -@ ${task.cpus} \
        -s \
        -f '${meta.sample_id}.markdup.stats.txt' \
        '${meta.sample_id}.fixmate.sorted.bam' \
        '${meta.sample_id}.markdup.bam'
    samtools index \
        -@ ${task.cpus} \
        '${meta.sample_id}.markdup.bam' \
        '${meta.sample_id}.markdup.bam.bai'
    samtools quickcheck '${meta.sample_id}.markdup.bam'
    samtools flagstat '${meta.sample_id}.markdup.bam' > '${meta.sample_id}.flagstat.txt'
    samtools stats '${meta.sample_id}.markdup.bam' > '${meta.sample_id}.samtools_stats.txt'
    """
}
