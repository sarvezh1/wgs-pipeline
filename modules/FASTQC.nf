process FASTQC {
    label 'process_low'
    tag "${meta.sample_id}:${meta.role}:${fastq.simpleName}"

    publishDir "${params.outdir ?: 'results'}/fastqc", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(fastq)
    val validation_passed

    output:
    tuple val(meta), path('*_fastqc.html'), path('*_fastqc.zip'), emit: reports

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        '${fastq}'
    """
}
