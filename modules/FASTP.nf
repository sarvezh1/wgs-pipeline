process FASTP {
    label 'process_medium'
    tag "${meta.sample_id}:${meta.role}"

    publishDir "${params.outdir ?: 'results'}/fastp", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_3'

    input:
    tuple val(meta), path(read1), path(read2)
    val validation_passed

    output:
    tuple val(meta), path("${meta.sample_id}_R1.fastp.fastq.gz"), path("${meta.sample_id}_R2.fastp.fastq.gz"), path("${meta.sample_id}_fastp.html"), path("${meta.sample_id}_fastp.json"), emit: cleaned

    script:
    """
    fastp \\
        --thread ${task.cpus} \\
        --in1 '${read1}' \\
        --in2 '${read2}' \\
        --out1 '${meta.sample_id}_R1.fastp.fastq.gz' \\
        --out2 '${meta.sample_id}_R2.fastp.fastq.gz' \\
        --html '${meta.sample_id}_fastp.html' \\
        --json '${meta.sample_id}_fastp.json' \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred 20 \\
        --unqualified_percent_limit 40 \\
        --length_required 15
    """
}
