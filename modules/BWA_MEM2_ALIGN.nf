process BWA_MEM2_ALIGN {
    label 'process_high'
    tag "${meta.sample_id}:${meta.role}"

    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    input:
    tuple val(meta), path(read1), path(read2), path(reference), path(index_files), val(index_prefix)

    output:
    tuple val(meta), path("${meta.sample_id}.sam"), path("${meta.sample_id}.alignment.log"), emit: aligned

    script:
    """
    set -euo pipefail
    bwa-mem2 mem \
        -t ${task.cpus} \
        -R '@RG\\tID:${meta.sample_id}\\tSM:${meta.sample_id}\\tPL:ILLUMINA' \
        '${index_prefix}' '${read1}' '${read2}' \
        -o '${meta.sample_id}.sam' \
        2> '${meta.sample_id}.alignment.log'
    """
}
