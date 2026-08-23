process DEEPVARIANT {
    label 'process_high'
    tag "${meta.sample_id}:${meta.role}"

    publishDir "${params.outdir ?: 'results'}/variants", mode: 'copy', overwrite: true
    container 'google/deepvariant:1.10.0'

    input:
    tuple val(meta), val(caller_name), path(caller_bam), path(caller_bai), path(reference_fa), path(reference_fai), path(reference_dict), val(interval)

    output:
    tuple val(meta), val(caller_name), path("${meta.sample_id}.deepvariant.vcf.gz"), path("${meta.sample_id}.deepvariant.vcf.gz.tbi"), emit: calls

    script:
    def interval_arg = interval ? "--regions='${interval}'" : ''
    """
    set -euo pipefail

    /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref='${reference_fa}' \
        --reads='${caller_bam}' \
        --output_vcf='${meta.sample_id}.deepvariant.vcf.gz' \
        ${interval_arg} \
        --num_shards=${task.cpus} \
        --intermediate_results_dir='deepvariant_intermediate' \
        --logging_dir='deepvariant_logging'

    test -s '${meta.sample_id}.deepvariant.vcf.gz'
    test -s '${meta.sample_id}.deepvariant.vcf.gz.tbi'
    """
}
