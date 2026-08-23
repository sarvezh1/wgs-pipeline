process DEEPSOMATIC {
    label 'process_high'
    tag "${meta.patient_id}:${meta.tumor_sample_id}/${meta.normal_sample_id}"

    publishDir "${params.outdir ?: 'results'}/somatic", mode: 'copy', overwrite: true
    container 'google/deepsomatic:1.10.0'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(reference_fa), path(reference_fai), path(reference_dict), val(expected_tumor_sample), val(expected_normal_sample), val(interval)

    output:
    tuple val(meta), val('deepsomatic'), path("${meta.patient_id}.deepsomatic.vcf.gz"), path("${meta.patient_id}.deepsomatic.vcf.gz.tbi"), emit: calls

    script:
    def interval_arg = interval ? "--regions='${interval}'" : ''
    """
    set -euo pipefail

    mkdir -p deepsomatic_intermediate deepsomatic_logging

    run_deepsomatic \
        --model_type=WGS \
        --ref='${reference_fa}' \
        --reads_tumor='${tumor_bam}' \
        --reads_normal='${normal_bam}' \
        --sample_name_tumor='${expected_tumor_sample}' \
        --sample_name_normal='${expected_normal_sample}' \
        --output_vcf='${meta.patient_id}.deepsomatic.vcf.gz' \
        ${interval_arg} \
        --num_shards=${task.cpus} \
        --intermediate_results_dir='deepsomatic_intermediate' \
        --logging_dir='deepsomatic_logging'

    test -s '${meta.patient_id}.deepsomatic.vcf.gz'
    test -s '${meta.patient_id}.deepsomatic.vcf.gz.tbi'
    """
}
