process CALLSET_COMPARISON {
    label 'process_low'
    tag "${meta.comparison_id}:${analysis_type}:comparison"

    publishDir "${params.outdir ?: 'results'}/comparison", mode: 'copy', overwrite: true, pattern: '*.tsv'
    publishDir "${params.outdir ?: 'results'}/comparison/vcfs", mode: 'copy', overwrite: true, pattern: '*.comparison.*.vcf.gz*'
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    tuple val(meta), val(callers), path(normalized_vcfs), path(normalized_indexes), path(reference_fa), path(reference_fai), path(reference_dict), path(comparison_resources, stageAs: 'comparison_resources/*'), val(analysis_type), val(assembly), val(bqsr_enabled), val(benchmark_enabled), val(benchmark_name), val(benchmark_version), val(truth_name), val(truth_index_name), val(regions_name)

    output:
    tuple val(meta), path("${meta.comparison_id}.concordance_summary.tsv"), path("${meta.comparison_id}.pairwise_jaccard.tsv"), path("${meta.comparison_id}.agreement_matrix.tsv"), path("${meta.comparison_id}.membership.tsv"), path("${meta.comparison_id}.discordant_variants.tsv"), path("${meta.comparison_id}.comparison_provenance.tsv"), path("${meta.comparison_id}.comparison.*.vcf.gz*"), emit: results
    path "${meta.comparison_id}.genotype_concordance.tsv", optional: true, emit: genotype_summary
    tuple val(meta), path("${meta.comparison_id}.genotype_details.tsv"), optional: true, emit: genotype_details

    script:
    def output_id = meta.comparison_id
    def expected_sample_arg = analysis_type == 'germline' ? "--expected-sample '${meta.sample_id}'" : ''
    def vcf_args = callers.indices.collect { index ->
        "--vcf '${callers[index]}=${normalized_vcfs[index]}' --index '${callers[index]}=${normalized_indexes[index]}'"
    }.join(' ')
    def caller_versions = analysis_type == 'germline'
        ? "--caller-version 'haplotypecaller=GATK 4.6.2.0' --caller-version 'deepvariant=1.10.0'"
        : "--caller-version 'mutect2=GATK 4.6.2.0' --caller-version 'strelka2=2.9.10' --caller-version 'deepsomatic=1.10.0'"
    def truth_setup = benchmark_enabled.toString() == 'true' ? """
    test -s 'comparison_resources/${truth_name}' || { echo 'CALLSET_COMPARISON ERROR: benchmark truth is missing' >&2; exit 1; }
    test -s 'comparison_resources/${truth_index_name}' || { echo 'CALLSET_COMPARISON ERROR: benchmark truth index is missing' >&2; exit 1; }
    test -s 'comparison_resources/${regions_name}' || { echo 'CALLSET_COMPARISON ERROR: benchmark regions are missing' >&2; exit 1; }
    bcftools norm \\
        --fasta-ref '${reference_fa}' \\
        --multiallelics -any \\
        --output-type z \\
        --output truth.normalized.vcf.gz \\
        'comparison_resources/${truth_name}'
    bcftools index --tbi truth.normalized.vcf.gz
    bcftools view \\
        --regions-file 'comparison_resources/${regions_name}' \\
        --output-type z \\
        --output truth.benchmark-region.vcf.gz \\
        truth.normalized.vcf.gz
    bcftools index --tbi truth.benchmark-region.vcf.gz
    truth_args="--truth truth.benchmark-region.vcf.gz --regions comparison_resources/${regions_name}"
    """ : "truth_args=''"
    """
    set -euo pipefail

    ${truth_setup}

    compare_callsets.pl \\
        --analysis '${analysis_type}' \\
        --sample '${output_id}' \\
        ${expected_sample_arg} \\
        --assembly '${assembly}' \\
        --bqsr '${bqsr_enabled}' \\
        --benchmark-enabled '${benchmark_enabled}' \\
        --benchmark-name '${benchmark_name}' \\
        --benchmark-version '${benchmark_version}' \\
        --output-prefix '${output_id}' \\
        ${caller_versions} \\
        ${vcf_args} \\
        \${truth_args}

    for vcf in '${output_id}'.comparison.*.vcf; do
        bgzip -c "\${vcf}" > "\${vcf}.gz"
        tabix -p vcf "\${vcf}.gz"
        rm "\${vcf}"
    done

    test -s '${output_id}.concordance_summary.tsv'
    test -s '${output_id}.pairwise_jaccard.tsv'
    test -s '${output_id}.agreement_matrix.tsv'
    test -s '${output_id}.membership.tsv'
    test -s '${output_id}.discordant_variants.tsv'
    test -s '${output_id}.comparison_provenance.tsv'
    """
}
