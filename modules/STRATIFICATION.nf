process STRATIFIED_EXACT_BENCHMARK {
    label 'process_low'
    tag "${meta.patient_id ?: meta.sample_id}:${caller_name}:stratified-benchmark"
    publishDir "${params.outdir ?: 'results'}/benchmark/stratified", mode: 'copy', overwrite: true
    container 'python:3.11@sha256:d0199e2a90bf7a206a485b115323a75bc946f30b463d704c5435a454aca084dd'

    input:
    tuple val(meta), val(caller_name), path(truth_alleles), path(query_alleles), path(benchmark_regions), path(stratification_resources), val(stratification_tsv_name), val(analysis_type), val(benchmark_name), val(benchmark_version), val(benchmark_assembly)

    output:
    tuple val(meta), val(caller_name), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark_stratified_summary.tsv"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark_stratified_errors.tsv"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark_stratification_provenance.tsv"), emit: results

    script:
    def output_id = meta.patient_id ?: meta.sample_id
    """
    set -euo pipefail
    awk -F '\t' 'BEGIN { OFS="\t" } NF && !/^#/ { n=split(\$2, parts, "/"); print \$1, parts[n] }' \
        '${stratification_tsv_name}' > staged.stratifications.tsv
    stratify_exact_benchmark.py \
        --manifest staged.stratifications.tsv --benchmark-regions '${benchmark_regions}' \
        --truth-alleles '${truth_alleles}' --query-alleles '${query_alleles}' \
        --output '${output_id}.${caller_name}.benchmark_stratified_summary.tsv' \
        --errors '${output_id}.${caller_name}.benchmark_stratified_errors.tsv' \
        --provenance '${output_id}.${caller_name}.benchmark_stratification_provenance.tsv' \
        --sample '${output_id}' --analysis '${analysis_type}' --caller '${caller_name}' \
        --benchmark-name '${benchmark_name}' --benchmark-version '${benchmark_version}' \
        --engine normalized_exact_allele --engine-version '1.0 (bcftools 1.20)' --assembly '${benchmark_assembly}'
    """
}

process STRATIFIED_GERMLINE_COMPARISON {
    label 'process_low'
    tag "${meta.comparison_id}:germline:stratified-comparison"
    publishDir "${params.outdir ?: 'results'}/comparison/stratified", mode: 'copy', overwrite: true
    container 'python:3.11@sha256:d0199e2a90bf7a206a485b115323a75bc946f30b463d704c5435a454aca084dd'

    input:
    tuple val(meta), path(membership), path(discordant), path(genotype_details), path(benchmark_regions), path(stratification_resources), val(stratification_tsv_name), val(benchmark_name), val(benchmark_version)

    output:
    tuple val(meta), path("${meta.comparison_id}.stratified_concordance_summary.tsv"), path("${meta.comparison_id}.stratified_jaccard.tsv"), path("${meta.comparison_id}.stratified_discordant_variants.tsv"), path("${meta.comparison_id}.stratified_comparison_provenance.tsv"), emit: results
    path "${meta.comparison_id}.stratified_genotype_concordance.tsv", optional: true, emit: genotypes

    script:
    def output_id = meta.comparison_id
    """
    set -euo pipefail
    awk -F '\t' 'BEGIN { OFS="\t" } NF && !/^#/ { n=split(\$2, parts, "/"); print \$1, parts[n] }' \
        '${stratification_tsv_name}' > staged.stratifications.tsv
    stratify_comparison.py \
        --manifest staged.stratifications.tsv --benchmark-regions '${benchmark_regions}' \
        --membership '${membership}' --discordant '${discordant}' --output-prefix '${output_id}' \
        --sample '${output_id}' --analysis germline --benchmark-name '${benchmark_name}' \
        --benchmark-version '${benchmark_version}' --genotype-details '${genotype_details}'
    """
}

process STRATIFIED_SOMATIC_COMPARISON {
    label 'process_low'
    tag "${meta.comparison_id}:somatic:stratified-comparison"
    publishDir "${params.outdir ?: 'results'}/comparison/stratified", mode: 'copy', overwrite: true
    container 'python:3.11@sha256:d0199e2a90bf7a206a485b115323a75bc946f30b463d704c5435a454aca084dd'

    input:
    tuple val(meta), path(membership), path(discordant), path(benchmark_regions), path(stratification_resources), val(stratification_tsv_name), val(benchmark_name), val(benchmark_version)

    output:
    tuple val(meta), path("${meta.comparison_id}.stratified_concordance_summary.tsv"), path("${meta.comparison_id}.stratified_jaccard.tsv"), path("${meta.comparison_id}.stratified_discordant_variants.tsv"), path("${meta.comparison_id}.stratified_comparison_provenance.tsv"), emit: results

    script:
    def output_id = meta.comparison_id
    """
    set -euo pipefail
    awk -F '\t' 'BEGIN { OFS="\t" } NF && !/^#/ { n=split(\$2, parts, "/"); print \$1, parts[n] }' \
        '${stratification_tsv_name}' > staged.stratifications.tsv
    stratify_comparison.py \
        --manifest staged.stratifications.tsv --benchmark-regions '${benchmark_regions}' \
        --membership '${membership}' --discordant '${discordant}' --output-prefix '${output_id}' \
        --sample '${output_id}' --analysis somatic --benchmark-name '${benchmark_name}' \
        --benchmark-version '${benchmark_version}'
    """
}
