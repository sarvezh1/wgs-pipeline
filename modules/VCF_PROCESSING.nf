process VCF_PROCESSING {
    label 'process_low'
    tag "${meta.patient_id ?: meta.sample_id}:${caller_name}"

    publishDir "${params.outdir ?: 'results'}/vcf_processing", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    tuple val(meta), val(caller_name), path(native_vcf), path(native_vcf_index), path(reference_fa), path(reference_fai), path(reference_dict), val(analysis_type)

    output:
    tuple val(meta), val(caller_name), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.normalized.vcf.gz"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.normalized.vcf.gz.tbi"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.bcftools_stats.txt"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.variant_summary.tsv"), emit: processed

    script:
    def output_id = meta.patient_id ?: meta.sample_id
    """
    set -euo pipefail

    test -s '${native_vcf}'
    test -s '${native_vcf_index}'

    bcftools norm \
        --fasta-ref '${reference_fa}' \
        --multiallelics -any \
        --output-type z \
        --output '${output_id}.${caller_name}.normalized.vcf.gz' \
        '${native_vcf}'

    bcftools index \
        --tbi \
        --force \
        '${output_id}.${caller_name}.normalized.vcf.gz'

    bcftools stats \
        '${output_id}.${caller_name}.normalized.vcf.gz' \
        > '${output_id}.${caller_name}.bcftools_stats.txt'

    total_records=\$(bcftools view -H '${output_id}.${caller_name}.normalized.vcf.gz' | wc -l)
    pass_records=\$(bcftools view -f PASS -H '${output_id}.${caller_name}.normalized.vcf.gz' | wc -l)
    snvs=\$(bcftools view -v snps -H '${output_id}.${caller_name}.normalized.vcf.gz' | wc -l)
    indels=\$(bcftools view -v indels -H '${output_id}.${caller_name}.normalized.vcf.gz' | wc -l)
    multiallelic_records=\$(bcftools query -f '%ALT\\n' '${output_id}.${caller_name}.normalized.vcf.gz' | awk -F, 'NF > 1 { count++ } END { print count + 0 }')

    printf 'caller\tanalysis_type\tsample_or_patient\ttotal_records\tpass_records\tsnvs\tindels\tmultiallelic_records\n' \
        > '${output_id}.${caller_name}.variant_summary.tsv'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        '${caller_name}' \
        '${analysis_type}' \
        '${output_id}' \
        "\${total_records}" \
        "\${pass_records}" \
        "\${snvs}" \
        "\${indels}" \
        "\${multiallelic_records}" \
        >> '${output_id}.${caller_name}.variant_summary.tsv'

    bcftools view -Ov '${output_id}.${caller_name}.normalized.vcf.gz' > /dev/null
    test -s '${output_id}.${caller_name}.normalized.vcf.gz'
    test -s '${output_id}.${caller_name}.normalized.vcf.gz.tbi'
    test -s '${output_id}.${caller_name}.bcftools_stats.txt'
    test -s '${output_id}.${caller_name}.variant_summary.tsv'
    """
}
