process GATK_HAPLOTYPECALLER {
    label 'process_high'
    tag "${meta.sample_id}:${meta.role}"

    publishDir "${params.outdir ?: 'results'}/variants", mode: 'copy', overwrite: true
    container 'broadinstitute/gatk:4.6.2.0'

    input:
    tuple val(meta), val(caller_name), path(caller_bam), path(caller_bai), path(reference_fa), path(reference_fai), path(reference_dict), val(interval)

    output:
    tuple val(meta), val(caller_name), path("${meta.sample_id}.haplotypecaller.vcf.gz"), path("${meta.sample_id}.haplotypecaller.vcf.gz.tbi"), emit: calls

    script:
    def interval_arg = interval ? "--intervals '${interval}'" : ''
    """
    set -euo pipefail

    gatk --java-options '-Xmx768m' HaplotypeCaller \
        --reference '${reference_fa}' \
        --input '${caller_bam}' \
        --output '${meta.sample_id}.haplotypecaller.vcf.gz' \
        ${interval_arg} \
        --native-pair-hmm-threads ${task.cpus} \
        --create-output-variant-index false

    gatk --java-options '-Xmx768m' IndexFeatureFile \
        --input '${meta.sample_id}.haplotypecaller.vcf.gz' \
        --output '${meta.sample_id}.haplotypecaller.vcf.gz.tbi'
    """
}

process GATK_HAPLOTYPECALLER_GVCF {
    label 'process_high'
    tag "${meta.sample_id}:${meta.role}:gvcf"

    publishDir "${params.outdir ?: 'results'}/variants/gvcf", mode: 'copy', overwrite: true
    container 'broadinstitute/gatk:4.6.2.0'

    input:
    tuple val(meta), val(caller_name), path(caller_bam), path(caller_bai), path(reference_fa), path(reference_fai), path(reference_dict), val(interval)

    output:
    tuple val(meta), val(caller_name), path("${meta.sample_id}.haplotypecaller.g.vcf.gz"), path("${meta.sample_id}.haplotypecaller.g.vcf.gz.tbi"), emit: gvcfs

    script:
    def interval_arg = interval ? "--intervals '${interval}'" : ''
    """
    set -euo pipefail

    gatk --java-options '-Xmx768m' HaplotypeCaller \
        --reference '${reference_fa}' \
        --input '${caller_bam}' \
        --output '${meta.sample_id}.haplotypecaller.g.vcf.gz' \
        --emit-ref-confidence GVCF \
        ${interval_arg} \
        --native-pair-hmm-threads ${task.cpus} \
        --create-output-variant-index false

    gatk --java-options '-Xmx768m' IndexFeatureFile \
        --input '${meta.sample_id}.haplotypecaller.g.vcf.gz' \
        --output '${meta.sample_id}.haplotypecaller.g.vcf.gz.tbi'

    awk '/^##GVCFBlock/ { block=1 } /^#CHROM/ { header=1 } /<NON_REF>/ { nonref=1 } block && header && nonref { found=1; exit } END { exit(found ? 0 : 1) }' \
        < <(gzip -dc '${meta.sample_id}.haplotypecaller.g.vcf.gz')
    test -s '${meta.sample_id}.haplotypecaller.g.vcf.gz.tbi'
    """
}

process GATK_GENOTYPE_GVCF {
    label 'process_medium'
    tag "${meta.sample_id}:genotype-gvcf"

    publishDir "${params.outdir ?: 'results'}/variants", mode: 'copy', overwrite: true
    container 'broadinstitute/gatk:4.6.2.0'

    input:
    tuple val(meta), val(caller_name), path(gvcf), path(gvcf_tbi), path(reference_fa), path(reference_fai), path(reference_dict), val(interval)

    output:
    tuple val(meta), val(caller_name), path("${meta.sample_id}.haplotypecaller.vcf.gz"), path("${meta.sample_id}.haplotypecaller.vcf.gz.tbi"), emit: calls

    script:
    def interval_arg = interval ? "--intervals '${interval}'" : ''
    """
    set -euo pipefail

    gatk --java-options '-Xmx768m' GenotypeGVCFs \
        --reference '${reference_fa}' \
        --variant '${gvcf}' \
        --output '${meta.sample_id}.haplotypecaller.vcf.gz' \
        ${interval_arg} \
        --create-output-variant-index false

    gatk --java-options '-Xmx768m' IndexFeatureFile \
        --input '${meta.sample_id}.haplotypecaller.vcf.gz' \
        --output '${meta.sample_id}.haplotypecaller.vcf.gz.tbi'

    test -s '${meta.sample_id}.haplotypecaller.vcf.gz'
    test -s '${meta.sample_id}.haplotypecaller.vcf.gz.tbi'
    """
}
