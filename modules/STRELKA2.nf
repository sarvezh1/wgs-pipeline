process STRELKA2_REGIONS {
    label 'process_low'
    tag "Strelka2 call regions"

    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    val calling_interval

    output:
    tuple val(calling_interval), path('strelka_regions'), emit: indexed_regions

    script:
    """
    set -euo pipefail
    mkdir -p strelka_regions
    if [[ -n '${calling_interval}' ]]; then
        region='${calling_interval}'
        region_contig="\${region%%:*}"
        region_range="\${region#*:}"
        region_start="\${region_range%-*}"
        region_end="\${region_range#*-}"
        printf '%s\t%s\t%s\n' \
            "\${region_contig}" \
            "\$((region_start - 1))" \
            "\${region_end}" \
            | bgzip -c > strelka_regions/calling_interval.bed.gz
        tabix -p bed strelka_regions/calling_interval.bed.gz
    else
        touch strelka_regions/.unrestricted
    fi
    """
}

process STRELKA2 {
    label 'process_high'
    tag "${meta.patient_id}:${meta.tumor_sample_id}/${meta.normal_sample_id}"

    publishDir "${params.outdir ?: 'results'}/somatic", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/strelka:2.9.10--h9ee0642_1'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(reference_fa), path(reference_fai), path(reference_dict), val(interval), path(call_regions_dir), val(expected_tumor_sample), val(expected_normal_sample)

    output:
    tuple val(meta), val('strelka2'), path("${meta.patient_id}.strelka2.somatic.snvs.vcf.gz"), path("${meta.patient_id}.strelka2.somatic.snvs.vcf.gz.tbi"), path("${meta.patient_id}.strelka2.somatic.indels.vcf.gz"), path("${meta.patient_id}.strelka2.somatic.indels.vcf.gz.tbi"), emit: raw_calls

    script:
    def call_regions_arg = interval ? "--callRegions='${call_regions_dir}/calling_interval.bed.gz'" : ''
    """
    set -euo pipefail

    configureStrelkaSomaticWorkflow.py \
        --normalBam='${normal_bam}' \
        --tumorBam='${tumor_bam}' \
        --referenceFasta='${reference_fa}' \
        ${call_regions_arg} \
        --runDir='strelka_run'

    python2 strelka_run/runWorkflow.py -m local -j ${task.cpus}

    cp strelka_run/results/variants/somatic.snvs.vcf.gz '${meta.patient_id}.strelka2.somatic.snvs.vcf.gz'
    cp strelka_run/results/variants/somatic.snvs.vcf.gz.tbi '${meta.patient_id}.strelka2.somatic.snvs.vcf.gz.tbi'
    cp strelka_run/results/variants/somatic.indels.vcf.gz '${meta.patient_id}.strelka2.somatic.indels.vcf.gz'
    cp strelka_run/results/variants/somatic.indels.vcf.gz.tbi '${meta.patient_id}.strelka2.somatic.indels.vcf.gz.tbi'

    test -s '${meta.patient_id}.strelka2.somatic.snvs.vcf.gz'
    test -s '${meta.patient_id}.strelka2.somatic.snvs.vcf.gz.tbi'
    test -s '${meta.patient_id}.strelka2.somatic.indels.vcf.gz'
    test -s '${meta.patient_id}.strelka2.somatic.indels.vcf.gz.tbi'
    """
}

process MERGE_STRELKA2_VCFS {
    label 'process_low'
    tag "${meta.patient_id}"

    publishDir "${params.outdir ?: 'results'}/somatic", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    tuple val(meta), val(caller_name), path(snv_vcf), path(snv_tbi), path(indel_vcf), path(indel_tbi)

    output:
    tuple val(meta), val(caller_name), path("${meta.patient_id}.strelka2.vcf.gz"), path("${meta.patient_id}.strelka2.vcf.gz.tbi"), emit: calls

    script:
    """
    set -euo pipefail

    bcftools concat \
        --allow-overlaps \
        --output-type z \
        --output '${meta.patient_id}.strelka2.vcf.gz' \
        '${snv_vcf}' '${indel_vcf}'
    bcftools index --tbi '${meta.patient_id}.strelka2.vcf.gz'
    bcftools view -H '${meta.patient_id}.strelka2.vcf.gz' >/dev/null

    test -s '${meta.patient_id}.strelka2.vcf.gz'
    test -s '${meta.patient_id}.strelka2.vcf.gz.tbi'
    """
}
