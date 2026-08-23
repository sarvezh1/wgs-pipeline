process MUTECT2 {
    label 'process_high'
    tag "${meta.patient_id}:${meta.tumor_sample_id}/${meta.normal_sample_id}"

    publishDir "${params.outdir ?: 'results'}/somatic", mode: 'copy', overwrite: true
    container 'broadinstitute/gatk:4.6.2.0'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(reference_fa), path(reference_fai), path(reference_dict), val(expected_tumor_sample), val(expected_normal_sample), val(interval)

    output:
    tuple val(meta), val('mutect2'), path("${meta.patient_id}.mutect2.unfiltered.vcf.gz"), path("${meta.patient_id}.mutect2.unfiltered.vcf.gz.tbi"), path("${meta.patient_id}.mutect2.unfiltered.vcf.gz.stats"), emit: unfiltered

    script:
    def interval_arg = interval ? "--intervals '${interval}'" : ''
    """
    set -euo pipefail

    gatk --java-options '-Xmx768m' GetSampleName \
        --input '${tumor_bam}' \
        --output tumor.sample_name.txt
    gatk --java-options '-Xmx768m' GetSampleName \
        --input '${normal_bam}' \
        --output normal.sample_name.txt

    tumor_sample=\$(tr -d '\\r\\n' < tumor.sample_name.txt)
    normal_sample=\$(tr -d '\\r\\n' < normal.sample_name.txt)
    [[ -n "\$tumor_sample" ]] || { echo 'MUTECT2 ERROR: tumor BAM has no SM sample name' >&2; exit 1; }
    [[ -n "\$normal_sample" ]] || { echo 'MUTECT2 ERROR: normal BAM has no SM sample name' >&2; exit 1; }
    [[ "\$tumor_sample" == '${expected_tumor_sample}' ]] || { echo "MUTECT2 ERROR: tumor SM '\$tumor_sample' does not match expected '${expected_tumor_sample}'" >&2; exit 1; }
    [[ "\$normal_sample" == '${expected_normal_sample}' ]] || { echo "MUTECT2 ERROR: normal SM '\$normal_sample' does not match expected '${expected_normal_sample}'" >&2; exit 1; }
    [[ "\$tumor_sample" != "\$normal_sample" ]] || { echo 'MUTECT2 ERROR: tumor and normal SM values are identical' >&2; exit 1; }

    gatk --java-options '-Xmx768m' Mutect2 \
        --reference '${reference_fa}' \
        --input '${tumor_bam}' \
        --input '${normal_bam}' \
        --tumor-sample "\$tumor_sample" \
        --normal-sample "\$normal_sample" \
        ${interval_arg} \
        --output '${meta.patient_id}.mutect2.unfiltered.vcf.gz' \
        --create-output-variant-index false

    gatk --java-options '-Xmx768m' IndexFeatureFile \
        --input '${meta.patient_id}.mutect2.unfiltered.vcf.gz' \
        --output '${meta.patient_id}.mutect2.unfiltered.vcf.gz.tbi'

    test -s '${meta.patient_id}.mutect2.unfiltered.vcf.gz'
    test -s '${meta.patient_id}.mutect2.unfiltered.vcf.gz.tbi'
    test -s '${meta.patient_id}.mutect2.unfiltered.vcf.gz.stats'
    """
}

process MUTECT2_PRODUCTION {
    label 'process_high'
    tag "${meta.patient_id}:${meta.tumor_sample_id}/${meta.normal_sample_id}:production"

    publishDir "${params.outdir ?: 'results'}/somatic", mode: 'copy', overwrite: true, pattern: '*.vcf.gz*'
    publishDir "${params.outdir ?: 'results'}/somatic/orientation_bias", mode: 'copy', overwrite: true, pattern: '*.f1r2.tar.gz'
    container 'broadinstitute/gatk:4.6.2.0'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(reference_fa), path(reference_fai), path(reference_dict), val(expected_tumor_sample), val(expected_normal_sample), val(interval), path(germline_resource), path(germline_resource_tbi), path(panel_of_normals), path(panel_of_normals_tbi)

    output:
    tuple val(meta), val('mutect2'), path("${meta.patient_id}.mutect2.unfiltered.vcf.gz"), path("${meta.patient_id}.mutect2.unfiltered.vcf.gz.tbi"), path("${meta.patient_id}.mutect2.unfiltered.vcf.gz.stats"), emit: unfiltered
    tuple val(meta), path("${meta.patient_id}.mutect2.f1r2.tar.gz"), emit: f1r2

    script:
    def interval_arg = interval ? "--intervals '${interval}'" : ''
    """
    set -euo pipefail

    gatk --java-options '-Xmx768m' GetSampleName \
        --input '${tumor_bam}' \
        --output tumor.sample_name.txt
    gatk --java-options '-Xmx768m' GetSampleName \
        --input '${normal_bam}' \
        --output normal.sample_name.txt

    tumor_sample=\$(tr -d '\\r\\n' < tumor.sample_name.txt)
    normal_sample=\$(tr -d '\\r\\n' < normal.sample_name.txt)
    [[ -n "\$tumor_sample" ]] || { echo 'MUTECT2 ERROR: tumor BAM has no SM sample name' >&2; exit 1; }
    [[ -n "\$normal_sample" ]] || { echo 'MUTECT2 ERROR: normal BAM has no SM sample name' >&2; exit 1; }
    [[ "\$tumor_sample" == '${expected_tumor_sample}' ]] || { echo "MUTECT2 ERROR: tumor SM '\$tumor_sample' does not match expected '${expected_tumor_sample}'" >&2; exit 1; }
    [[ "\$normal_sample" == '${expected_normal_sample}' ]] || { echo "MUTECT2 ERROR: normal SM '\$normal_sample' does not match expected '${expected_normal_sample}'" >&2; exit 1; }
    [[ "\$tumor_sample" != "\$normal_sample" ]] || { echo 'MUTECT2 ERROR: tumor and normal SM values are identical' >&2; exit 1; }

    gatk --java-options '-Xmx768m' Mutect2 \
        --reference '${reference_fa}' \
        --input '${tumor_bam}' \
        --input '${normal_bam}' \
        --tumor-sample "\$tumor_sample" \
        --normal-sample "\$normal_sample" \
        --germline-resource '${germline_resource}' \
        --panel-of-normals '${panel_of_normals}' \
        --f1r2-tar-gz '${meta.patient_id}.mutect2.f1r2.tar.gz' \
        ${interval_arg} \
        --output '${meta.patient_id}.mutect2.unfiltered.vcf.gz' \
        --create-output-variant-index false

    gatk --java-options '-Xmx768m' IndexFeatureFile \
        --input '${meta.patient_id}.mutect2.unfiltered.vcf.gz' \
        --output '${meta.patient_id}.mutect2.unfiltered.vcf.gz.tbi'

    test -s '${meta.patient_id}.mutect2.unfiltered.vcf.gz'
    test -s '${meta.patient_id}.mutect2.unfiltered.vcf.gz.tbi'
    test -s '${meta.patient_id}.mutect2.unfiltered.vcf.gz.stats'
    test -s '${meta.patient_id}.mutect2.f1r2.tar.gz'
    """
}

process LEARN_READ_ORIENTATION_MODEL {
    label 'process_low'
    tag "${meta.patient_id}:orientation-model"

    publishDir "${params.outdir ?: 'results'}/somatic/orientation_bias", mode: 'copy', overwrite: true
    container 'broadinstitute/gatk:4.6.2.0'

    input:
    tuple val(meta), path(f1r2_tar)

    output:
    tuple val(meta), path("${meta.patient_id}.mutect2.read-orientation-model.tar.gz"), emit: priors

    script:
    """
    set -euo pipefail
    test -s '${f1r2_tar}' || { echo 'ORIENTATION MODEL ERROR: required Mutect2 F1R2 archive is missing or empty' >&2; exit 1; }
    gatk --java-options '-Xmx768m' LearnReadOrientationModel \
        --input '${f1r2_tar}' \
        --output '${meta.patient_id}.mutect2.read-orientation-model.tar.gz'
    test -s '${meta.patient_id}.mutect2.read-orientation-model.tar.gz'
    """
}

process GET_PILEUP_SUMMARIES {
    label 'process_medium'
    tag "${meta.patient_id}:tumor-normal-pileups"

    publishDir "${params.outdir ?: 'results'}/somatic/contamination", mode: 'copy', overwrite: true
    container 'broadinstitute/gatk:4.6.2.0'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(reference_fa), path(reference_fai), path(reference_dict), val(interval), path(contamination_sites), path(contamination_sites_tbi)

    output:
    tuple val(meta), path("${meta.patient_id}.tumor.pileups.table"), path("${meta.patient_id}.normal.pileups.table"), emit: pileups

    script:
    def interval_args = interval \
        ? "--intervals '${contamination_sites}' --intervals '${interval}' --interval-set-rule INTERSECTION" \
        : "--intervals '${contamination_sites}'"
    """
    set -euo pipefail

    gatk --java-options '-Xmx768m' GetPileupSummaries \
        --reference '${reference_fa}' \
        --input '${tumor_bam}' \
        --variant '${contamination_sites}' \
        ${interval_args} \
        --output '${meta.patient_id}.tumor.pileups.table'

    gatk --java-options '-Xmx768m' GetPileupSummaries \
        --reference '${reference_fa}' \
        --input '${normal_bam}' \
        --variant '${contamination_sites}' \
        ${interval_args} \
        --output '${meta.patient_id}.normal.pileups.table'

    test -s '${meta.patient_id}.tumor.pileups.table'
    test -s '${meta.patient_id}.normal.pileups.table'
    """
}

process CALCULATE_CONTAMINATION {
    label 'process_low'
    tag "${meta.patient_id}:contamination"

    publishDir "${params.outdir ?: 'results'}/somatic/contamination", mode: 'copy', overwrite: true
    container 'broadinstitute/gatk:4.6.2.0'

    input:
    tuple val(meta), path(tumor_pileups), path(normal_pileups)

    output:
    tuple val(meta), path("${meta.patient_id}.contamination.table"), path("${meta.patient_id}.tumor-segmentation.table"), emit: estimates

    script:
    """
    set -euo pipefail
    gatk --java-options '-Xmx768m' CalculateContamination \
        --input '${tumor_pileups}' \
        --matched-normal '${normal_pileups}' \
        --output '${meta.patient_id}.contamination.table' \
        --tumor-segmentation '${meta.patient_id}.tumor-segmentation.table'
    test -s '${meta.patient_id}.contamination.table'
    test -s '${meta.patient_id}.tumor-segmentation.table'
    """
}

process FILTERMUTECTCALLS {
    label 'process_medium'
    tag "${meta.patient_id}"

    publishDir "${params.outdir ?: 'results'}/somatic", mode: 'copy', overwrite: true
    container 'broadinstitute/gatk:4.6.2.0'

    input:
    tuple val(meta), val(caller_name), path(unfiltered_vcf), path(unfiltered_tbi), path(unfiltered_stats), path(reference_fa), path(reference_fai), path(reference_dict)

    output:
    tuple val(meta), val(caller_name), path("${meta.patient_id}.mutect2.filtered.vcf.gz"), path("${meta.patient_id}.mutect2.filtered.vcf.gz.tbi"), emit: calls

    script:
    """
    set -euo pipefail

    gatk --java-options '-Xmx768m' FilterMutectCalls \
        --reference '${reference_fa}' \
        --variant '${unfiltered_vcf}' \
        --output '${meta.patient_id}.mutect2.filtered.vcf.gz' \
        --create-output-variant-index false

    gatk --java-options '-Xmx768m' IndexFeatureFile \
        --input '${meta.patient_id}.mutect2.filtered.vcf.gz' \
        --output '${meta.patient_id}.mutect2.filtered.vcf.gz.tbi'

    test -s '${meta.patient_id}.mutect2.filtered.vcf.gz'
    test -s '${meta.patient_id}.mutect2.filtered.vcf.gz.tbi'
    """
}

process FILTERMUTECTCALLS_PRODUCTION {
    label 'process_medium'
    tag "${meta.patient_id}:production"

    publishDir "${params.outdir ?: 'results'}/somatic", mode: 'copy', overwrite: true
    container 'broadinstitute/gatk:4.6.2.0'

    input:
    tuple val(meta), val(caller_name), path(unfiltered_vcf), path(unfiltered_tbi), path(unfiltered_stats), path(reference_fa), path(reference_fai), path(reference_dict), path(orientation_priors), path(contamination_table), path(tumor_segmentation)

    output:
    tuple val(meta), val(caller_name), path("${meta.patient_id}.mutect2.filtered.vcf.gz"), path("${meta.patient_id}.mutect2.filtered.vcf.gz.tbi"), emit: calls

    script:
    """
    set -euo pipefail
    test -s '${orientation_priors}' || { echo 'FILTERMUTECTCALLS ERROR: orientation-bias priors are required for the production profile' >&2; exit 1; }
    test -s '${contamination_table}' || { echo 'FILTERMUTECTCALLS ERROR: contamination table is required for the production profile' >&2; exit 1; }
    test -s '${tumor_segmentation}' || { echo 'FILTERMUTECTCALLS ERROR: tumor segmentation is required for the production profile' >&2; exit 1; }

    gatk --java-options '-Xmx768m' FilterMutectCalls \
        --reference '${reference_fa}' \
        --variant '${unfiltered_vcf}' \
        --stats '${unfiltered_stats}' \
        --orientation-bias-artifact-priors '${orientation_priors}' \
        --contamination-table '${contamination_table}' \
        --tumor-segmentation '${tumor_segmentation}' \
        --output '${meta.patient_id}.mutect2.filtered.vcf.gz' \
        --create-output-variant-index false

    gatk --java-options '-Xmx768m' IndexFeatureFile \
        --input '${meta.patient_id}.mutect2.filtered.vcf.gz' \
        --output '${meta.patient_id}.mutect2.filtered.vcf.gz.tbi'

    test -s '${meta.patient_id}.mutect2.filtered.vcf.gz'
    test -s '${meta.patient_id}.mutect2.filtered.vcf.gz.tbi'
    """
}
