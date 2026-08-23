process CNV_COLLECT_READ_COUNTS {
    label 'process_high'
    tag "${meta.sample_id}:CNV read counts"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir ?: 'results'}/copy_number/counts", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam), path(bai), path(intervals)

    output:
    tuple val(meta), path("${meta.sample_id}.counts.hdf5"), emit: counts

    script:
    """
    set -euo pipefail
    gatk CollectReadCounts -I '${bam}' -L '${intervals}' --interval-merging-rule OVERLAPPING_ONLY -O '${meta.sample_id}.counts.hdf5'
    test -s '${meta.sample_id}.counts.hdf5'
    """
}

process GERMLINE_CONTIG_PLOIDY_CASE {
    label 'process_high'
    tag "${meta.sample_id}:gCNV contig ploidy"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir ?: 'results'}/copy_number/germline/ploidy", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(counts), path(ploidy_model)

    output:
    tuple val(meta), path(counts), path("${meta.sample_id}.ploidy-calls"), emit: calls

    script:
    """
    set -euo pipefail
    gatk DetermineGermlineContigPloidy --model '${ploidy_model}' --input '${counts}' --output . --output-prefix '${meta.sample_id}.ploidy'
    test -d '${meta.sample_id}.ploidy-calls'
    """
}

process GERMLINE_CNV_CALLER_CASE {
    label 'process_high'
    tag "${meta.sample_id}:gCNV CASE"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir ?: 'results'}/copy_number/germline/calls", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(counts), path(ploidy_calls), path(gcnv_model)

    output:
    tuple val(meta), path(counts), path(ploidy_calls), path("${meta.sample_id}.gcnv-calls"), emit: calls

    script:
    """
    set -euo pipefail
    gatk GermlineCNVCaller --run-mode CASE --contig-ploidy-calls '${ploidy_calls}' --model '${gcnv_model}' --input '${counts}' --output . --output-prefix '${meta.sample_id}.gcnv'
    test -d '${meta.sample_id}.gcnv-calls'
    """
}

process POSTPROCESS_GERMLINE_CNV {
    label 'process_high'
    tag "${meta.sample_id}:gCNV postprocess"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir ?: 'results'}/copy_number/germline/final", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(counts), path(ploidy_calls), path(gcnv_calls), path(gcnv_model), val(allosomal_contigs)

    output:
    tuple val(meta), path("${meta.sample_id}.gcnv.intervals.vcf"), path("${meta.sample_id}.gcnv.segments.vcf"), path("${meta.sample_id}.gcnv.denoised_copy_ratios.tsv"), emit: results

    script:
    def allosomalArgs = allosomal_contigs.collect { contig -> "--allosomal-contig '${contig}'" }.join(' ')
    """
    set -euo pipefail
    gatk PostprocessGermlineCNVCalls \
        --model-shard-path '${gcnv_model}' \
        --calls-shard-path '${gcnv_calls}' \
        --contig-ploidy-calls '${ploidy_calls}' \
        --sample-index 0 \
        --autosomal-ref-copy-number 2 \
        ${allosomalArgs} \
        --output-genotyped-intervals '${meta.sample_id}.gcnv.intervals.vcf' \
        --output-genotyped-segments '${meta.sample_id}.gcnv.segments.vcf' \
        --output-denoised-copy-ratios '${meta.sample_id}.gcnv.denoised_copy_ratios.tsv'
    test -s '${meta.sample_id}.gcnv.intervals.vcf'
    test -s '${meta.sample_id}.gcnv.segments.vcf'
    """
}

process SOMATIC_CNV_DENOISE {
    label 'process_high'
    tag "${meta.sample_id}:somatic CNV denoise"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir ?: 'results'}/copy_number/somatic/copy_ratios", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(counts), path(read_count_pon)

    output:
    tuple val(meta), path("${meta.sample_id}.standardizedCR.tsv"), path("${meta.sample_id}.denoisedCR.tsv"), emit: ratios

    script:
    """
    set -euo pipefail
    gatk DenoiseReadCounts -I '${counts}' --count-panel-of-normals '${read_count_pon}' --standardized-copy-ratios '${meta.sample_id}.standardizedCR.tsv' --denoised-copy-ratios '${meta.sample_id}.denoisedCR.tsv'
    test -s '${meta.sample_id}.denoisedCR.tsv'
    """
}

process CNV_COLLECT_ALLELIC_COUNTS {
    label 'process_high'
    tag "${meta.sample_id}:CNV allelic counts"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir ?: 'results'}/copy_number/counts", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam), path(bai), path(reference_fa), path(reference_fai), path(reference_dict), path(common_sites), path(common_sites_tbi)

    output:
    tuple val(meta), path("${meta.sample_id}.allelic_counts.tsv"), emit: counts

    script:
    """
    set -euo pipefail
    gatk CollectAllelicCounts -I '${bam}' -R '${reference_fa}' -L '${common_sites}' -O '${meta.sample_id}.allelic_counts.tsv'
    test -s '${meta.sample_id}.allelic_counts.tsv'
    """
}

process SOMATIC_CNV_MODEL_SEGMENTS {
    label 'process_high'
    tag "${meta.patient_id}:somatic CNV segments"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir ?: 'results'}/copy_number/somatic/model", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(denoised_copy_ratios), path(tumor_allelic_counts), path(normal_allelic_counts)

    output:
    tuple val(meta), path("${meta.patient_id}.cr.seg"), path("${meta.patient_id}.modelFinal.seg"), path("${meta.patient_id}.hets.tsv"), emit: segments

    script:
    """
    set -euo pipefail
    gatk ModelSegments --denoised-copy-ratios '${denoised_copy_ratios}' --allelic-counts '${tumor_allelic_counts}' --normal-allelic-counts '${normal_allelic_counts}' --output . --output-prefix '${meta.patient_id}'
    test -s '${meta.patient_id}.cr.seg'
    test -s '${meta.patient_id}.modelFinal.seg'
    """
}

process SOMATIC_CNV_CALL_SEGMENTS {
    label 'process_medium'
    tag "${meta.patient_id}:somatic CNV calls"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir ?: 'results'}/copy_number/somatic/final", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(copy_ratio_segments), path(model_segments), path(hets)

    output:
    tuple val(meta), path(copy_ratio_segments), path(model_segments), path(hets), path("${meta.patient_id}.called.seg"), emit: results

    script:
    """
    set -euo pipefail
    gatk CallCopyRatioSegments -I '${copy_ratio_segments}' -O '${meta.patient_id}.called.seg'
    test -s '${meta.patient_id}.called.seg'
    """
}

process CNV_SUMMARY {
    label 'process_low'
    tag "${meta.sample_id ?: meta.patient_id}:CNV summary"
    container 'python:3.11@sha256:d0199e2a90bf7a206a485b115323a75bc946f30b463d704c5435a454aca084dd'
    publishDir "${params.outdir ?: 'results'}/copy_number/summaries", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(final_segments), val(input_format), val(analysis_type), val(assembly)

    output:
    tuple val(meta), path('*.cnv_summary.tsv'), emit: summary

    script:
    def identifier = analysis_type == 'germline' ? meta.sample_id : meta.patient_id
    """
    set -euo pipefail
    summarize_cnv.py --input '${final_segments}' --format '${input_format}' --id '${identifier}' --analysis '${analysis_type}' --assembly '${assembly}' --output '${identifier}.cnv_summary.tsv'
    """
}
