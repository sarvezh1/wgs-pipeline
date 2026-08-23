process ANALYSIS_REPORT {
    label 'process_low'
    tag "${report_id}:${analysis_type}:integrated report"

    publishDir "${params.outdir ?: 'results'}/report", mode: 'copy', overwrite: true
    container 'python:3.11@sha256:d0199e2a90bf7a206a485b115323a75bc946f30b463d704c5435a454aca084dd'

    input:
    tuple path(report_inputs, stageAs: 'report_inputs/*'), path(samplesheet), val(report_id), val(analysis_type), val(caller), val(assembly), val(resource_mode), val(pipeline_version), val(nextflow_version), val(run_bqsr), val(run_benchmark), val(run_vep), val(run_sv), val(run_cnv), val(run_coverage_qc), val(run_multiqc), val(run_vcf2circos)

    output:
    tuple val(report_id), path("${report_id}.${analysis_type}.wgs_report.html"), path('figures/*.svg', optional: true), emit: report

    script:
    def optional_flags = [
        run_bqsr ? '--run-bqsr' : '',
        run_benchmark ? '--run-benchmark' : '',
        run_vep ? '--run-vep' : '',
        run_sv ? '--run-sv' : '',
        run_cnv ? '--run-cnv' : '',
        run_coverage_qc ? '--run-coverage-qc' : '',
    ].findAll { flag -> flag }.join(' ')
    def multiqc_link = run_multiqc ? "--multiqc-link '../multiqc/${analysis_type}.multiqc.html'" : ''
    def vcf2circos_link = run_vcf2circos ? "--vcf2circos-link '../vcf2circos/${report_id}.${caller == 'compare' ? (analysis_type == 'germline' ? 'haplotypecaller' : 'mutect2') : caller}.vcf2circos.html'" : ''
    """
    set -euo pipefail
    generate_analysis_report.py \
        --input-root report_inputs \
        --output '${report_id}.${analysis_type}.wgs_report.html' \
        --samplesheet '${samplesheet}' \
        --report-id '${report_id}' \
        --analysis-type '${analysis_type}' \
        --caller '${caller}' \
        --assembly '${assembly}' \
        --resource-mode '${resource_mode}' \
        --pipeline-version '${pipeline_version}' \
        --nextflow-version '${nextflow_version}' \
        ${optional_flags} ${multiqc_link} ${vcf2circos_link}
    test -s '${report_id}.${analysis_type}.wgs_report.html'
    """
}
