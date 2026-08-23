process MULTIQC {
    label 'process_medium'
    tag "${analysis_type}"

    publishDir "${params.outdir ?: 'results'}/multiqc", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/multiqc:1.27.1--pyhdfd78af_0'

    input:
    tuple val(analysis_type), path(qc_files)

    output:
    tuple val(analysis_type), path("${analysis_type}.multiqc.html"), path("${analysis_type}.multiqc_data"), emit: report

    script:
    def report_title = "${analysis_type.capitalize()} WGS quality control"
    """
    set -euo pipefail

    mkdir qc_inputs
    for qc_file in ${qc_files.join(' ')}; do
        ln -s "\$(readlink -f "\${qc_file}")" "qc_inputs/\$(basename "\${qc_file}")"
    done

    cat > multiqc_config.yaml <<'YAML'
    title: "${report_title}"
    subtitle: "Shared pipeline QC aggregation"
    report_comment: "Reporting-only aggregation of native pipeline QC outputs."
    show_analysis_paths: false
    show_analysis_time: false
    export_plots: false
    YAML

    multiqc \
        --force \
        --data-dir \
        --filename '${analysis_type}.multiqc.html' \
        --outdir . \
        --config multiqc_config.yaml \
        qc_inputs

    test -s '${analysis_type}.multiqc.html'
    test -s '${analysis_type}.multiqc_data/multiqc_general_stats.txt'
    ! grep -F '/tmp/' '${analysis_type}.multiqc.html'
    ! grep -F '/work/' '${analysis_type}.multiqc.html'
    """
}
