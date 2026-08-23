process MANTA_REGIONS {
    label 'process_low'
    tag 'Manta call regions'
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    val calling_interval

    output:
    tuple val(calling_interval), path('manta_regions'), emit: regions

    script:
    """
    set -euo pipefail
    mkdir manta_regions
    if [[ -n '${calling_interval}' ]]; then
        region='${calling_interval}'
        contig="\${region%%:*}"
        range="\${region#*:}"
        start="\${range%-*}"
        end="\${range#*-}"
        printf '%s\t%s\t%s\n' "\$contig" "\$((start - 1))" "\$end" | bgzip -c > manta_regions/calling_interval.bed.gz
        tabix -p bed manta_regions/calling_interval.bed.gz
    else
        touch manta_regions/.unrestricted
    fi
    """
}

process MANTA_GERMLINE {
    label 'process_high'
    tag "${meta.sample_id}:manta"
    container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_3'
    publishDir "${params.outdir ?: 'results'}/structural_variants/germline", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam), path(bai), path(reference_fa), path(reference_fai), path(reference_dict), val(interval), path(region_dir)

    output:
    tuple val(meta), val('diploidSV'), path("${meta.sample_id}.manta.diploidSV.vcf.gz"), path("${meta.sample_id}.manta.diploidSV.vcf.gz.tbi"), path("${meta.sample_id}.manta.candidateSV.vcf.gz"), path("${meta.sample_id}.manta.candidateSV.vcf.gz.tbi"), emit: calls

    script:
    def regions = interval ? "--callRegions='${region_dir}/calling_interval.bed.gz'" : ''
    // Manta 1.6.0 hard-codes a 4 GiB reservation for graph merge; compact jobs use far less in practice.
    def mantaMemoryGb = Math.max(4, task.memory.toGiga().toInteger())
    """
    set -euo pipefail
    configManta.py --bam='${bam}' --referenceFasta='${reference_fa}' ${regions} --runDir=manta_run
    python2 manta_run/runWorkflow.py -m local -j ${task.cpus} -g ${mantaMemoryGb}
    cp manta_run/results/variants/diploidSV.vcf.gz '${meta.sample_id}.manta.diploidSV.vcf.gz'
    cp manta_run/results/variants/diploidSV.vcf.gz.tbi '${meta.sample_id}.manta.diploidSV.vcf.gz.tbi'
    cp manta_run/results/variants/candidateSV.vcf.gz '${meta.sample_id}.manta.candidateSV.vcf.gz'
    cp manta_run/results/variants/candidateSV.vcf.gz.tbi '${meta.sample_id}.manta.candidateSV.vcf.gz.tbi'
    test -s '${meta.sample_id}.manta.diploidSV.vcf.gz'
    test -s '${meta.sample_id}.manta.diploidSV.vcf.gz.tbi'
    """
}

process MANTA_SOMATIC {
    label 'process_high'
    tag "${meta.patient_id}:manta"
    container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_3'
    publishDir "${params.outdir ?: 'results'}/structural_variants/somatic", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(reference_fa), path(reference_fai), path(reference_dict), val(interval), path(region_dir)

    output:
    tuple val(meta), val('somaticSV'), path("${meta.patient_id}.manta.somaticSV.vcf.gz"), path("${meta.patient_id}.manta.somaticSV.vcf.gz.tbi"), path("${meta.patient_id}.manta.normal.diploidSV.vcf.gz"), path("${meta.patient_id}.manta.normal.diploidSV.vcf.gz.tbi"), path("${meta.patient_id}.manta.candidateSV.vcf.gz"), path("${meta.patient_id}.manta.candidateSV.vcf.gz.tbi"), emit: calls

    script:
    def regions = interval ? "--callRegions='${region_dir}/calling_interval.bed.gz'" : ''
    // Manta 1.6.0 hard-codes a 4 GiB reservation for graph merge; compact jobs use far less in practice.
    def mantaMemoryGb = Math.max(4, task.memory.toGiga().toInteger())
    """
    set -euo pipefail
    configManta.py --normalBam='${normal_bam}' --tumorBam='${tumor_bam}' --referenceFasta='${reference_fa}' ${regions} --runDir=manta_run
    python2 manta_run/runWorkflow.py -m local -j ${task.cpus} -g ${mantaMemoryGb}
    cp manta_run/results/variants/somaticSV.vcf.gz '${meta.patient_id}.manta.somaticSV.vcf.gz'
    cp manta_run/results/variants/somaticSV.vcf.gz.tbi '${meta.patient_id}.manta.somaticSV.vcf.gz.tbi'
    cp manta_run/results/variants/diploidSV.vcf.gz '${meta.patient_id}.manta.normal.diploidSV.vcf.gz'
    cp manta_run/results/variants/diploidSV.vcf.gz.tbi '${meta.patient_id}.manta.normal.diploidSV.vcf.gz.tbi'
    cp manta_run/results/variants/candidateSV.vcf.gz '${meta.patient_id}.manta.candidateSV.vcf.gz'
    cp manta_run/results/variants/candidateSV.vcf.gz.tbi '${meta.patient_id}.manta.candidateSV.vcf.gz.tbi'
    test -s '${meta.patient_id}.manta.somaticSV.vcf.gz'
    test -s '${meta.patient_id}.manta.somaticSV.vcf.gz.tbi'
    """
}

process SV_SUMMARY {
    label 'process_low'
    tag "${meta.sample_id ?: meta.patient_id}:Manta summary"
    container 'python:3.11@sha256:d0199e2a90bf7a206a485b115323a75bc946f30b463d704c5435a454aca084dd'
    publishDir "${params.outdir ?: 'results'}/structural_variants/summaries", mode: 'copy', overwrite: true

    input:
    tuple val(meta), val(output_semantics), path(final_vcf), path(final_tbi), path(reference_fai), val(analysis_type), val(assembly)

    output:
    tuple val(meta), path('*.sv_summary.tsv'), path('*.sv_class_summary.tsv'), path('*.sv_size_summary.tsv'), path('*.sv_chromosome_summary.tsv'), emit: summaries

    script:
    def identifier = analysis_type == 'germline' ? meta.sample_id : meta.patient_id
    def sample = meta.sample_id ?: meta.tumor_sample_id
    def patient = meta.patient_id
    """
    set -euo pipefail
    summarize_sv.py --vcf '${final_vcf}' --fai '${reference_fai}' --sample '${sample}' --patient '${patient}' --analysis '${analysis_type}' --semantics '${output_semantics}' --assembly '${assembly}' --prefix '${identifier}.manta'
    """
}
