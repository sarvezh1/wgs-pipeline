process VEP {
    label 'process_high'
    tag "${meta.patient_id ?: meta.sample_id}:${caller_name}"

    publishDir "${params.outdir ?: 'results'}/vep", mode: 'copy', overwrite: true
    container 'ensemblorg/ensembl-vep:release_116.0'
    containerOptions '--user $(id -u):$(id -g)'

    input:
    tuple val(meta), val(caller_name), path(normalized_vcf), path(normalized_vcf_index), path(reference_fa), path(reference_fai), val(vep_mode), path(annotation_resources)

    output:
    tuple val(meta), val(caller_name), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.vep.vcf.gz"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.vep.vcf.gz.tbi"), emit: annotated

    script:
    def output_id = meta.patient_id ?: meta.sample_id
    def annotation_args

    if (vep_mode == 'cache') {
        def cache_dir = annotation_resources[0]
        annotation_args = """
            --cache \\
            --offline \\
            --dir_cache '${cache_dir}' \\
            --cache_version 116 \\
            --species homo_sapiens \\
            --assembly GRCh38 \\
            --check_existing
        """.stripIndent().trim()
    } else {
        def gtf = annotation_resources.find { resource -> !resource.name.endsWith('.tbi') }
        annotation_args = "--gtf '${gtf}'"
    }

    """
    set -euo pipefail

    test -s '${normalized_vcf}'
    test -s '${normalized_vcf_index}'

    vep \
        ${annotation_args} \
        --fasta '${reference_fa}' \
        --format vcf \
        --vcf \
        --symbol \
        --canonical \
        --biotype \
        --numbers \
        --variant_class \
        --force_overwrite \
        --no_stats \
        --fork 1 \
        --buffer_size 500 \
        --input_file '${normalized_vcf}' \
        --output_file '${output_id}.${caller_name}.vep.vcf.gz' \
        --compress_output bgzip

    tabix -f -p vcf '${output_id}.${caller_name}.vep.vcf.gz'

    gzip -cd '${output_id}.${caller_name}.vep.vcf.gz' | \
        awk '/^##INFO=<ID=CSQ/ { found=1 } END { exit(found ? 0 : 1) }'

    test -s '${output_id}.${caller_name}.vep.vcf.gz'
    test -s '${output_id}.${caller_name}.vep.vcf.gz.tbi'
    """
}
