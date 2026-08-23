process HAPPY_BENCHMARK {
    label 'process_high'
    tag "${meta.patient_id ?: meta.sample_id}:${caller_name}:hap.py"

    publishDir "${params.outdir ?: 'results'}/benchmark", mode: 'copy', overwrite: true, pattern: '*.benchmark_summary.tsv'
    publishDir "${params.outdir ?: 'results'}/benchmark", mode: 'copy', overwrite: true, pattern: '*.benchmark_provenance.tsv'
    publishDir "${params.outdir ?: 'results'}/benchmark/stratified", mode: 'copy', overwrite: true, pattern: '*.benchmark_stratified_*.tsv'
    publishDir "${params.outdir ?: 'results'}/benchmark/stratified", mode: 'copy', overwrite: true, pattern: '*.benchmark_stratification_provenance.tsv'
    publishDir "${params.outdir ?: 'results'}/benchmark/native_engine_outputs", mode: 'copy', overwrite: true, pattern: '*.happy.*'
    container 'mgibio/hap.py@sha256:492574403bf5cc648cf9c297957cc2eac60cbbc7da1a968984bafad060536bdc'

    input:
    tuple val(meta), val(caller_name), path(query_vcf), path(query_vcf_index), path(truth_vcf), path(truth_vcf_index), path(benchmark_regions), path(reference_fa), path(reference_fai), path(reference_dict), path(config_manifest), val(benchmark_name), val(benchmark_version), val(benchmark_mode), val(benchmark_assembly), path(stratification_resources), val(stratification_tsv_name)

    output:
    tuple val(meta), val(caller_name), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark_summary.tsv"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark_provenance.tsv"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.happy.*"), emit: results
    path "${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark_stratified_summary.tsv", optional: true, emit: stratified_summary
    path "${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark_stratified_errors.tsv", optional: true, emit: stratified_errors
    path "${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark_stratification_provenance.tsv", optional: true, emit: stratified_provenance

    script:
    def output_id = meta.patient_id ?: meta.sample_id
    def stratification_setup = stratification_tsv_name ? """
    awk -F '\\t' 'BEGIN { OFS="\\t" }
        NF == 0 || /^#/ { next }
        { n=split(\$2, parts, "/"); print \$1, parts[n] }
    ' '${stratification_tsv_name}' > staged.stratifications.tsv
    stratification_arg='--stratification staged.stratifications.tsv'
    """ : "stratification_arg=''"
    def stratified_parse = stratification_tsv_name ? """
    parse_happy_stratified.py \\
        --extended '${output_id}.${caller_name}.happy.extended.csv' \\
        --happy-vcf '${output_id}.${caller_name}.happy.vcf.gz' \\
        --manifest staged.stratifications.tsv \\
        --benchmark-regions '${benchmark_regions}' \\
        --output '${output_id}.${caller_name}.benchmark_stratified_summary.tsv' \\
        --errors '${output_id}.${caller_name}.benchmark_stratified_errors.tsv' \\
        --provenance '${output_id}.${caller_name}.benchmark_stratification_provenance.tsv' \\
        --sample '${output_id}' --caller '${caller_name}' \\
        --benchmark-name '${benchmark_name}' --benchmark-version '${benchmark_version}' \\
        --engine-version 'hap.py 0.3.15; RTG Tools 3.12.1' --assembly '${benchmark_assembly}'
    """ : ''
    """
    set -euo pipefail

    test -s '${query_vcf}'
    test -s '${query_vcf_index}'
    test -s '${truth_vcf}'
    test -s '${truth_vcf_index}'
    test -s '${benchmark_regions}'
    ${stratification_setup}

    export RTG_MEM=1g
    hap_py_gender=\$(awk -F '\t' 'NR==1 { for (i=1; i<=NF; i++) column[\$i]=i; next } NR==2 { print \$column["hap_py_gender_argument"]; exit }' '${config_manifest}')
    case "\$hap_py_gender" in
        male|female) gender_arg="--gender \$hap_py_gender" ;;
        not_passed) gender_arg='' ;;
        *) echo "BENCHMARK ERROR: invalid validated hap.py gender argument '\$hap_py_gender'" >&2; exit 1 ;;
    esac
    /opt/hap.py/bin/hap.py \
        '${truth_vcf}' \
        '${query_vcf}' \
        --false-positives '${benchmark_regions}' \
        --reference '${reference_fa}' \
        --report-prefix '${output_id}.${caller_name}.happy' \
        --engine vcfeval \
        --engine-vcfeval-path /opt/hap.py/libexec/rtg-tools-install/rtg \
        --threads ${task.cpus} \
        --write-vcf \
        --write-counts \
        \${gender_arg} \
        \${stratification_arg}

    test -s '${output_id}.${caller_name}.happy.summary.csv'
    test -s '${output_id}.${caller_name}.happy.extended.csv'
    benchmark_region_bases=\$(benchmark_region_bases.sh '${benchmark_regions}')
    parse_happy_summary.py \
        '${output_id}.${caller_name}.happy.summary.csv' \
        '${output_id}.${caller_name}.benchmark_summary.tsv' \
        '${output_id}' \
        'germline' \
        '${caller_name}' \
        '${benchmark_name}' \
        '${benchmark_version}' \
        '${benchmark_mode}' \
        "\${benchmark_region_bases}" \
        '${benchmark_assembly}'
    ${stratified_parse}

    printf 'sample_or_patient\tanalysis_type\tcaller\t' > '${output_id}.${caller_name}.benchmark_provenance.tsv'
    head -n 1 '${config_manifest}' >> '${output_id}.${caller_name}.benchmark_provenance.tsv'
    printf '%s\t%s\t%s\t' '${output_id}' 'germline' '${caller_name}' >> '${output_id}.${caller_name}.benchmark_provenance.tsv'
    tail -n 1 '${config_manifest}' >> '${output_id}.${caller_name}.benchmark_provenance.tsv'
    test -s '${output_id}.${caller_name}.benchmark_summary.tsv'
    test -s '${output_id}.${caller_name}.benchmark_provenance.tsv'
    """
}

process EXACT_ALLELE_BENCHMARK {
    label 'process_low'
    tag "${meta.patient_id ?: meta.sample_id}:${caller_name}:exact-${benchmark_mode}"

    publishDir "${params.outdir ?: 'results'}/benchmark", mode: 'copy', overwrite: true, pattern: '*.benchmark_summary.tsv'
    publishDir "${params.outdir ?: 'results'}/benchmark", mode: 'copy', overwrite: true, pattern: '*.benchmark_provenance.tsv'
    publishDir "${params.outdir ?: 'results'}/benchmark/stratified", mode: 'copy', overwrite: true, pattern: '*.benchmark_stratified_*.tsv'
    publishDir "${params.outdir ?: 'results'}/benchmark/stratified", mode: 'copy', overwrite: true, pattern: '*.benchmark_stratification_provenance.tsv'
    publishDir "${params.outdir ?: 'results'}/benchmark/native_engine_outputs", mode: 'copy', overwrite: true, pattern: '*.benchmark.*'
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    tuple val(meta), val(caller_name), path(query_vcf), path(query_vcf_index), path(truth_vcf), path(truth_vcf_index), path(benchmark_regions), path(reference_fa), path(reference_fai), path(reference_dict), path(config_manifest), val(analysis_type), val(benchmark_name), val(benchmark_version), val(benchmark_mode), val(benchmark_assembly)

    output:
    tuple val(meta), val(caller_name), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark_summary.tsv"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark_provenance.tsv"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark.*"), emit: results
    tuple val(meta), val(caller_name), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark.truth.alleles.tsv"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.benchmark.query.alleles.tsv"), emit: alleles

    script:
    def output_id = meta.patient_id ?: meta.sample_id
    def query_filter = analysis_type == 'somatic' ? '--apply-filters PASS' : ''
    """
    set -euo pipefail

    test -s '${query_vcf}'
    test -s '${query_vcf_index}'
    test -s '${truth_vcf}'
    test -s '${truth_vcf_index}'
    test -s '${benchmark_regions}'

    bcftools norm \
        --fasta-ref '${reference_fa}' \
        --multiallelics -any \
        --output-type z \
        --output truth.normalized.vcf.gz \
        '${truth_vcf}'
    bcftools index --tbi truth.normalized.vcf.gz
    bcftools view \
        --regions-file '${benchmark_regions}' \
        --output-type z \
        --output '${output_id}.${caller_name}.benchmark.truth.vcf.gz' \
        truth.normalized.vcf.gz
    bcftools index --tbi '${output_id}.${caller_name}.benchmark.truth.vcf.gz'

    bcftools norm \
        --fasta-ref '${reference_fa}' \
        --multiallelics -any \
        --output-type z \
        --output query.normalized.vcf.gz \
        '${query_vcf}'
    bcftools index --tbi query.normalized.vcf.gz
    bcftools view \
        ${query_filter} \
        --regions-file '${benchmark_regions}' \
        --output-type z \
        --output '${output_id}.${caller_name}.benchmark.query.vcf.gz' \
        query.normalized.vcf.gz
    bcftools index --tbi '${output_id}.${caller_name}.benchmark.query.vcf.gz'

    for source in truth query; do
        vcf="${output_id}.${caller_name}.benchmark.\$source.vcf.gz"
        alleles="${output_id}.${caller_name}.benchmark.\$source.alleles.tsv"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "\${vcf}" \
            | LC_ALL=C sort -u > "\${alleles}"
        awk -F '\\t' 'length(\$3)==1 && length(\$4)==1' \
            "\${alleles}" \
            > "\${source}.snv.tsv"
        awk -F '\\t' 'length(\$3)!=1 || length(\$4)!=1' \
            "\${alleles}" \
            > "\${source}.indel.tsv"
    done

    benchmark_region_bases=\$(benchmark_region_bases.sh '${benchmark_regions}')
    printf 'sample_or_patient\tanalysis_type\tcaller\tbenchmark_name\tbenchmark_version\tbenchmark_mode\tvariant_class\ttruth_tp\tquery_tp\tfp\tfn\tprecision\trecall\tf1\tbenchmark_region_bases\tengine\tengine_version\treference_assembly\n' \
        > '${output_id}.${caller_name}.benchmark_summary.tsv'

    write_metrics() {
        local variant_class=\$1 truth_keys=\$2 query_keys=\$3
        local truth_count query_count tp fp fn
        truth_count=\$(wc -l < "\${truth_keys}")
        query_count=\$(wc -l < "\${query_keys}")
        tp=\$(comm -12 "\${truth_keys}" "\${query_keys}" | wc -l)
        fp=\$((query_count-tp))
        fn=\$((truth_count-tp))
        awk -v sample='${output_id}' -v analysis='${analysis_type}' -v caller='${caller_name}' \
            -v name='${benchmark_name}' -v version='${benchmark_version}' -v mode='${benchmark_mode}' \
            -v class="\${variant_class}" -v tp="\${tp}" -v fp="\${fp}" -v fn="\${fn}" \
            -v bases="\${benchmark_region_bases}" -v assembly='${benchmark_assembly}' 'BEGIN {
                precision=(tp+fp)==0 ? 0 : tp/(tp+fp)
                recall=(tp+fn)==0 ? 0 : tp/(tp+fn)
                f1=(precision+recall)==0 ? 0 : 2*precision*recall/(precision+recall)
                OFS="\\t"
                print sample,analysis,caller,name,version,mode,class,tp,tp,fp,fn,
                    sprintf("%.12g",precision),sprintf("%.12g",recall),sprintf("%.12g",f1),bases,
                    "normalized_exact_allele","1.0 (bcftools 1.20)",assembly
            }' >> '${output_id}.${caller_name}.benchmark_summary.tsv'
    }

    write_metrics OVERALL \
        '${output_id}.${caller_name}.benchmark.truth.alleles.tsv' \
        '${output_id}.${caller_name}.benchmark.query.alleles.tsv'
    write_metrics SNV truth.snv.tsv query.snv.tsv
    write_metrics INDEL truth.indel.tsv query.indel.tsv

    printf 'sample_or_patient\tanalysis_type\tcaller\t' > '${output_id}.${caller_name}.benchmark_provenance.tsv'
    head -n 1 '${config_manifest}' >> '${output_id}.${caller_name}.benchmark_provenance.tsv'
    printf '%s\t%s\t%s\t' '${output_id}' '${analysis_type}' '${caller_name}' >> '${output_id}.${caller_name}.benchmark_provenance.tsv'
    tail -n 1 '${config_manifest}' >> '${output_id}.${caller_name}.benchmark_provenance.tsv'

    test -s '${output_id}.${caller_name}.benchmark_summary.tsv'
    test -s '${output_id}.${caller_name}.benchmark_provenance.tsv'
    """
}
