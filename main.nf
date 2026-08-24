nextflow.enable.dsl=2

include { INPUT_CHECK } from './modules/INPUT_CHECK'
include { FASTQC } from './modules/FASTQC'
include { FASTP } from './modules/FASTP'
include { FIXTURE_REFERENCE_PREP; FIXTURE_BWA_INDEX } from './modules/REFERENCE_PREP'
include { BWA_MEM2_ALIGN } from './modules/BWA_MEM2_ALIGN'
include { SAMTOOLS_SORT_INDEX } from './modules/SAMTOOLS_SORT_INDEX'
include { SAMTOOLS_MARKDUP_QC } from './modules/SAMTOOLS_MARKDUP_QC'
include { GATK_BQSR } from './modules/GATK_BQSR'
include { SAMTOOLS_BQSR_INDEX } from './modules/SAMTOOLS_BQSR_SUPPORT'
include { GATK_HAPLOTYPECALLER; GATK_HAPLOTYPECALLER_GVCF; GATK_GENOTYPE_GVCF } from './modules/GATK_HAPLOTYPECALLER'
include { DEEPVARIANT } from './modules/DEEPVARIANT'
include { MUTECT2; MUTECT2_PRODUCTION; LEARN_READ_ORIENTATION_MODEL; GET_PILEUP_SUMMARIES; CALCULATE_CONTAMINATION; FILTERMUTECTCALLS; FILTERMUTECTCALLS_PRODUCTION } from './modules/GATK_MUTECT2'
include { STRELKA2_REGIONS; STRELKA2; MERGE_STRELKA2_VCFS } from './modules/STRELKA2'
include { DEEPSOMATIC } from './modules/DEEPSOMATIC'
include { VCF_PROCESSING } from './modules/VCF_PROCESSING'
include { VEP } from './modules/VEP'
include { VCF2CIRCOS_PREP; VCF2CIRCOS; VCF2CIRCOS_EMPTY_TRACKS; VCF2CIRCOS_GENOME } from './modules/VCF2CIRCOS'
include { MULTIQC } from './modules/MULTIQC'
include { HAPPY_BENCHMARK; EXACT_ALLELE_BENCHMARK } from './modules/BENCHMARK'
include { CALLSET_COMPARISON } from './modules/CALLSET_COMPARISON'
include { STRATIFIED_EXACT_BENCHMARK; STRATIFIED_GERMLINE_COMPARISON; STRATIFIED_SOMATIC_COMPARISON } from './modules/STRATIFICATION'
include { MANTA_REGIONS; MANTA_GERMLINE; MANTA_SOMATIC; SV_SUMMARY } from './modules/MANTA'
include { MOSDEPTH; COVERAGE_SUMMARY } from './modules/MOSDEPTH'
include { CNV_COLLECT_READ_COUNTS; GERMLINE_CONTIG_PLOIDY_CASE; GERMLINE_CNV_CALLER_CASE; POSTPROCESS_GERMLINE_CNV; SOMATIC_CNV_DENOISE; CNV_COLLECT_ALLELIC_COUNTS; SOMATIC_CNV_MODEL_SEGMENTS; SOMATIC_CNV_CALL_SEGMENTS; CNV_SUMMARY } from './modules/GATK_CNV'
include { ANALYSIS_REPORT } from './modules/ANALYSIS_REPORT'

workflow {
    ['run_sv', 'run_cnv', 'run_coverage_qc'].each { option_name ->
        def option_value = params[option_name]?.toString()?.trim()?.toLowerCase() ?: 'false'
        if (!(option_value in ['true', 'false'])) {
            throw new IllegalArgumentException("WGS modality configuration error: --${option_name} must be true or false")
        }
    }
    run_sv = params.run_sv?.toString()?.trim()?.toLowerCase() == 'true'
    run_cnv = params.run_cnv?.toString()?.trim()?.toLowerCase() == 'true'
    run_coverage_qc = params.run_coverage_qc?.toString()?.trim()?.toLowerCase() == 'true'
    coverage_thresholds = params.coverage_thresholds?.toString()?.trim() ?: '10,20,30'
    performance_flag = params.run_performance_report?.toString()?.trim()?.toLowerCase() ?: 'false'
    if (!(performance_flag in ['true', 'false'])) {
        throw new IllegalArgumentException("Performance configuration error: --run_performance_report must be true or false")
    }
    analysis_report_flag = params.run_analysis_report?.toString()?.trim()?.toLowerCase() ?: 'false'
    if (!(analysis_report_flag in ['true', 'false'])) {
        throw new IllegalArgumentException("Integrated report configuration error: --run_analysis_report must be true or false")
    }
    run_analysis_report = analysis_report_flag == 'true'
    if (params.performance_executor?.toString() == 'slurm') {
        def integer_options = [
            hpc_queue_size: params.hpc_queue_size,
            hpc_low_cpus: params.hpc_low_cpus,
            hpc_medium_cpus: params.hpc_medium_cpus,
            hpc_high_cpus: params.hpc_high_cpus
        ]
        integer_options.each { option_name, value ->
            def text = value?.toString()?.trim() ?: ''
            if (!(text ==~ /[1-9][0-9]*/)) {
                throw new IllegalArgumentException("HPC resource configuration error: --${option_name} must be a positive integer")
            }
        }
        def memory_options = [
            hpc_low_memory: params.hpc_low_memory,
            hpc_medium_memory: params.hpc_medium_memory,
            hpc_high_memory: params.hpc_high_memory
        ]
        memory_options.each { option_name, value ->
            def text = value?.toString()?.trim() ?: ''
            def match = text =~ /^([0-9]+(?:\.[0-9]+)?)\s*(B|KB|MB|GB|TB)$/
            if (!match.matches() || new BigDecimal(match.group(1)) <= 0) {
                throw new IllegalArgumentException("HPC resource configuration error: --${option_name} must be a positive Nextflow memory value such as '16 GB'")
            }
        }
        def duration_options = [
            hpc_low_time: params.hpc_low_time,
            hpc_medium_time: params.hpc_medium_time,
            hpc_high_time: params.hpc_high_time
        ]
        duration_options.each { option_name, value ->
            def text = value?.toString()?.trim() ?: ''
            def match = text =~ /^([0-9]+(?:\.[0-9]+)?)\s*(ms|s|m|h|d)$/
            if (!match.matches() || new BigDecimal(match.group(1)) <= 0) {
                throw new IllegalArgumentException("HPC resource configuration error: --${option_name} must be a positive duration such as '24h'")
            }
        }
    }
    resource_mode = params.resource_mode?.toString()?.trim() ?: 'fixture'
    assembly = params.assembly?.toString()?.trim() ?: 'GRCh38'
    if (!(resource_mode in ['fixture', 'production'])) {
        throw new IllegalArgumentException("Resource configuration error: --resource_mode must be 'fixture' or 'production'")
    }
    run_bqsr = params.run_bqsr?.toString()?.toBoolean() ?: false
    run_vep = params.run_vep?.toString()?.toBoolean() ?: false
    run_vcf2circos = params.run_vcf2circos?.toString()?.toBoolean() ?: false
    vcf2circos_mode = params.vcf2circos_mode?.toString()?.trim()?.toLowerCase() ?: 'regional'
    if (!(vcf2circos_mode in ['regional', 'genome'])) {
        throw new IllegalArgumentException("VCF2Circos configuration error: --vcf2circos_mode must be 'regional' or 'genome'")
    }
    def vcf2circos_bin_text = params.vcf2circos_bin_size?.toString()?.trim() ?: '1000000'
    if (!(vcf2circos_bin_text ==~ /[1-9][0-9]*/)) {
        throw new IllegalArgumentException("VCF2Circos configuration error: --vcf2circos_bin_size must be a positive integer")
    }
    vcf2circos_bin_size = vcf2circos_bin_text.toInteger()
    run_multiqc = params.run_multiqc?.toString()?.toBoolean() ?: false
    run_benchmark = params.run_benchmark?.toString()?.toBoolean() ?: false
    selected_caller = params.caller?.toString()?.trim() ?: ''
    gatk_germline_mode = params.gatk_germline_mode?.toString()?.trim() ?: 'vcf'
    gatk_somatic_profile = params.gatk_somatic_profile?.toString()?.trim() ?: 'basic'
    if (gatk_somatic_profile == 'basic' && (params.mutect2_germline_resource || params.mutect2_pon || params.mutect2_contamination_sites)) {
        throw new IllegalArgumentException("Mutect2 resource configuration error: select --gatk_somatic_profile production to use germline, PoN, or contamination resources")
    }
    vep_mode = params.vep_mode?.toString() ?: 'cache'
    input_sheet = file(params.input ?: '')
    input_path = input_sheet.toAbsolutePath().toString()
    resource_dir_path = params.resource_dir ? file(params.resource_dir).toAbsolutePath().toString() : ''
    if (resource_mode == 'fixture') {
        reference_path = file(params.reference ?: "${projectDir}/test_data/reference/GRCh38_chr20.fa.gz").toAbsolutePath().toString()
        reference_fai_path = ''
        reference_dict_path = ''
        bwa_index_prefix = ''
        resource_manifest_path = file(params.resource_manifest ?: "${projectDir}/test_data/reference/resource_manifest.tsv").toAbsolutePath().toString()
    } else {
        def configured_reference = params.reference ?: (resource_dir_path ? "${resource_dir_path}/genome/GRCh38.fa" : '')
        reference_path = configured_reference ? file(configured_reference).toAbsolutePath().toString() : ''
        def reference_stem = reference_path.replaceFirst(/\.(fa|fasta)$/, '')
        reference_fai_path = params.reference_fai ? file(params.reference_fai).toAbsolutePath().toString() : (reference_path ? "${reference_path}.fai" : '')
        reference_dict_path = params.reference_dict ? file(params.reference_dict).toAbsolutePath().toString() : (reference_stem ? "${reference_stem}.dict" : '')
        def reference_name = reference_stem ? file(reference_stem).name : 'GRCh38'
        bwa_index_prefix = params.bwa_index ? file(params.bwa_index).toAbsolutePath().toString() : (resource_dir_path ? "${resource_dir_path}/bwa_mem2/${reference_name}" : '')
        resource_manifest_path = params.resource_manifest ? file(params.resource_manifest).toAbsolutePath().toString() : (resource_dir_path ? "${resource_dir_path}/metadata/manifest.tsv" : '')
    }
    calling_interval = params.calling_interval?.toString()?.trim() ?: ''
    if (run_vcf2circos && vcf2circos_mode == 'genome') {
        if (resource_mode != 'production') {
            throw new IllegalArgumentException("VCF2Circos genome mode requires --resource_mode production and a full reference FAI")
        }
        if (calling_interval) {
            throw new IllegalArgumentException("VCF2Circos genome mode requires unrestricted production calling; remove --calling_interval")
        }
    }
    if (params.known_sites) {
        known_sites_paths = params.known_sites.toString().split(',').collect { known_site -> file(known_site.trim()).toAbsolutePath().toString() }
    } else if (resource_mode == 'fixture') {
        known_sites_paths = [file("${projectDir}/test_data/reference/known_sites/dbsnp138.GRCh38.chr20_fixture_intervals.vcf.gz").toAbsolutePath().toString()]
    } else {
        known_sites_paths = []
    }
    known_sites_joined = known_sites_paths.join('|')
    bqsr_interval = params.bqsr_interval ?: 'auto'
    vep_cache_path = params.vep_cache ? file(params.vep_cache).toAbsolutePath().toString() : (resource_mode == 'production' && resource_dir_path ? "${resource_dir_path}/annotation/vep" : '')
    vep_gtf_path = params.vep_gtf ? file(params.vep_gtf).toAbsolutePath().toString() : ''
    vcf2circos_static_path = params.vcf2circos_static ? file(params.vcf2circos_static).toAbsolutePath().toString() : ''
    mutect2_germline_resource_path = params.mutect2_germline_resource \
        ? file(params.mutect2_germline_resource).toAbsolutePath().toString() \
        : (resource_dir_path ? "${resource_dir_path}/somatic/germline_resource.vcf.gz" : '')
    mutect2_pon_path = params.mutect2_pon \
        ? file(params.mutect2_pon).toAbsolutePath().toString() \
        : (resource_dir_path ? "${resource_dir_path}/somatic/panel_of_normals.vcf.gz" : '')
    mutect2_contamination_sites_path = params.mutect2_contamination_sites \
        ? file(params.mutect2_contamination_sites).toAbsolutePath().toString() \
        : (resource_dir_path ? "${resource_dir_path}/somatic/contamination_sites.vcf.gz" : '')
    benchmark_mode = params.benchmark_mode?.toString()?.trim() ?: 'regression'
    benchmark_sex = params.benchmark_sex?.toString()?.trim()?.toLowerCase() ?: 'auto'
    benchmark_assembly = params.benchmark_assembly?.toString()?.trim() ?: assembly
    if (resource_mode == 'fixture') {
        if ((params.analysis_type ?: '') == 'germline') {
            default_benchmark_truth = "${projectDir}/test_data/reference/benchmarks/HG002_GRCh38_giab_v4.2.1_chr20_12446968_12455543.bgzip.vcf.gz"
            default_benchmark_regions = "${projectDir}/test_data/reference/germline_test_interval.bed"
            default_benchmark_name = 'HG002_GIAB'
            default_benchmark_version = 'v4.2.1'
        } else {
            default_benchmark_truth = "${projectDir}/test_data/reference/benchmarks/HG008_GRCh38_somatic_smvar_v0.3_chr20_21900000_21980636.bgzip.vcf.gz"
            default_benchmark_regions = "${projectDir}/test_data/reference/somatic_test_interval.bed"
            default_benchmark_name = 'HG008_T_somatic_small_variant_truth'
            default_benchmark_version = 'v0.3'
        }
    } else {
        default_benchmark_truth = resource_dir_path ? "${resource_dir_path}/benchmarking/truth.vcf.gz" : ''
        default_benchmark_regions = resource_dir_path ? "${resource_dir_path}/benchmarking/confident_regions.bed" : ''
        default_benchmark_name = ''
        default_benchmark_version = ''
    }
    benchmark_truth_path = params.benchmark_truth \
        ? file(params.benchmark_truth).toAbsolutePath().toString() \
        : (run_benchmark ? file(default_benchmark_truth).toAbsolutePath().toString() : '')
    benchmark_truth_index_path = params.benchmark_truth_index \
        ? file(params.benchmark_truth_index).toAbsolutePath().toString() \
        : (benchmark_truth_path ? "${benchmark_truth_path}.tbi" : '')
    benchmark_regions_path = params.benchmark_regions \
        ? file(params.benchmark_regions).toAbsolutePath().toString() \
        : (run_benchmark ? file(default_benchmark_regions).toAbsolutePath().toString() : '')
    benchmark_name = params.benchmark_name?.toString()?.trim() ?: default_benchmark_name
    benchmark_version = params.benchmark_version?.toString()?.trim() ?: default_benchmark_version
    benchmark_stratifications_path = params.benchmark_stratifications \
        ? file(params.benchmark_stratifications).toAbsolutePath().toString() \
        : ''
    stratification_resource_paths = []
    stratification_tsv_name = ''
    if (benchmark_stratifications_path) {
        def stratification_tsv_file = file(benchmark_stratifications_path, checkIfExists: false)
        stratification_tsv_name = stratification_tsv_file.name
        stratification_resource_paths << stratification_tsv_file
        if (stratification_tsv_file.exists()) {
            stratification_tsv_file.eachLine { line ->
                def stripped = line.trim()
                if (stripped && !stripped.startsWith('#')) {
                    def fields = line.split('\t', -1)
                    if (fields.size() == 2 && fields[1].trim()) {
                        def configured_bed = fields[1].trim()
                        def resolved_bed = configured_bed.startsWith('/') \
                            ? file(configured_bed, checkIfExists: false) \
                            : file("${stratification_tsv_file.parent}/${configured_bed}", checkIfExists: false)
                        stratification_resource_paths << resolved_bed
                    }
                }
            }
        }
    }
    benchmark_engine_container = benchmark_mode == 'germline_happy' \
        ? 'mgibio/hap.py@sha256:492574403bf5cc648cf9c297957cc2eac60cbbc7da1a968984bafad060536bdc' \
        : 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    ploidy_config_path = params.ploidy_config ? file(params.ploidy_config).toAbsolutePath().toString() : ''
    cnv_intervals_path = params.cnv_intervals \
        ? file(params.cnv_intervals).toAbsolutePath().toString() \
        : (resource_dir_path ? "${resource_dir_path}/cnv/intervals.interval_list" : '')
    cnv_ploidy_model_path = params.cnv_ploidy_model \
        ? file(params.cnv_ploidy_model).toAbsolutePath().toString() \
        : (resource_dir_path ? "${resource_dir_path}/cnv/germline/contig-ploidy-model" : '')
    cnv_germline_model_path = params.cnv_germline_model \
        ? file(params.cnv_germline_model).toAbsolutePath().toString() \
        : (resource_dir_path ? "${resource_dir_path}/cnv/germline/gcnv-model" : '')
    cnv_read_count_pon_path = params.cnv_read_count_pon \
        ? file(params.cnv_read_count_pon).toAbsolutePath().toString() \
        : (resource_dir_path ? "${resource_dir_path}/cnv/somatic/read-count-pon.hdf5" : '')
    cnv_common_sites_path = params.cnv_common_sites \
        ? file(params.cnv_common_sites).toAbsolutePath().toString() \
        : (resource_dir_path ? "${resource_dir_path}/cnv/somatic/common-sites.vcf.gz" : '')

    validation = INPUT_CHECK(
        input_path,
        reference_path,
        params.analysis_type ?: '',
        params.caller ?: '',
        params.outdir ?: '',
        run_bqsr,
        bqsr_interval,
        calling_interval,
        resource_mode,
        assembly,
        projectDir.toAbsolutePath().toString(),
        reference_fai_path,
        reference_dict_path,
        bwa_index_prefix,
        known_sites_joined,
        resource_manifest_path,
        run_vep,
        vep_mode,
        vep_cache_path,
        vep_gtf_path,
        run_vcf2circos && vcf2circos_mode == 'regional',
        vcf2circos_static_path,
        gatk_germline_mode,
        gatk_somatic_profile,
        mutect2_germline_resource_path,
        mutect2_pon_path,
        mutect2_contamination_sites_path,
        run_benchmark,
        benchmark_mode,
        benchmark_truth_path,
        benchmark_truth_index_path,
        benchmark_regions_path,
        benchmark_assembly,
        benchmark_name,
        benchmark_version,
        benchmark_sex,
        benchmark_stratifications_path,
        benchmark_engine_container,
        run_sv,
        run_cnv,
        run_coverage_qc,
        coverage_thresholds,
        ploidy_config_path,
        cnv_intervals_path,
        cnv_ploidy_model_path,
        cnv_germline_model_path,
        cnv_read_count_pon_path,
        cnv_common_sites_path
    )

    samplesheet = input_sheet
    sheet_dir = samplesheet.parent
    def report_lines = samplesheet.readLines().findAll { line -> line.trim() && !line.startsWith('#') }
    def report_header = report_lines ? report_lines[0].split('\t', -1).toList() : []
    def patient_index = report_header.indexOf('patient_id')
    def sample_index = report_header.indexOf('sample_id')
    def report_fields = report_lines.size() > 1 ? report_lines[1].split('\t', -1) : []
    report_identifier = patient_index >= 0 && patient_index < report_fields.size() && report_fields[patient_index].trim() \
        ? report_fields[patient_index].trim() \
        : (sample_index >= 0 && sample_index < report_fields.size() ? report_fields[sample_index].trim() : 'analysis')

    paired_reads = channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample_meta = [
                sample_id: row.sample_id.toString(),
                patient_id: row.patient_id.toString(),
                role: row.role.toString()
            ]
            def read1_value = row.fastq_r1.toString()
            def read2_value = row.fastq_r2.toString()
            def read1_path = read1_value.startsWith('/') ? file(read1_value) : file("${sheet_dir}/${read1_value}")
            def read2_path = read2_value.startsWith('/') ? file(read2_value) : file("${sheet_dir}/${read2_value}")
            tuple(sample_meta, read1_path, read2_path)
        }

    fastq_reads = paired_reads.flatMap { sample_meta, read1, read2 ->
        [tuple(sample_meta, read1), tuple(sample_meta, read2)]
    }

    validation_gate = validation.passed.map { true }
    fastqc_reports = FASTQC(fastq_reads, validation_gate).reports
    fastp_results = FASTP(paired_reads, validation_gate).cleaned
    cleaned_reads = fastp_results
    report_artifacts = validation.resource_manifest
        .mix(validation.gatk_path_manifest, validation.benchmark_config_manifest, validation.wgs_modality_provenance)
        .mix(fastqc_reports.map { _meta, _html, archive -> archive })
        .mix(fastp_results.map { _meta, _read1, _read2, _html, summary_json -> summary_json })

    if (resource_mode == 'fixture') {
        fixture_reference_source = channel.fromPath(reference_path, checkIfExists: false)
        gated_fixture_reference = fixture_reference_source.combine(validation_gate).map { fixture_source, _passed -> fixture_source }
        prepared_fixture_reference = FIXTURE_REFERENCE_PREP(gated_fixture_reference).reference
        reference_bundle = FIXTURE_BWA_INDEX(prepared_fixture_reference).bundle
    } else {
        def bwa_suffixes = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']
        def production_index_files = bwa_suffixes.collect { suffix -> file("${bwa_index_prefix}${suffix}", checkIfExists: false) }
        def bwa_alt_file = file("${bwa_index_prefix}.alt", checkIfExists: false)
        if (bwa_alt_file.exists()) production_index_files << bwa_alt_file
        production_reference_values = channel.value(tuple(
            file(reference_path, checkIfExists: false),
            file(reference_fai_path, checkIfExists: false),
            file(reference_dict_path, checkIfExists: false),
            production_index_files,
            file(bwa_index_prefix).name
        ))
        reference_bundle = production_reference_values.combine(validation_gate).map { reference_fa, reference_fai, reference_dict, bwa_files, bwa_prefix_name, _passed ->
            tuple(reference_fa, reference_fai, reference_dict, bwa_files, bwa_prefix_name)
        }
    }

    reference = reference_bundle.map { reference_fa, reference_fai, reference_dict, _bwa_files, _bwa_prefix_name ->
        tuple(reference_fa, reference_fai, reference_dict)
    }
    bwa_reference = reference_bundle.map { reference_fa, _reference_fai, _reference_dict, bwa_files, bwa_prefix_name ->
        tuple(reference_fa, bwa_files, bwa_prefix_name)
    }
    caller_calls = channel.empty()
    alignment_reads = cleaned_reads.map { meta, read1, read2, _fastp_html, _fastp_json ->
        tuple(meta, read1, read2)
    }
    alignment_inputs = alignment_reads.combine(bwa_reference)
    aligned_reads = BWA_MEM2_ALIGN(alignment_inputs).aligned
    sort_index_results = SAMTOOLS_SORT_INDEX(aligned_reads)
    sorted_reads = sort_index_results.aligned
    markdup_results = SAMTOOLS_MARKDUP_QC(sorted_reads)
    duplicate_marked = markdup_results.marked

    if (run_bqsr) {
        known_sites_resources = known_sites_paths.collectMany { known_path ->
            [file(known_path, checkIfExists: false), file("${known_path}.tbi", checkIfExists: false)]
        }
        bqsr_inputs = duplicate_marked.combine(reference).map { meta, markdup_bam, markdup_bai, markdup_stats, reference_fa, reference_fai, reference_dict ->
            tuple(
                meta,
                markdup_bam,
                markdup_bai,
                markdup_stats,
                reference_fa,
                reference_fai,
                reference_dict,
                known_sites_resources,
                bqsr_interval == 'auto' ? calling_interval : bqsr_interval
            )
        }
        recalibrated = GATK_BQSR(bqsr_inputs).recalibrated
        indexed_recalibrated = SAMTOOLS_BQSR_INDEX(recalibrated)
        caller_ready_bam = indexed_recalibrated.caller_ready.map { meta, recal_bam, recal_bai, _recal_table ->
            tuple(meta, recal_bam, recal_bai)
        }
    } else {
        caller_ready_bam = duplicate_marked.map { meta, markdup_bam, markdup_bai, _markdup_stats ->
            tuple(meta, markdup_bam, markdup_bai)
        }
    }

    if (run_coverage_qc) {
        def reporting_thresholds = ([1] + coverage_thresholds.split(',').collect { threshold -> threshold.trim().toInteger() }).unique().sort().join(',')
        mosdepth_inputs = caller_ready_bam.combine(reference).map { meta, caller_bam, caller_bai, _reference_fa, reference_fai, _reference_dict ->
            tuple(meta, caller_bam, caller_bai, reference_fai, calling_interval, reporting_thresholds)
        }
        mosdepth_results = MOSDEPTH(mosdepth_inputs)
        coverage_summary_inputs = mosdepth_results.coverage.map { meta, summary, global_dist, region_dist, regions, regions_index, threshold_bed, coverage_regions ->
            tuple(meta, summary, global_dist, region_dist, regions, regions_index, threshold_bed, coverage_regions, params.analysis_type.toString(), assembly)
        }
        coverage_summaries = COVERAGE_SUMMARY(coverage_summary_inputs).summaries
        report_artifacts = report_artifacts.mix(coverage_summaries.flatMap { _meta, summary, thresholds, chromosomes -> [summary, thresholds, chromosomes] })
    }

    if (run_sv) {
        manta_region_channel = channel.of(calling_interval).combine(validation_gate).map { configured_interval, _passed -> configured_interval }
        manta_regions = MANTA_REGIONS(manta_region_channel).regions
        manta_resources = reference.combine(manta_regions)
        if ((params.analysis_type ?: '') == 'germline') {
            manta_inputs = caller_ready_bam.combine(manta_resources).map { meta, bam, bai, reference_fa, reference_fai, reference_dict, interval, region_dir ->
                tuple(meta, bam, bai, reference_fa, reference_fai, reference_dict, interval, region_dir)
            }
            manta_calls = MANTA_GERMLINE(manta_inputs).calls
            manta_summary_inputs = manta_calls.combine(reference).map { meta, semantics, final_vcf, final_tbi, _candidate_vcf, _candidate_tbi, _reference_fa, reference_fai, _reference_dict ->
                tuple(meta, semantics, final_vcf, final_tbi, reference_fai, 'germline', assembly)
            }
        } else {
            manta_pairs = caller_ready_bam.map { meta, bam, bai -> tuple(meta.patient_id, meta, bam, bai) }.groupTuple().map { patient_id, metas, bams, bais ->
                def tumor_index = metas.findIndexOf { sample_meta -> sample_meta.role == 'tumor' }
                def normal_index = metas.findIndexOf { sample_meta -> sample_meta.role == 'normal' }
                if (tumor_index < 0 || normal_index < 0 || metas.size() != 2) throw new IllegalArgumentException("Manta pairing error for patient '${patient_id}': expected exactly one tumor and one normal")
                tuple([patient_id: patient_id, tumor_sample_id: metas[tumor_index].sample_id, normal_sample_id: metas[normal_index].sample_id], bams[tumor_index], bais[tumor_index], bams[normal_index], bais[normal_index])
            }
            manta_inputs = manta_pairs.combine(manta_resources).map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, reference_fa, reference_fai, reference_dict, interval, region_dir ->
                tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, reference_fa, reference_fai, reference_dict, interval, region_dir)
            }
            manta_calls = MANTA_SOMATIC(manta_inputs).calls
            manta_summary_inputs = manta_calls.combine(reference).map { meta, semantics, final_vcf, final_tbi, _normal_vcf, _normal_tbi, _candidate_vcf, _candidate_tbi, _reference_fa, reference_fai, _reference_dict ->
                tuple(meta, semantics, final_vcf, final_tbi, reference_fai, 'somatic', assembly)
            }
        }
        sv_summaries = SV_SUMMARY(manta_summary_inputs).summaries
        report_artifacts = report_artifacts.mix(sv_summaries.flatMap { _meta, summary, classes, sizes, chromosomes -> [summary, classes, sizes, chromosomes] })
    }

    if (run_cnv) {
        cnv_intervals = file(cnv_intervals_path, checkIfExists: false)
        cnv_counts = CNV_COLLECT_READ_COUNTS(caller_ready_bam.map { meta, bam, bai -> tuple(meta, bam, bai, cnv_intervals) }).counts
        if ((params.analysis_type ?: '') == 'germline') {
            ploidy_model = file(cnv_ploidy_model_path, checkIfExists: false)
            gcnv_model = file(cnv_germline_model_path, checkIfExists: false)
            ploidy_calls = GERMLINE_CONTIG_PLOIDY_CASE(cnv_counts.map { meta, counts -> tuple(meta, counts, ploidy_model) }).calls
            gcnv_calls = GERMLINE_CNV_CALLER_CASE(ploidy_calls.map { meta, counts, calls -> tuple(meta, counts, calls, gcnv_model) }).calls
            allosomal_contigs = []
            if (reference_fai_path && file(reference_fai_path, checkIfExists: false).exists()) {
                def contig_names = file(reference_fai_path).readLines().collect { reference_line -> reference_line.split('\t')[0] }
                allosomal_contigs = ['chrX', 'chrY', 'X', 'Y'].findAll { contig_name -> contig_name in contig_names }
            }
            germline_cnv_final = POSTPROCESS_GERMLINE_CNV(gcnv_calls.map { meta, counts, ploidy, calls -> tuple(meta, counts, ploidy, calls, gcnv_model, allosomal_contigs) }).results
            cnv_summary_inputs = germline_cnv_final.map { meta, _interval_vcf, segment_vcf, _ratios -> tuple(meta, segment_vcf, 'germline_vcf', 'germline', assembly) }
        } else {
            read_count_pon = file(cnv_read_count_pon_path, checkIfExists: false)
            common_sites = file(cnv_common_sites_path, checkIfExists: false)
            common_sites_tbi = file("${cnv_common_sites_path}.tbi", checkIfExists: false)
            tumor_ratios = SOMATIC_CNV_DENOISE(cnv_counts.filter { meta, _counts -> meta.role == 'tumor' }.map { meta, counts -> tuple(meta, counts, read_count_pon) }).ratios
            allelic_counts = CNV_COLLECT_ALLELIC_COUNTS(caller_ready_bam.combine(reference).map { meta, bam, bai, reference_fa, reference_fai, reference_dict -> tuple(meta, bam, bai, reference_fa, reference_fai, reference_dict, common_sites, common_sites_tbi) }).counts
            paired_allelic_counts = allelic_counts.map { meta, counts -> tuple(meta.patient_id, meta, counts) }.groupTuple()
            keyed_tumor_ratios = tumor_ratios.map { meta, _standardized, denoised -> tuple(meta.patient_id, meta, denoised) }
            somatic_model_inputs = keyed_tumor_ratios.join(paired_allelic_counts).map { _patient, tumor_meta, denoised, metas, counts ->
                def tumor_index = metas.findIndexOf { sample_meta -> sample_meta.role == 'tumor' }
                def normal_index = metas.findIndexOf { sample_meta -> sample_meta.role == 'normal' }
                tuple([patient_id: tumor_meta.patient_id, tumor_sample_id: tumor_meta.sample_id, normal_sample_id: metas[normal_index].sample_id], denoised, counts[tumor_index], counts[normal_index])
            }
            somatic_cnv_final = SOMATIC_CNV_CALL_SEGMENTS(SOMATIC_CNV_MODEL_SEGMENTS(somatic_model_inputs).segments).results
            cnv_summary_inputs = somatic_cnv_final.map { meta, _cr, _model, _hets, called -> tuple(meta, called, 'somatic_seg', 'somatic', assembly) }
        }
        cnv_summaries = CNV_SUMMARY(cnv_summary_inputs).summary
        report_artifacts = report_artifacts.mix(cnv_summaries.map { _meta, summary -> summary })
    }

    if ((params.analysis_type ?: '') == 'germline' && selected_caller in ['haplotypecaller', 'compare']) {
        haplotypecaller_inputs = caller_ready_bam.combine(reference).map { meta, caller_bam, caller_bai, reference_fa, reference_fai, reference_dict ->
            tuple(
                meta,
                'haplotypecaller',
                caller_bam,
                caller_bai,
                reference_fa,
                reference_fai,
                reference_dict,
                calling_interval
            )
        }
        if (gatk_germline_mode == 'vcf') {
            haplotypecaller_calls = GATK_HAPLOTYPECALLER(haplotypecaller_inputs).calls
        } else {
            haplotypecaller_gvcfs = GATK_HAPLOTYPECALLER_GVCF(haplotypecaller_inputs).gvcfs
            genotype_gvcf_inputs = haplotypecaller_gvcfs.combine(reference).map { meta, caller_name, gvcf, gvcf_tbi, reference_fa, reference_fai, reference_dict ->
                tuple(meta, caller_name, gvcf, gvcf_tbi, reference_fa, reference_fai, reference_dict, calling_interval)
            }
            haplotypecaller_calls = GATK_GENOTYPE_GVCF(genotype_gvcf_inputs).calls
        }
        if (selected_caller == 'haplotypecaller') caller_calls = haplotypecaller_calls
    }

    if ((params.analysis_type ?: '') == 'germline' && selected_caller in ['deepvariant', 'compare']) {
        deepvariant_inputs = caller_ready_bam.combine(reference).map { meta, caller_bam, caller_bai, reference_fa, reference_fai, reference_dict ->
            tuple(
                meta,
                'deepvariant',
                caller_bam,
                caller_bai,
                reference_fa,
                reference_fai,
                reference_dict,
                calling_interval
            )
        }
        deepvariant_calls = DEEPVARIANT(deepvariant_inputs).calls
        if (selected_caller == 'deepvariant') caller_calls = deepvariant_calls
    }

    if ((params.analysis_type ?: '') == 'somatic' && selected_caller in ['mutect2', 'compare']) {
        mutect2_reference = reference
        somatic_pairs = caller_ready_bam
            .map { meta, caller_bam, caller_bai -> tuple(meta.patient_id, meta, caller_bam, caller_bai) }
            .groupTuple()
            .map { patient_id, metas, bams, bais ->
                def tumor_index = metas.findIndexOf { sample_meta -> sample_meta.role == 'tumor' }
                def normal_index = metas.findIndexOf { sample_meta -> sample_meta.role == 'normal' }
                if (tumor_index < 0 || normal_index < 0 || tumor_index == normal_index || metas.size() != 2) {
                    throw new IllegalArgumentException("MUTECT2 pairing error for patient '${patient_id}': expected exactly one tumor and one normal")
                }
                def tumor_meta = metas[tumor_index]
                def normal_meta = metas[normal_index]
                tuple(
                    [
                        patient_id: patient_id,
                        tumor_sample_id: tumor_meta.sample_id,
                        normal_sample_id: normal_meta.sample_id,
                        tumor_role: tumor_meta.role,
                        normal_role: normal_meta.role
                    ],
                    bams[tumor_index],
                    bais[tumor_index],
                    bams[normal_index],
                    bais[normal_index],
                    tumor_meta.sample_id,
                    normal_meta.sample_id
                )
            }
        if (gatk_somatic_profile == 'basic') {
            mutect2_inputs = somatic_pairs.combine(mutect2_reference).map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, expected_tumor_sample, expected_normal_sample, reference_fa, reference_fai, reference_dict ->
                tuple(
                    meta,
                    tumor_bam,
                    tumor_bai,
                    normal_bam,
                    normal_bai,
                    reference_fa,
                    reference_fai,
                    reference_dict,
                    expected_tumor_sample,
                    expected_normal_sample,
                    calling_interval
                )
            }
            mutect2_unfiltered = MUTECT2(mutect2_inputs).unfiltered
            filtermutect_inputs = mutect2_unfiltered.combine(mutect2_reference).map { meta, caller_name, unfiltered_vcf, unfiltered_tbi, unfiltered_stats, reference_fa, reference_fai, reference_dict ->
                tuple(meta, caller_name, unfiltered_vcf, unfiltered_tbi, unfiltered_stats, reference_fa, reference_fai, reference_dict)
            }
            mutect2_calls = FILTERMUTECTCALLS(filtermutect_inputs).calls
        } else {
            def germline_resource = file(mutect2_germline_resource_path, checkIfExists: false)
            def germline_resource_tbi = file("${mutect2_germline_resource_path}.tbi", checkIfExists: false)
            def panel_of_normals = file(mutect2_pon_path, checkIfExists: false)
            def panel_of_normals_tbi = file("${mutect2_pon_path}.tbi", checkIfExists: false)
            def contamination_sites = file(mutect2_contamination_sites_path, checkIfExists: false)
            def contamination_sites_tbi = file("${mutect2_contamination_sites_path}.tbi", checkIfExists: false)

            mutect2_production_inputs = somatic_pairs.combine(mutect2_reference).map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, expected_tumor_sample, expected_normal_sample, reference_fa, reference_fai, reference_dict ->
                tuple(
                    meta, tumor_bam, tumor_bai, normal_bam, normal_bai,
                    reference_fa, reference_fai, reference_dict,
                    expected_tumor_sample, expected_normal_sample, calling_interval,
                    germline_resource, germline_resource_tbi, panel_of_normals, panel_of_normals_tbi
                )
            }
            mutect2_production_results = MUTECT2_PRODUCTION(mutect2_production_inputs)
            orientation_priors = LEARN_READ_ORIENTATION_MODEL(mutect2_production_results.f1r2).priors

            pileup_inputs = somatic_pairs.combine(mutect2_reference).map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, _expected_tumor_sample, _expected_normal_sample, reference_fa, reference_fai, reference_dict ->
                tuple(
                    meta, tumor_bam, tumor_bai, normal_bam, normal_bai,
                    reference_fa, reference_fai, reference_dict, calling_interval,
                    contamination_sites, contamination_sites_tbi
                )
            }
            pileups = GET_PILEUP_SUMMARIES(pileup_inputs).pileups
            contamination_estimates = CALCULATE_CONTAMINATION(pileups).estimates

            keyed_unfiltered = mutect2_production_results.unfiltered.map { meta, caller_name, unfiltered_vcf, unfiltered_tbi, unfiltered_stats ->
                tuple(meta.patient_id, meta, caller_name, unfiltered_vcf, unfiltered_tbi, unfiltered_stats)
            }
            keyed_orientation = orientation_priors.map { meta, priors -> tuple(meta.patient_id, priors) }
            keyed_contamination = contamination_estimates.map { meta, contamination, segmentation ->
                tuple(meta.patient_id, contamination, segmentation)
            }
            production_filter_resources = keyed_unfiltered
                .join(keyed_orientation)
                .join(keyed_contamination)
            production_filter_inputs = production_filter_resources.combine(mutect2_reference).map { _patient_id, meta, caller_name, unfiltered_vcf, unfiltered_tbi, unfiltered_stats, priors, contamination, segmentation, reference_fa, reference_fai, reference_dict ->
                tuple(
                    meta, caller_name, unfiltered_vcf, unfiltered_tbi, unfiltered_stats,
                    reference_fa, reference_fai, reference_dict,
                    priors, contamination, segmentation
                )
            }
            mutect2_calls = FILTERMUTECTCALLS_PRODUCTION(production_filter_inputs).calls
        }
        if (selected_caller == 'mutect2') caller_calls = mutect2_calls
    }

    if ((params.analysis_type ?: '') == 'somatic' && selected_caller in ['strelka2', 'compare']) {
        gated_strelka2_interval = channel.of(calling_interval)
            .combine(validation_gate)
            .map { configured_interval, _passed -> configured_interval }
        strelka2_regions = STRELKA2_REGIONS(gated_strelka2_interval).indexed_regions
        strelka2_resources = reference.combine(strelka2_regions)
        strelka2_pairs = caller_ready_bam
            .map { meta, caller_bam, caller_bai -> tuple(meta.patient_id, meta, caller_bam, caller_bai) }
            .groupTuple()
            .map { patient_id, metas, bams, bais ->
                def tumor_index = metas.findIndexOf { sample_meta -> sample_meta.role == 'tumor' }
                def normal_index = metas.findIndexOf { sample_meta -> sample_meta.role == 'normal' }
                if (tumor_index < 0 || normal_index < 0 || tumor_index == normal_index || metas.size() != 2) {
                    throw new IllegalArgumentException("STRELKA2 pairing error for patient '${patient_id}': expected exactly one tumor and one normal")
                }
                def tumor_meta = metas[tumor_index]
                def normal_meta = metas[normal_index]
                tuple(
                    [
                        patient_id: patient_id,
                        tumor_sample_id: tumor_meta.sample_id,
                        normal_sample_id: normal_meta.sample_id,
                        tumor_role: tumor_meta.role,
                        normal_role: normal_meta.role
                    ],
                    bams[tumor_index],
                    bais[tumor_index],
                    bams[normal_index],
                    bais[normal_index],
                    tumor_meta.sample_id,
                    normal_meta.sample_id
                )
            }
        strelka2_inputs = strelka2_pairs.combine(strelka2_resources).map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, expected_tumor_sample, expected_normal_sample, reference_fa, reference_fai, reference_dict, configured_interval, call_regions_dir ->
            tuple(
                meta,
                tumor_bam,
                tumor_bai,
                normal_bam,
                normal_bai,
                reference_fa,
                reference_fai,
                reference_dict,
                configured_interval,
                call_regions_dir,
                expected_tumor_sample,
                expected_normal_sample
            )
        }
        strelka2_raw = STRELKA2(strelka2_inputs).raw_calls
        strelka2_calls = MERGE_STRELKA2_VCFS(strelka2_raw).calls
        if (selected_caller == 'strelka2') caller_calls = strelka2_calls
    }

    if ((params.analysis_type ?: '') == 'somatic' && selected_caller in ['deepsomatic', 'compare']) {
        deepsomatic_pairs = caller_ready_bam
            .map { meta, caller_bam, caller_bai -> tuple(meta.patient_id, meta, caller_bam, caller_bai) }
            .groupTuple()
            .map { patient_id, metas, bams, bais ->
                def tumor_index = metas.findIndexOf { sample_meta -> sample_meta.role == 'tumor' }
                def normal_index = metas.findIndexOf { sample_meta -> sample_meta.role == 'normal' }
                if (tumor_index < 0 || normal_index < 0 || tumor_index == normal_index || metas.size() != 2) {
                    throw new IllegalArgumentException("DEEPSOMATIC pairing error for patient '${patient_id}': expected exactly one tumor and one normal")
                }
                def tumor_meta = metas[tumor_index]
                def normal_meta = metas[normal_index]
                tuple(
                    [
                        patient_id: patient_id,
                        tumor_sample_id: tumor_meta.sample_id,
                        normal_sample_id: normal_meta.sample_id,
                        tumor_role: tumor_meta.role,
                        normal_role: normal_meta.role
                    ],
                    bams[tumor_index],
                    bais[tumor_index],
                    bams[normal_index],
                    bais[normal_index],
                    tumor_meta.sample_id,
                    normal_meta.sample_id
                )
            }
        deepsomatic_inputs = deepsomatic_pairs.combine(reference).map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, expected_tumor_sample, expected_normal_sample, reference_fa, reference_fai, reference_dict ->
            tuple(
                meta,
                tumor_bam,
                tumor_bai,
                normal_bam,
                normal_bai,
                reference_fa,
                reference_fai,
                reference_dict,
                expected_tumor_sample,
                expected_normal_sample,
                calling_interval
            )
        }
        deepsomatic_calls = DEEPSOMATIC(deepsomatic_inputs).calls
        if (selected_caller == 'deepsomatic') caller_calls = deepsomatic_calls
    }

    if (selected_caller == 'compare' && (params.analysis_type ?: '') in ['germline', 'somatic']) {
        if ((params.analysis_type ?: '') == 'germline') {
            caller_calls = haplotypecaller_calls.mix(deepvariant_calls)
        } else if ((params.analysis_type ?: '') == 'somatic') {
            caller_calls = mutect2_calls.mix(strelka2_calls, deepsomatic_calls)
        }
    }

    vcf_processing_inputs = caller_calls.combine(reference).map { meta, caller_name, native_vcf, native_vcf_index, reference_fa, reference_fai, reference_dict ->
        tuple(
            meta,
            caller_name,
            native_vcf,
            native_vcf_index,
            reference_fa,
            reference_fai,
            reference_dict,
            params.analysis_type.toString()
        )
    }
    processed_calls = VCF_PROCESSING(vcf_processing_inputs).processed
    report_artifacts = report_artifacts.mix(processed_calls.map { _meta, _caller_name, _vcf, _index, _stats, summary -> summary })

    if (selected_caller == 'compare' && (params.analysis_type ?: '') in ['germline', 'somatic']) {
        def comparison_callers = (params.analysis_type ?: '') == 'germline' \
            ? ['haplotypecaller', 'deepvariant'] \
            : ['mutect2', 'strelka2', 'deepsomatic']
        grouped_comparison_calls = processed_calls
            .map { meta, caller_name, normalized_vcf, normalized_index, _stats, _summary ->
                def comparison_id = params.analysis_type.toString() == 'germline' ? meta.sample_id : meta.patient_id
                tuple(comparison_id, meta, caller_name, normalized_vcf, normalized_index)
            }
            .groupTuple()
            .map { comparison_id, metas, caller_names, normalized_vcfs, normalized_indexes ->
                if (caller_names.size() != comparison_callers.size() || caller_names.toSet() != comparison_callers.toSet()) {
                    throw new IllegalArgumentException("Comparison routing error for '${comparison_id}': expected ${comparison_callers.join(',')}, observed ${caller_names.join(',')}")
                }
                def ordered_vcfs = comparison_callers.collect { caller_name -> normalized_vcfs[caller_names.indexOf(caller_name)] }
                def ordered_indexes = comparison_callers.collect { caller_name -> normalized_indexes[caller_names.indexOf(caller_name)] }
                def comparison_meta = metas[caller_names.indexOf(comparison_callers[0])] + [comparison_id: comparison_id]
                tuple(comparison_meta, comparison_callers, ordered_vcfs, ordered_indexes)
            }

        comparison_context = reference
            .combine(validation.passed)
            .map { reference_fa, reference_fai, reference_dict, validation_sentinel ->
                def staged_resources = [validation_sentinel]
                def configured_truth_name = ''
                def configured_truth_index_name = ''
                def configured_regions_name = ''
                if (run_benchmark) {
                    def truth_resource = file(benchmark_truth_path, checkIfExists: false)
                    def truth_index_resource = file(benchmark_truth_index_path, checkIfExists: false)
                    def regions_resource = file(benchmark_regions_path, checkIfExists: false)
                    staged_resources.addAll([truth_resource, truth_index_resource, regions_resource])
                    configured_truth_name = truth_resource.name
                    configured_truth_index_name = truth_index_resource.name
                    configured_regions_name = regions_resource.name
                }
                tuple(
                    reference_fa, reference_fai, reference_dict, staged_resources,
                    params.analysis_type.toString(), assembly, run_bqsr.toString(), run_benchmark.toString(),
                    run_benchmark ? benchmark_name : '', run_benchmark ? benchmark_version : '',
                    configured_truth_name, configured_truth_index_name, configured_regions_name
                )
            }
        comparison_inputs = grouped_comparison_calls.combine(comparison_context).map { meta, callers, normalized_vcfs, normalized_indexes, reference_fa, reference_fai, reference_dict, comparison_resources, configured_analysis_type, configured_assembly, configured_bqsr, configured_benchmark, configured_benchmark_name, configured_benchmark_version, truth_name, truth_index_name, regions_name ->
            tuple(
                meta, callers, normalized_vcfs, normalized_indexes,
                reference_fa, reference_fai, reference_dict, comparison_resources,
                configured_analysis_type, configured_assembly, configured_bqsr, configured_benchmark,
                configured_benchmark_name, configured_benchmark_version,
                truth_name, truth_index_name, regions_name
            )
        }
        CALLSET_COMPARISON(comparison_inputs)
        report_artifacts = report_artifacts.mix(CALLSET_COMPARISON.out.results.flatMap { _meta, summary, pairwise, matrix, membership, discordant, provenance, _vcfs -> [summary, pairwise, matrix, membership, discordant, provenance] })
        if ((params.analysis_type ?: '') == 'germline') {
            report_artifacts = report_artifacts.mix(CALLSET_COMPARISON.out.genotype_summary)
        }
        if (run_benchmark && stratification_tsv_name) {
            comparison_base = CALLSET_COMPARISON.out.results.map { meta, _summary, _jaccard, _matrix, membership, discordant, _provenance, _vcfs ->
                tuple(meta.comparison_id, meta, membership, discordant)
            }
            if ((params.analysis_type ?: '') == 'germline') {
                comparison_genotypes = CALLSET_COMPARISON.out.genotype_details.map { meta, genotype_details ->
                    tuple(meta.comparison_id, genotype_details)
                }
                stratified_comparison_source = comparison_base
                    .join(comparison_genotypes)
                    .map { _comparison_id, meta, membership, discordant, genotype_details ->
                        tuple(meta, membership, discordant, genotype_details)
                    }
                stratified_comparison_context = channel.value(tuple(
                    file(benchmark_regions_path, checkIfExists: false), stratification_resource_paths,
                    stratification_tsv_name, benchmark_name, benchmark_version
                ))
                stratified_comparison_inputs = stratified_comparison_source
                    .combine(stratified_comparison_context)
                    .map { meta, membership, discordant, genotype_details, benchmark_regions, stratification_resources, configured_stratification_tsv_name, configured_benchmark_name, configured_benchmark_version ->
                        tuple(meta, membership, discordant, genotype_details, benchmark_regions,
                              stratification_resources, configured_stratification_tsv_name,
                              configured_benchmark_name, configured_benchmark_version)
                    }
                STRATIFIED_GERMLINE_COMPARISON(stratified_comparison_inputs)
                report_artifacts = report_artifacts.mix(STRATIFIED_GERMLINE_COMPARISON.out.results.flatMap { _meta, summary, jaccard_summary, discordant, provenance -> [summary, jaccard_summary, discordant, provenance] })
                    .mix(STRATIFIED_GERMLINE_COMPARISON.out.genotypes)
            } else {
                stratified_comparison_source = comparison_base.map { _comparison_id, meta, membership, discordant ->
                    tuple(meta, membership, discordant)
                }
                stratified_comparison_context = channel.value(tuple(
                    file(benchmark_regions_path, checkIfExists: false), stratification_resource_paths,
                    stratification_tsv_name, benchmark_name, benchmark_version
                ))
                stratified_comparison_inputs = stratified_comparison_source
                    .combine(stratified_comparison_context)
                    .map { meta, membership, discordant, benchmark_regions, stratification_resources, configured_stratification_tsv_name, configured_benchmark_name, configured_benchmark_version ->
                        tuple(meta, membership, discordant, benchmark_regions,
                              stratification_resources, configured_stratification_tsv_name,
                              configured_benchmark_name, configured_benchmark_version)
                    }
                STRATIFIED_SOMATIC_COMPARISON(stratified_comparison_inputs)
                report_artifacts = report_artifacts.mix(STRATIFIED_SOMATIC_COMPARISON.out.results.flatMap { _meta, summary, jaccard_summary, discordant, provenance -> [summary, jaccard_summary, discordant, provenance] })
            }
        }
    }

    if (run_benchmark) {
        benchmark_context = reference
            .combine(validation.benchmark_config_manifest)
            .map { reference_fa, reference_fai, reference_dict, benchmark_config_manifest ->
                tuple(
                    reference_fa,
                    reference_fai,
                    reference_dict,
                    benchmark_config_manifest,
                    file(benchmark_truth_path, checkIfExists: false),
                    file(benchmark_truth_index_path, checkIfExists: false),
                    file(benchmark_regions_path, checkIfExists: false),
                    benchmark_name,
                    benchmark_version,
                    benchmark_mode,
                    benchmark_assembly,
                    stratification_resource_paths,
                    stratification_tsv_name
                )
            }
        benchmark_inputs = processed_calls.combine(benchmark_context)
        if (benchmark_mode == 'germline_happy') {
            happy_inputs = benchmark_inputs.map { meta, caller_name, query_vcf, query_vcf_index, _bcftools_stats, _variant_summary, reference_fa, reference_fai, reference_dict, benchmark_config_manifest, truth_vcf, truth_vcf_index, benchmark_regions, configured_benchmark_name, configured_benchmark_version, configured_benchmark_mode, configured_benchmark_assembly, stratification_resources, configured_stratification_tsv_name ->
                tuple(
                    meta, caller_name, query_vcf, query_vcf_index,
                    truth_vcf, truth_vcf_index, benchmark_regions,
                    reference_fa, reference_fai, reference_dict, benchmark_config_manifest,
                    configured_benchmark_name, configured_benchmark_version,
                    configured_benchmark_mode, configured_benchmark_assembly,
                    stratification_resources, configured_stratification_tsv_name
                )
            }
            HAPPY_BENCHMARK(happy_inputs)
            report_artifacts = report_artifacts.mix(HAPPY_BENCHMARK.out.results.flatMap { _meta, _caller_name, summary, provenance, _native -> [summary, provenance] })
                .mix(HAPPY_BENCHMARK.out.stratified_summary, HAPPY_BENCHMARK.out.stratified_provenance)
        } else {
            exact_inputs = benchmark_inputs.map { meta, caller_name, query_vcf, query_vcf_index, _bcftools_stats, _variant_summary, reference_fa, reference_fai, reference_dict, benchmark_config_manifest, truth_vcf, truth_vcf_index, benchmark_regions, configured_benchmark_name, configured_benchmark_version, configured_benchmark_mode, configured_benchmark_assembly, _stratification_resources, _configured_stratification_tsv_name ->
                tuple(
                    meta, caller_name, query_vcf, query_vcf_index,
                    truth_vcf, truth_vcf_index, benchmark_regions,
                    reference_fa, reference_fai, reference_dict, benchmark_config_manifest,
                    params.analysis_type.toString(), configured_benchmark_name,
                    configured_benchmark_version, configured_benchmark_mode,
                    configured_benchmark_assembly
                )
            }
            EXACT_ALLELE_BENCHMARK(exact_inputs)
            report_artifacts = report_artifacts.mix(EXACT_ALLELE_BENCHMARK.out.results.flatMap { _meta, _caller_name, summary, provenance, _native -> [summary, provenance] })
            if (stratification_tsv_name) {
                stratified_exact_context = channel.value(tuple(
                    file(benchmark_regions_path, checkIfExists: false), stratification_resource_paths,
                    stratification_tsv_name, params.analysis_type.toString(), benchmark_name,
                    benchmark_version, benchmark_assembly
                ))
                stratified_exact_inputs = EXACT_ALLELE_BENCHMARK.out.alleles
                    .combine(stratified_exact_context)
                    .map { meta, caller_name, truth_alleles, query_alleles, benchmark_regions, stratification_resources, configured_stratification_tsv_name, configured_analysis_type, configured_benchmark_name, configured_benchmark_version, configured_benchmark_assembly ->
                        tuple(meta, caller_name, truth_alleles, query_alleles, benchmark_regions,
                              stratification_resources, configured_stratification_tsv_name,
                              configured_analysis_type, configured_benchmark_name,
                              configured_benchmark_version, configured_benchmark_assembly)
                    }
                STRATIFIED_EXACT_BENCHMARK(stratified_exact_inputs)
                report_artifacts = report_artifacts.mix(STRATIFIED_EXACT_BENCHMARK.out.results.flatMap { _meta, _caller_name, summary, _errors, provenance -> [summary, provenance] })
            }
        }
    }

    if (run_vep) {
        if (vep_mode == 'cache') {
            vep_annotation_resources = channel.value(tuple('cache', [file(vep_cache_path, checkIfExists: false)]))
        } else {
            vep_annotation_resources = channel.value(tuple('gtf', [
                file(vep_gtf_path, checkIfExists: false),
                file("${vep_gtf_path}.tbi", checkIfExists: false)
            ]))
        }

        vep_inputs = processed_calls
            .combine(reference)
            .combine(vep_annotation_resources)
            .map { meta, caller_name, normalized_vcf, normalized_vcf_index, _bcftools_stats, _variant_summary, reference_fa, reference_fai, _reference_dict, annotation_mode, annotation_resources ->
                tuple(
                    meta,
                    caller_name,
                    normalized_vcf,
                    normalized_vcf_index,
                    reference_fa,
                    reference_fai,
                    annotation_mode,
                    annotation_resources
                )
            }
        annotated_calls = VEP(vep_inputs).annotated
    }

    if (run_vcf2circos) {
        if (vcf2circos_mode == 'regional') {
            if (run_vep) {
                visualization_inputs = annotated_calls.map { meta, caller_name, annotated_vcf, annotated_vcf_index ->
                    tuple(meta, caller_name, annotated_vcf, annotated_vcf_index, 'vep', params.analysis_type.toString(), calling_interval)
                }
            } else {
                visualization_inputs = processed_calls.map { meta, caller_name, normalized_vcf, normalized_vcf_index, _bcftools_stats, _variant_summary ->
                    tuple(meta, caller_name, normalized_vcf, normalized_vcf_index, 'normalized', params.analysis_type.toString(), calling_interval)
                }
            }

            visualization_prepared = VCF2CIRCOS_PREP(visualization_inputs).prepared
            vcf2circos_resources = channel.value(file(vcf2circos_static_path, checkIfExists: false))
            vcf2circos_inputs = visualization_prepared
                .combine(vcf2circos_resources)
                .map { meta, caller_name, visualization_vcf, visualization_vcf_index, visualization_summary, input_kind, configured_interval, static_resources ->
                    tuple(meta, caller_name, visualization_vcf, visualization_vcf_index, visualization_summary, input_kind, configured_interval, static_resources)
                }
            vcf2circos_results = VCF2CIRCOS(vcf2circos_inputs).visualizations
            report_artifacts = report_artifacts.mix(vcf2circos_results.flatMap { _meta, _caller_name, html_report, summary -> [html_report, summary] })
        } else {
            genome_reference_fai = channel.value(file(reference_fai_path, checkIfExists: false))
            genome_empty_tracks = VCF2CIRCOS_EMPTY_TRACKS()
            if (run_coverage_qc) {
                genome_coverage_labels = mosdepth_results.coverage
                    .map { meta, _summary, _global_dist, _region_dist, _regions, _regions_index, _threshold_bed, _coverage_regions -> "${meta.role}:${meta.sample_id}" }
                    .collect()
                genome_coverage_files = mosdepth_results.coverage
                    .map { _meta, _summary, _global_dist, _region_dist, regions, _regions_index, _threshold_bed, _coverage_regions -> regions }
                    .collect()
            } else {
                genome_coverage_labels = channel.value([])
                genome_coverage_files = genome_empty_tracks.coverage.map { empty_track -> [empty_track] }
            }
            genome_sv_files = run_sv ? manta_calls.map { values -> values[2] }.collect() : genome_empty_tracks.sv.map { empty_track -> [empty_track] }
            if (run_cnv && params.analysis_type.toString() == 'germline') {
                genome_cnv_files = germline_cnv_final.map { _meta, _interval_vcf, segment_vcf, _ratios -> segment_vcf }.collect()
                genome_cnv_format = 'germline_vcf'
            } else if (run_cnv) {
                genome_cnv_files = somatic_cnv_final.map { _meta, _cr, _model, _hets, called -> called }.collect()
                genome_cnv_format = 'somatic_seg'
            } else {
                genome_cnv_files = genome_empty_tracks.cnv.map { empty_track -> [empty_track] }
                genome_cnv_format = 'NA'
            }
            genome_inputs = processed_calls
                .map { meta, caller_name, normalized_vcf, _normalized_vcf_index, _bcftools_stats, _variant_summary -> tuple(meta, caller_name, normalized_vcf) }
                .combine(genome_reference_fai)
                .combine(channel.value(params.analysis_type.toString()))
                .combine(channel.value(vcf2circos_bin_size))
                .combine(genome_coverage_labels)
                .combine(genome_coverage_files)
                .combine(genome_sv_files)
                .combine(genome_cnv_files)
                .combine(channel.value(genome_cnv_format))
            vcf2circos_results = VCF2CIRCOS_GENOME(genome_inputs).visualizations
            report_artifacts = report_artifacts.mix(vcf2circos_results.flatMap { _meta, _caller_name, html_report, summary, provenance -> [html_report, summary, provenance] })
        }
    }

    if (run_multiqc) {
        multiqc_artifacts = fastqc_reports.map { _meta, _fastqc_html, fastqc_zip -> fastqc_zip }
            .mix(fastp_results.map { _meta, _read1, _read2, _fastp_html, fastp_json -> fastp_json })
            .mix(markdup_results.marked.map { _meta, _bam, _bai, markdup_stats -> markdup_stats })
            .mix(markdup_results.qc.flatMap { _meta, flagstat, samtools_stats -> [flagstat, samtools_stats] })
            .mix(processed_calls.map { _meta, _caller_name, _vcf, _tbi, bcftools_stats, _variant_summary -> bcftools_stats })

        if (run_bqsr) {
            multiqc_artifacts = multiqc_artifacts.mix(recalibrated.map { _meta, _bam, recal_table -> recal_table })
        }
        if (run_coverage_qc) {
            multiqc_artifacts = multiqc_artifacts.mix(mosdepth_results.coverage.flatMap { _meta, summary, global_dist, region_dist, _regions, _regions_index, _threshold_bed, _coverage_regions ->
                [summary, global_dist, region_dist]
            })
        }

        multiqc_inputs = multiqc_artifacts.collect().map { files ->
            tuple(params.analysis_type.toString(), files)
        }
        multiqc_results = MULTIQC(multiqc_inputs).report
        report_artifacts = report_artifacts.mix(multiqc_results.map { _analysis_type, html_report, _data -> html_report })
    }

    if (run_analysis_report) {
        report_inputs = report_artifacts.collect()
        report_process_input = report_inputs.map { files -> tuple(
            files,
            samplesheet,
            report_identifier,
            params.analysis_type.toString(),
            selected_caller,
            assembly,
            resource_mode,
            workflow.manifest.version?.toString() ?: 'NA',
            nextflow.version.toString(),
            run_bqsr,
            run_benchmark,
            run_vep,
            run_sv,
            run_cnv,
            run_coverage_qc,
            run_multiqc,
            run_vcf2circos,
            vcf2circos_mode
        ) }
        ANALYSIS_REPORT(report_process_input)
    }

    workflow.onComplete = {
        def completion_performance_flag = params.run_performance_report?.toString()?.trim()?.toLowerCase() ?: 'false'
        if (completion_performance_flag == 'true') {
            def performance_dir = file("${params.outdir}/performance")
            performance_dir.mkdirs()
            def trace_path = file("${performance_dir}/nextflow/trace.txt")
            def elapsed_seconds = workflow.duration.toMillis() / 1000.0
            def benchmark_enabled = params.run_benchmark?.toString()?.toBoolean() ?: false
            def stratification_enabled = benchmark_enabled && !!(params.benchmark_stratifications?.toString()?.trim() ?: '')
            def command = [
                'python3', "${projectDir}/bin/parse_nextflow_trace.py",
                '--trace', trace_path.toString(),
                '--output-dir', performance_dir.toString(),
                '--analysis', params.analysis_type?.toString() ?: '',
                '--caller', params.caller?.toString() ?: '',
                '--resource-mode', params.resource_mode?.toString() ?: 'fixture',
                '--reference-assembly', params.assembly?.toString() ?: 'NA',
                '--bqsr', (params.run_bqsr?.toString()?.toBoolean() ?: false).toString(),
                '--benchmark', benchmark_enabled.toString(),
                '--benchmark-mode', benchmark_enabled ? benchmark_mode.toString() : 'NA',
                '--benchmark-name', benchmark_enabled ? (benchmark_name ?: 'NA') : 'NA',
                '--benchmark-version', benchmark_enabled ? (benchmark_version ?: 'NA') : 'NA',
                '--comparison', (params.caller?.toString() == 'compare').toString(),
                '--stratification', stratification_enabled.toString(),
                '--sv', (params.run_sv?.toString()?.toBoolean() ?: false).toString(),
                '--cnv', (params.run_cnv?.toString()?.toBoolean() ?: false).toString(),
                '--coverage-qc', (params.run_coverage_qc?.toString()?.toBoolean() ?: false).toString(),
                '--ploidy-configured', (!!(params.ploidy_config?.toString()?.trim() ?: '')).toString(),
                '--vep', (params.run_vep?.toString()?.toBoolean() ?: false).toString(),
                '--visualization', (params.run_vcf2circos?.toString()?.toBoolean() ?: false).toString(),
                '--multiqc', (params.run_multiqc?.toString()?.toBoolean() ?: false).toString(),
                '--executor', params.performance_executor?.toString() ?: 'NA',
                '--container-engine', params.performance_container_engine?.toString() ?: 'NA',
                '--nextflow-version', nextflow.version.toString(),
                '--pipeline-version', workflow.manifest.version?.toString() ?: 'NA',
                '--workflow-start', workflow.start.toString(),
                '--workflow-completion', workflow.complete.toString(),
                '--workflow-elapsed-seconds', elapsed_seconds.toString(),
                '--workflow-success', workflow.success.toString()
            ]
            def process = new ProcessBuilder(command).inheritIO().start()
            def exit_code = process.waitFor()
            if (exit_code != 0) {
                throw new IllegalStateException("Performance report generation failed with exit status ${exit_code}; inspect the native trace at ${trace_path}")
            }
        }
        def completion_report_flag = params.run_analysis_report?.toString()?.trim()?.toLowerCase() ?: 'false'
        if (completion_report_flag == 'true') {
            def report_output = file("${params.outdir}/report/${report_identifier}.${params.analysis_type}.wgs_report.html")
            def report_command = [
                'python3', "${projectDir}/bin/generate_analysis_report.py",
                '--input-root', file(params.outdir).toString(),
                '--output', report_output.toString(),
                '--samplesheet', input_sheet.toString(),
                '--report-id', report_identifier.toString(),
                '--analysis-type', params.analysis_type.toString(),
                '--caller', selected_caller.toString(),
                '--assembly', assembly.toString(),
                '--resource-mode', resource_mode.toString(),
                '--pipeline-version', workflow.manifest.version?.toString() ?: 'NA',
                '--nextflow-version', nextflow.version.toString(),
                '--benchmark-name', run_benchmark ? benchmark_name.toString() : '',
                '--benchmark-version', run_benchmark ? benchmark_version.toString() : ''
            ]
            if (run_bqsr) report_command << '--run-bqsr'
            if (run_benchmark) report_command << '--run-benchmark'
            if (run_vep) report_command << '--run-vep'
            if (run_sv) report_command << '--run-sv'
            if (run_cnv) report_command << '--run-cnv'
            if (run_coverage_qc) report_command << '--run-coverage-qc'
            if (run_multiqc) report_command.addAll(['--multiqc-link', "../multiqc/${params.analysis_type}.multiqc.html"])
            if (run_vcf2circos) {
                def display_caller = selected_caller == 'compare' ? (params.analysis_type.toString() == 'germline' ? 'haplotypecaller' : 'mutect2') : selected_caller
                def visualization_suffix = vcf2circos_mode == 'genome' ? '.genome.vcf2circos.html' : '.vcf2circos.html'
                report_command.addAll(['--vcf2circos-link', "../vcf2circos/${report_identifier}.${display_caller}${visualization_suffix}"])
            }
            def report_process = new ProcessBuilder(report_command).inheritIO().start()
            def report_exit_code = report_process.waitFor()
            if (report_exit_code != 0) {
                throw new IllegalStateException("Integrated report finalization failed with exit status ${report_exit_code}")
            }
        }
    }
}
