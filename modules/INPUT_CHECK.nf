process INPUT_CHECK {
    label 'process_low'
    tag "input validation"

    publishDir "${params.outdir ?: 'results'}/provenance", mode: 'copy', overwrite: true, pattern: '*.tsv'

    input:
    val input_path
    val reference_path
    val analysis_type
    val caller
    val outdir
    val run_bqsr
    val bqsr_interval
    val calling_interval
    val resource_mode
    val assembly
    val project_dir
    val reference_fai_path
    val reference_dict_path
    val bwa_index_prefix
    val known_sites_joined
    val resource_manifest_path
    val run_vep
    val vep_mode
    val vep_cache_path
    val vep_gtf_path
    val run_vcf2circos
    val vcf2circos_static_path
    val gatk_germline_mode
    val gatk_somatic_profile
    val mutect2_germline_resource_path
    val mutect2_pon_path
    val mutect2_contamination_sites_path
    val run_benchmark
    val benchmark_mode
    val benchmark_truth_path
    val benchmark_truth_index_path
    val benchmark_regions_path
    val benchmark_assembly
    val benchmark_name
    val benchmark_version
    val benchmark_sex
    val benchmark_stratifications_path
    val benchmark_engine_container
    val run_sv
    val run_cnv
    val run_coverage_qc
    val coverage_thresholds
    val ploidy_config_path
    val cnv_intervals_path
    val cnv_ploidy_model_path
    val cnv_germline_model_path
    val cnv_read_count_pon_path
    val cnv_common_sites_path

    output:
    path '.input_check.ok', emit: passed
    path 'resource_manifest.tsv', emit: resource_manifest
    path 'gatk_path_manifest.tsv', emit: gatk_path_manifest
    path 'benchmark_config_manifest.tsv', emit: benchmark_config_manifest
    path 'wgs_modality_provenance.tsv', emit: wgs_modality_provenance
    path 'ploidy_summary.tsv', emit: ploidy_summary

    script:
    """
    set -euo pipefail
    CHECK_INPUT='${input_path}'
    CHECK_REFERENCE='${reference_path}'
    CHECK_ANALYSIS='${analysis_type}'
    CHECK_CALLER='${caller}'
    CHECK_OUTDIR='${outdir}'
    CHECK_CALLING_INTERVAL='${calling_interval}'

    fail() {
        echo "INPUT_CHECK ERROR: \$1" >&2
        exit 1
    }

    case "\$CHECK_ANALYSIS" in
        germline|somatic) ;;
        *) fail "unsupported analysis_type '\$CHECK_ANALYSIS'. Expected germline or somatic" ;;
    esac
    case "\$CHECK_CALLER" in
        deepvariant|haplotypecaller|deepsomatic|mutect2|strelka2|compare) ;;
        *) fail "unsupported caller '\$CHECK_CALLER'" ;;
    esac
    case "\$CHECK_ANALYSIS:\$CHECK_CALLER" in
        germline:deepvariant|germline:haplotypecaller|germline:compare|somatic:deepsomatic|somatic:mutect2|somatic:strelka2|somatic:compare) ;;
        germline:*) fail "caller '\$CHECK_CALLER' is incompatible with analysis_type germline. Expected deepvariant, haplotypecaller or compare" ;;
        somatic:*) fail "caller '\$CHECK_CALLER' is incompatible with analysis_type somatic. Expected deepsomatic, mutect2, strelka2 or compare" ;;
    esac
    case '${gatk_germline_mode}' in
        vcf|gvcf) ;;
        *) fail "unsupported gatk_germline_mode '${gatk_germline_mode}'. Expected vcf or gvcf" ;;
    esac
    case '${gatk_somatic_profile}' in
        basic|production) ;;
        *) fail "unsupported gatk_somatic_profile '${gatk_somatic_profile}'. Expected basic or production" ;;
    esac
    if [[ '${gatk_germline_mode}' == gvcf && "\$CHECK_ANALYSIS:\$CHECK_CALLER" != germline:haplotypecaller && "\$CHECK_ANALYSIS:\$CHECK_CALLER" != germline:compare ]]; then
        fail "gatk_germline_mode gvcf requires analysis_type germline and caller haplotypecaller or compare"
    fi
    if [[ '${gatk_somatic_profile}' == production && "\$CHECK_ANALYSIS:\$CHECK_CALLER" != somatic:mutect2 && "\$CHECK_ANALYSIS:\$CHECK_CALLER" != somatic:compare ]]; then
        fail "gatk_somatic_profile production requires analysis_type somatic and caller mutect2 or compare"
    fi

    [[ -n "\$CHECK_INPUT" && -f "\$CHECK_INPUT" ]] || fail "missing samplesheet: '\$CHECK_INPUT'"
    [[ -n "\$CHECK_REFERENCE" && -f "\$CHECK_REFERENCE" ]] || fail "missing reference FASTA: '\$CHECK_REFERENCE'"
    [[ -n "\$CHECK_OUTDIR" ]] || fail "missing required parameter: outdir"
    if [[ -n "\$CHECK_CALLING_INTERVAL" ]]; then
        if [[ "\$CHECK_CALLING_INTERVAL" =~ ^([A-Za-z0-9_.-]+):([0-9]+)-([0-9]+)\$ ]]; then
            interval_start="\${BASH_REMATCH[2]}"
            interval_end="\${BASH_REMATCH[3]}"
            (( interval_start >= 1 && interval_end >= interval_start )) || \
                fail "invalid calling_interval '\$CHECK_CALLING_INTERVAL': require 1-based start <= end"
        else
            fail "invalid calling_interval '\$CHECK_CALLING_INTERVAL': expected contig:start-end"
        fi
    fi
    if [[ "${run_bqsr}" == true ]]; then
        [[ -n "${bqsr_interval}" ]] || fail "missing BQSR interval"
        if [[ '${resource_mode}' == fixture ]]; then
            if [[ '${bqsr_interval}' == auto ]]; then
                fixture_bqsr_interval="\$CHECK_CALLING_INTERVAL"
            else
                fixture_bqsr_interval='${bqsr_interval}'
            fi
            case "\$fixture_bqsr_interval" in
                chr20:12446968-12455543|chr20:21900000-21980636) ;;
                *) fail "fixture BQSR is restricted to the two known-sites fixture intervals; set calling_interval/bqsr_interval to a documented fixture interval or use production resources" ;;
            esac
        fi
    fi

    header=\$(head -n 1 "\$CHECK_INPUT")
    [[ "\$header" == \$'sample_id\\tpatient_id\\trole\\tfastq_r1\\tfastq_r2' ]] || \
        fail "malformed samplesheet: expected header sample_id\\tpatient_id\\trole\\tfastq_r1\\tfastq_r2"

    sheet_dir=\$(cd "\$(dirname "\$CHECK_INPUT")" && pwd)
    declare -A samples tumor_by_patient normal_by_patient patient_ids
    row_count=0
    tumor_count=0
    normal_count=0

    while IFS=\$'\\t' read -r sample_id patient_id role fastq_r1 fastq_r2 extra; do
        [[ -z "\$extra" ]] || fail "malformed samplesheet: row has more than 5 tab-delimited columns"
        ((row_count+=1))
        [[ -n "\$sample_id" && -n "\$patient_id" && -n "\$role" && -n "\$fastq_r1" && -n "\$fastq_r2" ]] || \
            fail "malformed samplesheet: row \$row_count has an empty field"
        [[ -z "\${samples[\$sample_id]+x}" ]] || fail "duplicate sample identifier '\$sample_id'"
        samples[\$sample_id]=1
        patient_ids[\$patient_id]=1

        case "\$CHECK_ANALYSIS:\$role" in
            germline:sample) ;;
            germline:*) fail "malformed samplesheet: germline role must be sample (row \$row_count)" ;;
            somatic:tumor) ((tumor_count+=1)); tumor_by_patient["\$patient_id"]="\${tumor_by_patient["\$patient_id"]-0}"; ((tumor_by_patient["\$patient_id"]+=1)) ;;
            somatic:normal) ((normal_count+=1)); normal_by_patient["\$patient_id"]="\${normal_by_patient["\$patient_id"]-0}"; ((normal_by_patient["\$patient_id"]+=1)) ;;
            somatic:*) fail "malformed samplesheet: somatic role must be tumor or normal (row \$row_count)" ;;
        esac

        for fastq in "\$fastq_r1" "\$fastq_r2"; do
            [[ "\$fastq" != *\$'\\t'* && "\$fastq" != *\$'\\n'* ]] || fail "invalid paired-end FASTQ definition for sample '\$sample_id'"
            [[ "\$fastq" == *.fastq || "\$fastq" == *.fq || "\$fastq" == *.fastq.gz || "\$fastq" == *.fq.gz ]] || \
                fail "invalid paired-end FASTQ definition for sample '\$sample_id': expected .fastq/.fq, optionally .gz"
            if [[ "\$fastq" = /* ]]; then fq_path="\$fastq"; else fq_path="\$sheet_dir/\$fastq"; fi
            [[ -f "\$fq_path" ]] || fail "missing FASTQ file for sample '\$sample_id': '\$fastq'"
        done
    done < <(tail -n +2 "\$CHECK_INPUT")

    ((row_count > 0)) || fail "malformed samplesheet: no sample rows found"
    if [[ "\$CHECK_ANALYSIS" == somatic ]]; then
        ((tumor_count > 0 && normal_count > 0)) || fail "somatic samplesheet must contain tumor and normal rows"
        for patient in "\${!patient_ids[@]}"; do
            [[ "\${tumor_by_patient[\$patient]-0}" == 1 && "\${normal_by_patient[\$patient]-0}" == 1 ]] || \
                fail "somatic samplesheet must contain exactly one tumor and one normal for patient '\$patient'"
        done
    fi

    validate_resource_config.sh \
        '${resource_mode}' \
        '${assembly}' \
        '${project_dir}' \
        '${reference_path}' \
        '${reference_fai_path}' \
        '${reference_dict_path}' \
        '${bwa_index_prefix}' \
        '${known_sites_joined}' \
        '${resource_manifest_path}' \
        '${run_bqsr}' \
        '${run_vep}' \
        '${vep_mode}' \
        '${vep_cache_path}' \
        '${vep_gtf_path}' \
        '${run_vcf2circos}' \
        '${vcf2circos_static_path}' \
        '${analysis_type}' \
        '${caller}' \
        '${gatk_germline_mode}' \
        '${gatk_somatic_profile}' \
        '${mutect2_germline_resource_path}' \
        '${mutect2_pon_path}' \
        '${mutect2_contamination_sites_path}'

    validate_benchmark_config.sh \
        '${run_benchmark}' \
        '${benchmark_mode}' \
        '${resource_mode}' \
        '${analysis_type}' \
        '${caller}' \
        '${benchmark_name}' \
        '${benchmark_version}' \
        '${benchmark_sex}' \
        '${benchmark_assembly}' \
        '${assembly}' \
        '${project_dir}' \
        '${reference_path}' \
        '${reference_fai_path}' \
        '${benchmark_truth_path}' \
        '${benchmark_truth_index_path}' \
        '${benchmark_regions_path}' \
        '${benchmark_stratifications_path}' \
        '${resource_manifest_path}' \
        '${benchmark_engine_container}'

    validate_wgs_modalities.py \
        --samplesheet '${input_path}' \
        --analysis '${analysis_type}' \
        --assembly '${assembly}' \
        --resource-mode '${resource_mode}' \
        --resource-manifest resource_manifest.tsv \
        --project-dir '${project_dir}' \
        --reference-fai '${reference_fai_path}' \
        --run-sv '${run_sv}' \
        --run-cnv '${run_cnv}' \
        --run-coverage '${run_coverage_qc}' \
        --coverage-thresholds '${coverage_thresholds}' \
        --ploidy-config '${ploidy_config_path}' \
        --cnv-intervals '${cnv_intervals_path}' \
        --cnv-ploidy-model '${cnv_ploidy_model_path}' \
        --cnv-germline-model '${cnv_germline_model_path}' \
        --cnv-read-count-pon '${cnv_read_count_pon_path}' \
        --cnv-common-sites '${cnv_common_sites_path}' \
        --provenance-out wgs_modality_provenance.tsv \
        --ploidy-out ploidy_summary.tsv

    orientation_model=false
    contamination_estimation=false
    filtering_path=''
    if [[ '${analysis_type}:${gatk_somatic_profile}' == somatic:production && ( '${caller}' == mutect2 || '${caller}' == compare ) ]]; then
        orientation_model=true
        contamination_estimation=true
        filtering_path='Mutect2+LearnReadOrientationModel+GetPileupSummaries+CalculateContamination+FilterMutectCalls'
    elif [[ '${analysis_type}:${caller}' == somatic:mutect2 || '${analysis_type}:${caller}' == somatic:compare ]]; then
        filtering_path='Mutect2+FilterMutectCalls'
    elif [[ '${analysis_type}:${gatk_germline_mode}' == germline:gvcf && ( '${caller}' == haplotypecaller || '${caller}' == compare ) ]]; then
        filtering_path='HaplotypeCaller(GVCF)+GenotypeGVCFs;raw_unfiltered_callset'
    elif [[ '${analysis_type}:${caller}' == germline:haplotypecaller || '${analysis_type}:${caller}' == germline:compare ]]; then
        filtering_path='HaplotypeCaller(VCF);raw_unfiltered_callset'
    fi
    printf 'analysis_type\tcaller\tgatk_version\tgatk_germline_mode\tgatk_somatic_profile\torientation_model\tcontamination_estimation\tfiltering_path\n' > gatk_path_manifest.tsv
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        '${analysis_type}' '${caller}' '4.6.2.0' '${gatk_germline_mode}' '${gatk_somatic_profile}' \
        "\$orientation_model" "\$contamination_estimation" "\$filtering_path" >> gatk_path_manifest.tsv

    mkdir -p "\$CHECK_OUTDIR"
    touch .input_check.ok
    echo "INPUT_CHECK: valid \$CHECK_ANALYSIS samplesheet (\$row_count sample(s))"
    """
}
