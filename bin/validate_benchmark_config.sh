#!/usr/bin/env bash
set -euo pipefail

fail() {
    echo "BENCHMARK_CHECK ERROR: $1" >&2
    exit 1
}

run_benchmark=$1
benchmark_mode=$2
resource_mode=$3
analysis_type=$4
caller=$5
benchmark_name=$6
benchmark_version=$7
benchmark_sex=$8
benchmark_assembly=$9
pipeline_assembly=${10}
project_dir=${11}
reference=${12}
reference_fai=${13}
truth_vcf=${14}
truth_vcf_index=${15}
benchmark_regions=${16}
benchmark_stratifications=${17}
resource_manifest=${18}
engine_container=${19}

config_header=$'benchmark_name\tbenchmark_version\tbenchmark_mode\tbenchmark_assembly\tbenchmark_sex_requested\tbenchmark_sex_resolved\tsex_resolution_source\thap_py_gender_argument\tsex_aware_benchmarking\tengine\tengine_version\tengine_container\ttruth_vcf\ttruth_vcf_index\ttruth_version\ttruth_source\ttruth_checksum\ttruth_index_checksum\tbenchmark_regions\tregions_version\tregions_source\tregions_checksum\tbenchmark_stratifications\tstratifications_version\tstratifications_source\tstratifications_checksum\tstratification_overlap_semantics\teffective_base_semantics'
printf '%s\n' "$config_header" > benchmark_config_manifest.tsv

[[ "$run_benchmark" == true || "$run_benchmark" == false ]] || \
    fail "run_benchmark must be true or false"
case "$benchmark_sex" in
    auto|male|female|none) ;;
    *) fail "unsupported benchmark_sex '$benchmark_sex'. Expected auto, male, female, or none" ;;
esac
if [[ "$run_benchmark" == false ]]; then
    [[ -z "$benchmark_stratifications" ]] || fail "benchmark_stratifications requires run_benchmark true"
    [[ "$benchmark_sex" == auto || "$benchmark_sex" == none ]] || \
        fail "benchmark_sex $benchmark_sex applies only to enabled germline hap.py benchmarking"
    exit 0
fi

case "$benchmark_mode" in
    regression|germline_happy|somatic) ;;
    *) fail "unsupported benchmark_mode '$benchmark_mode'. Expected regression, germline_happy, or somatic" ;;
esac

if [[ "$benchmark_mode" == regression && "$resource_mode" != fixture ]]; then
    fail "benchmark_mode regression is restricted to resource_mode fixture"
fi
if [[ "$benchmark_mode" == germline_happy && "$analysis_type" != germline ]]; then
    fail "benchmark_mode germline_happy requires analysis_type germline"
fi
if [[ "$benchmark_mode" == somatic && "$analysis_type" != somatic ]]; then
    fail "benchmark_mode somatic requires analysis_type somatic"
fi
if [[ "$benchmark_mode" != germline_happy && "$benchmark_sex" != auto && "$benchmark_sex" != none ]]; then
    fail "benchmark_sex $benchmark_sex applies only to germline_happy benchmarking"
fi
[[ -n "$benchmark_name" ]] || fail "benchmark_name is required when benchmarking is enabled"
[[ -n "$benchmark_assembly" ]] || fail "benchmark_assembly is required when benchmarking is enabled"
[[ "$benchmark_assembly" == "$pipeline_assembly" ]] || \
    fail "benchmark assembly '$benchmark_assembly' conflicts with configured reference assembly '$pipeline_assembly'"

for item in "$truth_vcf" "$truth_vcf_index" "$benchmark_regions"; do
    [[ -r "$item" && -s "$item" ]] || fail "missing or unreadable benchmark resource: '$item'"
done
[[ "$truth_vcf" == *.vcf.gz ]] || fail "benchmark truth must be a bgzipped .vcf.gz file: '$truth_vcf'"
[[ "$truth_vcf_index" == "${truth_vcf}.tbi" ]] || \
    fail "benchmark truth index must be the tabix companion '${truth_vcf}.tbi'"
if [[ "$resource_mode" == production ]]; then
    for item in "$truth_vcf" "$truth_vcf_index" "$benchmark_regions" "$benchmark_stratifications"; do
        [[ -z "$item" || "$item" != "$project_dir/test_data/"* ]] || \
            fail "production benchmarking cannot use repository fixture resource '$item'"
    done
fi

first_header=$(awk 'NR == 1 { print; exit }' < <(gzip -dc "$truth_vcf" 2>/dev/null))
[[ "$first_header" == '##fileformat=VCF'* ]] || \
    fail "benchmark truth is not a readable bgzipped VCF: '$truth_vcf'"

if [[ -n "$reference_fai" && -r "$reference_fai" ]]; then
    cut -f1,2 "$reference_fai" > .benchmark_reference_contigs.tsv
else
    awk '
        /^>/ {
            if (name != "") print name "\t" sequence_length
            name=$0
            sub(/^>/, "", name)
            sub(/[[:space:]].*$/, "", name)
            sequence_length=0
            next
        }
        { gsub(/[[:space:]]/, ""); sequence_length += length($0) }
        END { if (name != "") print name "\t" sequence_length }
    ' < <(gzip -dc "$reference") > .benchmark_reference_contigs.tsv
fi
[[ -s .benchmark_reference_contigs.tsv ]] || fail "could not derive reference contigs for benchmark validation"

awk '
    /^##contig=<ID=/ {
        value=$0
        sub(/^##contig=<ID=/, "", value)
        split(value, fields, /[,>]/)
        name=fields[1]
        sequence_length=""
        if (match($0, /length=[0-9]+/)) {
            sequence_length=substr($0, RSTART+7, RLENGTH-7)
        }
        print name "\t" sequence_length
    }
    !/^#/ { if (!seen[$1]++) print $1 "\t"; exit }
' < <(gzip -dc "$truth_vcf" 2>/dev/null) > .benchmark_truth_contigs.tsv
[[ -s .benchmark_truth_contigs.tsv ]] || fail "benchmark truth declares no contigs and contains no records"

awk -F '\t' '
    NR==FNR { ref[$1]=$2; next }
    $1 in ref {
        shared=1
        if ($2 != "" && ref[$1] != $2) mismatch=1
    }
    END { exit(shared && !mismatch ? 0 : 1) }
' .benchmark_reference_contigs.tsv .benchmark_truth_contigs.tsv || \
    fail "benchmark truth/reference contigs or declared lengths are incompatible"

read_bed() {
    local bed=$1
    if [[ "$bed" == *.gz ]]; then gzip -dc "$bed"; else command cat "$bed"; fi
}

validate_bed() {
    local label=$1 bed=$2 output=$3
    [[ -r "$bed" && -s "$bed" ]] || fail "missing or unreadable $label BED: '$bed'"
    read_bed "$bed" | awk 'BEGIN { OFS="\t" }
        /^#/ || /^track([[:space:]]|$)/ || /^browser([[:space:]]|$)/ { next }
        NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 <= $2 { exit 1 }
        { print $1, $2, $3; rows++ }
        END { if (rows == 0) exit 1 }
    ' > "$output" || fail "$label BED is malformed or contains no positive-width regions: '$bed'"
    awk -F '\t' '
        NR==FNR { ref[$1]=1; next }
        !($1 in ref) { exit 1 }
        { rows++ }
        END { if (rows == 0) exit 1 }
    ' .benchmark_reference_contigs.tsv "$output" || \
        fail "$label BED contigs are incompatible with the configured reference"
    awk -F '\t' '
        NR==FNR { ref[$1]=$2; next }
        $2 < 0 || $3 > ref[$1] { exit 1 }
        { bases += $3-$2; rows++ }
        END { if (rows == 0 || bases == 0) exit 1 }
    ' .benchmark_reference_contigs.tsv "$output" || \
        fail "$label BED has no valid non-empty overlap with the configured reference"
}

validate_bed "benchmark region" "$benchmark_regions" .benchmark_regions.tsv
awk -F '\t' '
    NR==FNR { truth[$1]=1; next }
    truth[$1] { shared=1 }
    END { exit(shared ? 0 : 1) }
' .benchmark_truth_contigs.tsv .benchmark_regions.tsv || \
    fail "benchmark truth and benchmark BED have no shared contigs"

if [[ -n "$benchmark_stratifications" ]]; then
    [[ -r "$benchmark_stratifications" && -s "$benchmark_stratifications" ]] || \
        fail "missing or unreadable benchmark stratification TSV: '$benchmark_stratifications'"
    stratification_dir=$(cd "$(dirname "$benchmark_stratifications")" && pwd)
    declare -A stratification_names stratification_basenames
    stratification_rows=0
    while IFS=$'\t' read -r region_name bed_path extra; do
        [[ -z "$region_name$bed_path$extra" ]] && continue
        [[ "$region_name" == \#* ]] && continue
        [[ -n "$region_name" && -n "$bed_path" && -z "$extra" ]] || \
            fail "malformed stratification TSV; expected exactly region_name<TAB>BED_path per row"
        [[ -z "${stratification_names[$region_name]+x}" ]] || \
            fail "duplicate stratification name '$region_name'"
        stratification_names[$region_name]=1
        if [[ "$bed_path" = /* ]]; then resolved_bed=$bed_path; else resolved_bed="$stratification_dir/$bed_path"; fi
        bed_basename=${resolved_bed##*/}
        [[ -z "${stratification_basenames[$bed_basename]+x}" ]] || \
            fail "stratification BED basenames must be unique for portable staging: '$bed_basename'"
        stratification_basenames[$bed_basename]=1
        validate_bed "stratification '$region_name'" "$resolved_bed" ".benchmark_stratification_${stratification_rows}.tsv"
        awk -F '\t' '
            NR==FNR { n[$1]++; start[$1,n[$1]]=$2; stop[$1,n[$1]]=$3; next }
            {
                for (i=1; i<=n[$1]; i++) {
                    if ($2 < stop[$1,i] && $3 > start[$1,i]) { shared=1; exit }
                }
            }
            END { exit(shared ? 0 : 1) }
        ' .benchmark_regions.tsv ".benchmark_stratification_${stratification_rows}.tsv" || \
            fail "stratification '$region_name' has no non-empty overlap with benchmark regions"
        ((stratification_rows+=1))
    done < "$benchmark_stratifications"
    ((stratification_rows > 0)) || fail "benchmark stratification TSV contains no resource rows"
fi

manifest_value() {
    local target=$1 column=$2
    awk -F '\t' -v target="$target" -v target_name="${target##*/}" -v column="$column" '
        NR > 1 && ($6 == target || $6 == target_name || substr($6, length($6)-length(target_name)+1) == target_name) { print $column; exit }
    ' "$resource_manifest"
}

truth_version=$(manifest_value "$truth_vcf" 4)
truth_source=$(manifest_value "$truth_vcf" 5)
truth_checksum=$(manifest_value "$truth_vcf" 7)
truth_index_checksum=$(manifest_value "$truth_vcf_index" 7)
regions_version=$(manifest_value "$benchmark_regions" 4)
regions_source=$(manifest_value "$benchmark_regions" 5)
regions_checksum=$(manifest_value "$benchmark_regions" 7)

# Benchmark identity is established from explicit configuration plus trusted
# resource-manifest metadata. Filenames are used only as mismatch sentinels;
# they never establish sample sex.
truth_basename=${truth_vcf##*/}
regions_basename=${benchmark_regions##*/}
declared_v5=false
declared_v421=false
[[ "$benchmark_version" == v5.0q || "$benchmark_version" == 'GIAB HG002 v5.0q' ]] && declared_v5=true
[[ "$benchmark_version" == v4.2.1 || "$benchmark_version" == 'GIAB HG002 v4.2.1' ]] && declared_v421=true
truth_meta_v5=false; regions_meta_v5=false; truth_meta_v421=false; regions_meta_v421=false
[[ "$truth_version" == *'HG002 v5.0q'* ]] && truth_meta_v5=true
[[ "$regions_version" == *'HG002 v5.0q'* ]] && regions_meta_v5=true
[[ "$truth_version" == *'HG002 v4.2.1'* ]] && truth_meta_v421=true
[[ "$regions_version" == *'HG002 v4.2.1'* ]] && regions_meta_v421=true
truth_name_v5=false; regions_name_v5=false; truth_name_v421=false; regions_name_v421=false
[[ "$truth_basename" == 'HG002_GRCh38_v5.0q_smvar.vcf.gz' ]] && truth_name_v5=true
[[ "$regions_basename" == 'HG002_GRCh38_v5.0q_smvar.benchmark.bed' ]] && regions_name_v5=true
[[ "$truth_basename" == *'v4.2.1'* ]] && truth_name_v421=true
[[ "$regions_basename" == *'v4.2.1'* ]] && regions_name_v421=true

if $truth_name_v5 && { $regions_name_v421 || $regions_meta_v421; }; then
    fail "HG002 v5.0q truth cannot be combined with v4.2.1 confident regions"
fi
if { $truth_name_v421 || $truth_meta_v421; } && { $regions_name_v5 || $regions_meta_v5; }; then
    fail "HG002 v4.2.1 truth cannot be combined with v5.0q confident regions"
fi
if $declared_v5 && { $truth_name_v421 || $truth_meta_v421 || $regions_name_v421 || $regions_meta_v421; }; then
    fail "HG002 v4.2.1 resources cannot be represented as benchmark version v5.0q"
fi
if $declared_v421 && { $truth_name_v5 || $truth_meta_v5 || $regions_name_v5 || $regions_meta_v5; }; then
    fail "HG002 v5.0q resources cannot be represented as benchmark version v4.2.1"
fi

primary_v5=false
if [[ "$resource_mode" == production && "$benchmark_mode" == germline_happy && \
      "$benchmark_name" == HG002_GIAB && "$benchmark_assembly" == GRCh38 ]] && $declared_v5; then
    $truth_meta_v5 && $regions_meta_v5 || \
        fail "HG002 v5.0q production identity requires matching trusted GIAB HG002 v5.0q manifest metadata for truth and confident regions"
    primary_v5=true
fi
if $truth_meta_v5 || $regions_meta_v5 || $truth_name_v5 || $regions_name_v5; then
    $primary_v5 || fail "HG002 v5.0q resources require production germline_happy identity HG002_GIAB / v5.0q / GRCh38"
fi

trusted_sex=''
trusted_sex_source=''
if $primary_v5; then
    trusted_sex=male
    trusted_sex_source='GIAB HG002 v5.0q benchmark identity and resource manifest'
elif [[ "$benchmark_name" == HG002_GIAB ]] && $declared_v421; then
    trusted_sex=male
    trusted_sex_source='declared GIAB HG002 v4.2.1 benchmark identity'
fi

benchmark_sex_resolved=none
sex_resolution_source='not_applicable'
hap_py_gender_argument='not_passed'
sex_aware_benchmarking=false
if [[ "$benchmark_mode" == germline_happy ]]; then
    case "$benchmark_sex" in
        auto)
            if $primary_v5; then
                benchmark_sex_resolved=male
                sex_resolution_source="$trusted_sex_source"
            else
                sex_resolution_source='legacy_semantics_preserved'
            fi
            ;;
        none)
            sex_resolution_source='explicit_none'
            ;;
        male|female)
            if [[ -n "$trusted_sex" && "$benchmark_sex" != "$trusted_sex" ]]; then
                fail "benchmark_sex $benchmark_sex conflicts with trusted $trusted_sex benchmark metadata ($trusted_sex_source)"
            fi
            benchmark_sex_resolved=$benchmark_sex
            sex_resolution_source='explicit_declaration'
            ;;
    esac
    if [[ "$benchmark_sex_resolved" == male || "$benchmark_sex_resolved" == female ]]; then
        hap_py_gender_argument=$benchmark_sex_resolved
        sex_aware_benchmarking=true
    fi
fi
if $primary_v5 && [[ "$benchmark_sex_resolved" != male ]]; then
    fail "HG002 v5.0q production benchmarking requires male-aware hap.py semantics; benchmark_sex '$benchmark_sex' resolved to '$benchmark_sex_resolved'"
fi
stratifications_version=''
stratifications_source=''
stratifications_checksum=''
if [[ -n "$benchmark_stratifications" ]]; then
    stratifications_version=$(manifest_value "$benchmark_stratifications" 4)
    stratifications_source=$(manifest_value "$benchmark_stratifications" 5)
    stratifications_checksum=$(manifest_value "$benchmark_stratifications" 7)
fi

case "$benchmark_mode" in
    germline_happy)
        engine='hap.py/vcfeval'
        engine_version='hap.py 0.3.15; RTG Tools 3.12.1'
        ;;
    regression|somatic)
        engine='normalized_exact_allele'
        engine_version='1.0 (bcftools 1.20)'
        ;;
esac

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$benchmark_name" "$benchmark_version" "$benchmark_mode" "$benchmark_assembly" \
    "$benchmark_sex" "$benchmark_sex_resolved" "$sex_resolution_source" "$hap_py_gender_argument" "$sex_aware_benchmarking" \
    "$engine" "$engine_version" "$engine_container" \
    "$truth_vcf" "$truth_vcf_index" "$truth_version" "$truth_source" "$truth_checksum" "$truth_index_checksum" \
    "$benchmark_regions" "$regions_version" "$regions_source" "$regions_checksum" \
    "$benchmark_stratifications" "$stratifications_version" "$stratifications_source" "$stratifications_checksum" \
    'strata may overlap; membership is evaluated independently and counts must not be summed as a partition' \
    'merged BED union intersect benchmark-region merged union; no within-stratum double-counting' \
    >> benchmark_config_manifest.tsv

emit_resource() {
    local name=$1 type=$2 path=$3 required_for=$4
    local version source checksum
    version=$(manifest_value "$path" 4)
    source=$(manifest_value "$path" 5)
    checksum=$(manifest_value "$path" 7)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$name" "$type" "$benchmark_assembly" "$version" "$source" "$path" "$checksum" "$required_for" "$resource_mode" \
        >> resource_manifest.tsv
}

emit_resource benchmark_truth_vcf truth_vcf "$truth_vcf" benchmarking
emit_resource benchmark_truth_vcf_index tabix_index "$truth_vcf_index" benchmarking
emit_resource benchmark_regions benchmark_regions_bed "$benchmark_regions" benchmarking
if [[ -n "$benchmark_stratifications" ]]; then
    emit_resource benchmark_stratifications stratification_tsv "$benchmark_stratifications" benchmarking
    stratification_dir=$(cd "$(dirname "$benchmark_stratifications")" && pwd)
    while IFS=$'\t' read -r region_name bed_path extra; do
        [[ -z "$region_name$bed_path$extra" || "$region_name" == \#* ]] && continue
        if [[ "$bed_path" = /* ]]; then resolved_bed=$bed_path; else resolved_bed="$stratification_dir/$bed_path"; fi
        declared_assembly=$(manifest_value "$resolved_bed" 3)
        [[ -z "$declared_assembly" || "$declared_assembly" == "$benchmark_assembly" ]] || \
            fail "stratification '$region_name' declares assembly '$declared_assembly' in the resource manifest, expected '$benchmark_assembly'"
        emit_resource "stratum_${region_name}" stratification_bed "$resolved_bed" stratified_benchmarking
    done < "$benchmark_stratifications"
fi

echo "BENCHMARK_CHECK: valid $benchmark_mode configuration for $analysis_type/$caller"
echo "BENCHMARK_CHECK: assembly checks use declared metadata and detectable contig/length compatibility; they do not prove biological resource identity"
