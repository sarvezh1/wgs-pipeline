#!/usr/bin/env bash
set -euo pipefail

fail() {
    echo "RESOURCE_CHECK ERROR: $1" >&2
    exit 1
}

mode=$1
assembly=$2
project_dir=$3
reference=$4
reference_fai=$5
reference_dict=$6
bwa_prefix=$7
known_sites_joined=$8
manifest=$9
run_bqsr=${10}
run_vep=${11}
vep_mode=${12}
vep_cache=${13}
vep_gtf=${14}
run_vcf2circos=${15}
vcf2circos_static=${16}
analysis_type=${17}
caller=${18}
gatk_germline_mode=${19}
gatk_somatic_profile=${20}
mutect2_germline_resource=${21}
mutect2_pon=${22}
mutect2_contamination_sites=${23}

[[ "$mode" == fixture || "$mode" == production ]] || \
    fail "unsupported resource_mode '$mode'. Expected fixture or production"
[[ -n "$assembly" ]] || fail "assembly must not be empty"
[[ -r "$reference" && -s "$reference" ]] || fail "missing or unreadable reference FASTA: '$reference'"

fixture_root="$project_dir/test_data/reference"
case "$mode" in
    fixture)
        [[ "$reference" == "$fixture_root/"* ]] || \
            fail "fixture mode requires a bundled reference below '$fixture_root'"
        [[ "$reference" == *.gz ]] || fail "fixture reference must be the bundled gzipped FASTA"
        ;;
    production)
        [[ "$reference" != "$project_dir/test_data/"* ]] || \
            fail "production mode cannot use a repository test_data reference; configure external production resources"
        [[ "$reference" != *.gz ]] || \
            fail "production reference must be an uncompressed FASTA with reusable .fai, dictionary, and BWA-MEM2 indexes"
        [[ "$(head -c 1 "$reference")" == '>' ]] || fail "production reference does not begin with a FASTA header: '$reference'"
        for item in "$reference_fai" "$reference_dict"; do
            [[ -r "$item" && -s "$item" ]] || fail "missing or unreadable production reference companion: '$item'"
        done
        [[ "$reference_fai" == "${reference}.fai" ]] || \
            fail "reference_fai must be the FASTA companion '${reference}.fai'"
        reference_base=${reference%.*}
        [[ "$reference_dict" == "${reference_base}.dict" ]] || \
            fail "reference_dict must be the GATK companion '${reference_base}.dict'"

        expected_suffixes=(.0123 .amb .ann .bwt.2bit.64 .pac)
        for suffix in "${expected_suffixes[@]}"; do
            [[ -r "${bwa_prefix}${suffix}" && -s "${bwa_prefix}${suffix}" ]] || \
                fail "incomplete BWA-MEM2 index: missing '${bwa_prefix}${suffix}'"
        done

        if awk '$1 ~ /_alt$/ { found=1 } END { exit(found ? 0 : 1) }' "$reference_fai"; then
            [[ -r "${bwa_prefix}.alt" && -s "${bwa_prefix}.alt" ]] || \
                fail "reference contains ALT contigs but BWA-MEM2 ALT metadata is missing: '${bwa_prefix}.alt'"
        fi

        awk 'NF < 5 || $2 !~ /^[0-9]+$/ { exit 1 } END { if (NR == 0) exit 1 }' "$reference_fai" || \
            fail "malformed FASTA index: '$reference_fai'"
        awk -F '\t' '
            /^@SQ/ {
                sn=""; ln=""
                for (i=2; i<=NF; i++) {
                    if ($i ~ /^SN:/) sn=substr($i,4)
                    if ($i ~ /^LN:/) ln=substr($i,4)
                }
                if (sn == "" || ln !~ /^[0-9]+$/) exit 1
                print sn "\t" ln
            }
        ' "$reference_dict" > .dict_sequences.tsv || fail "malformed sequence dictionary: '$reference_dict'"
        [[ -s .dict_sequences.tsv ]] || fail "sequence dictionary contains no @SQ records: '$reference_dict'"
        cut -f1,2 "$reference_fai" > .fai_sequences.tsv
        cmp -s .fai_sequences.tsv .dict_sequences.tsv || \
            fail "reference identity mismatch: FASTA index and sequence dictionary contigs/lengths differ"

        fai_length=$(awk '{ total += $2 } END { print total+0 }' "$reference_fai")
        fai_contigs=$(awk 'END { print NR+0 }' "$reference_fai")
        read -r bwa_length bwa_contigs _ < "${bwa_prefix}.ann" || fail "cannot read BWA-MEM2 .ann metadata"
        [[ "$bwa_length" == "$fai_length" && "$bwa_contigs" == "$fai_contigs" ]] || \
            fail "BWA-MEM2 index metadata does not match FASTA index total length/contig count"
        ;;
esac

expected_manifest_header=$'resource_name\tresource_type\tassembly\tversion\tsource\tpath_or_filename\tchecksum\trequired_for'
[[ -r "$manifest" && -s "$manifest" ]] || fail "missing or unreadable resource manifest: '$manifest'"
[[ "$(head -n 1 "$manifest")" == "$expected_manifest_header" ]] || \
    fail "malformed manifest header in '$manifest'; expected: $expected_manifest_header"
awk -F '\t' '
    NR == 1 { next }
    NF != 8 || $1 == "" || $2 == "" || $6 == "" { exit 1 }
    $7 != "" && $7 !~ /^sha256:[0-9a-fA-F]{64}$/ { exit 1 }
    END { if (NR < 2) exit 1 }
' "$manifest" || fail "malformed manifest row or checksum in '$manifest'"
if [[ "$mode" == production ]]; then
    awk -F '\t' -v expected="$assembly" 'NR > 1 && $3 != "" && $3 != expected { exit 1 }' "$manifest" || \
        fail "manifest assembly metadata conflicts with configured assembly '$assembly'"
else
    manifest_dir=$(cd "$(dirname "$manifest")" && pwd)
    while IFS=$'\t' read -r manifest_path checksum; do
        fixture_path="$manifest_dir/$manifest_path"
        [[ -r "$fixture_path" ]] || fail "fixture manifest resource is missing: '$fixture_path'"
        observed="sha256:$(sha256sum "$fixture_path" | awk '{print $1}')"
        [[ "$observed" == "$checksum" ]] || fail "fixture checksum mismatch for '$fixture_path'"
    done < <(awk -F '\t' 'NR > 1 && $7 != "" { print $6 "\t" $7 }' "$manifest")
fi

IFS='|' read -r -a known_sites <<< "$known_sites_joined"
if [[ "$run_bqsr" == true ]]; then
    ((${#known_sites[@]} > 0)) && [[ -n "${known_sites[0]}" ]] || \
        fail "BQSR requested but no known-sites resources were configured"
    for known in "${known_sites[@]}"; do
        [[ -r "$known" && -s "$known" ]] || fail "missing or unreadable BQSR known-sites VCF: '$known'"
        [[ -r "${known}.tbi" && -s "${known}.tbi" ]] || fail "missing BQSR known-sites index: '${known}.tbi'"
        if [[ "$mode" == production && "$known" == "$project_dir/test_data/"* ]]; then
            fail "production BQSR cannot use a repository fixture known-sites VCF: '$known'"
        fi
        if [[ "$mode" == production ]]; then
            gzip -dc "$known" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' > .vcf_contigs.tsv
            if [[ -s .vcf_contigs.tsv ]] && ! awk 'NR==FNR { ref[$1]=1; next } ref[$1] { found=1 } END { exit(found ? 0 : 1) }' \
                "$reference_fai" .vcf_contigs.tsv; then
                fail "known-sites VCF and production reference have no shared declared contigs; check assembly and contig naming"
            fi
        fi
    done
fi

if [[ "$run_vep" == true ]]; then
    case "$vep_mode" in
        cache)
            [[ "$assembly" == GRCh38 ]] || fail "VEP release-116 cache mode is configured for GRCh38, not '$assembly'"
            required_cache="$vep_cache/homo_sapiens/116_GRCh38"
            [[ -d "$required_cache" ]] || fail "VEP cache mode requires release-116 GRCh38 cache at '$required_cache'"
            find "$required_cache" -mindepth 1 -print -quit | grep -q . || fail "VEP cache directory is empty: '$required_cache'"
            ;;
        gtf)
            [[ -r "$vep_gtf" && -s "$vep_gtf" ]] || fail "VEP GTF mode requires a readable bgzipped GTF: '$vep_gtf'"
            [[ -r "${vep_gtf}.tbi" && -s "${vep_gtf}.tbi" ]] || fail "VEP GTF mode requires tabix index '${vep_gtf}.tbi'"
            ;;
        *) fail "vep_mode must be 'cache' or 'gtf'" ;;
    esac
fi

if [[ "$run_vcf2circos" == true ]]; then
    [[ "$assembly" == GRCh38 ]] || fail "VCF2Circos static resources are configured for hg38/GRCh38, not '$assembly'"
    required_static=(
        morbid.txt
        Assembly/hg38/chr.hg38.sorted.txt
        Assembly/hg38/cytoband_hg38_chr_infos.txt.gz
        Assembly/hg38/genes.hg38.sorted.txt
        Assembly/hg38/transcripts.hg38.sorted.txt
        Assembly/hg38/exons.hg38.sorted.txt
    )
    for relative in "${required_static[@]}"; do
        [[ -r "$vcf2circos_static/$relative" && -s "$vcf2circos_static/$relative" ]] || \
            fail "incomplete VCF2Circos hg38 resources: missing '$vcf2circos_static/$relative'"
    done
fi

validate_somatic_vcf() {
    local label=$1 vcf=$2 require_af=$3
    [[ -r "$vcf" && -s "$vcf" ]] || fail "missing or unreadable $label VCF: '$vcf'"
    [[ -r "${vcf}.tbi" && -s "${vcf}.tbi" ]] || fail "missing $label index: '${vcf}.tbi'"
    [[ "$vcf" != "$project_dir/test_data/"* ]] || \
        fail "production Mutect2 cannot use a repository fixture resource for $label: '$vcf'"
    first_header=$(awk 'NR == 1 { print; exit }' < <(gzip -dc "$vcf"))
    [[ "$first_header" == '##fileformat=VCF'* ]] || fail "$label is not a readable bgzipped VCF: '$vcf'"
    if [[ "$require_af" == true ]]; then
        awk '/^##INFO=<ID=AF,/ { found=1 } !/^#/ { exit } END { exit(found ? 0 : 1) }' < <(gzip -dc "$vcf") || \
            fail "$label VCF must declare the INFO/AF population allele-frequency field"
    fi
    awk '
        /^##contig=<ID=/ {
            value=$0
            sub(/^##contig=<ID=/, "", value)
            sub(/[,>].*$/, "", value)
            print value
        }
        !/^#/ { exit }
    ' < <(gzip -dc "$vcf") > ".${label//[^A-Za-z0-9]/_}.contigs.tsv"
    contig_file=".${label//[^A-Za-z0-9]/_}.contigs.tsv"
    if [[ -s "$contig_file" ]] && ! awk 'NR==FNR { ref[$1]=1; next } ref[$1] { found=1 } END { exit(found ? 0 : 1) }' \
        "$reference_fai" "$contig_file"; then
        fail "$label VCF and production reference have no shared declared contigs; check assembly and contig naming"
    fi
}

case "$gatk_germline_mode" in
    vcf|gvcf) ;;
    *) fail "unsupported gatk_germline_mode '$gatk_germline_mode'. Expected vcf or gvcf" ;;
esac
case "$gatk_somatic_profile" in
    basic|production) ;;
    *) fail "unsupported gatk_somatic_profile '$gatk_somatic_profile'. Expected basic or production" ;;
esac
if [[ "$gatk_germline_mode" == gvcf && "$analysis_type:$caller" != germline:haplotypecaller && "$analysis_type:$caller" != germline:compare ]]; then
    fail "gatk_germline_mode gvcf requires analysis_type germline and caller haplotypecaller or compare"
fi
if [[ "$gatk_somatic_profile" == production ]]; then
    [[ "$analysis_type:$caller" == somatic:mutect2 || "$analysis_type:$caller" == somatic:compare ]] || \
        fail "gatk_somatic_profile production requires analysis_type somatic and caller mutect2 or compare"
    [[ "$mode" == production ]] || \
        fail "gatk_somatic_profile production requires resource_mode production; fixture resources are not substituted"
    validate_somatic_vcf "Mutect2 germline resource" "$mutect2_germline_resource" true
    validate_somatic_vcf "Mutect2 panel of normals" "$mutect2_pon" false
    validate_somatic_vcf "Mutect2 contamination sites" "$mutect2_contamination_sites" true
fi

manifest_value() {
    local target=$1 column=$2
    awk -F '\t' -v target="$target" -v column="$column" '
        NR > 1 && ($6 == target || $6 == target_name || substr($6, length($6)-length(target_name)+1) == target_name) { print $column; exit }
    ' target_name="${target##*/}" "$manifest"
}

emit_resource() {
    local name=$1 type=$2 path=$3 required_for=$4
    local version source checksum
    version=$(manifest_value "$path" 4)
    source=$(manifest_value "$path" 5)
    checksum=$(manifest_value "$path" 7)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$name" "$type" "$assembly" "$version" "$source" "$path" "$checksum" "$required_for" "$mode" >> resource_manifest.tsv
}

printf 'resource_name\tresource_type\tassembly\tversion\tsource\tpath_or_filename\tchecksum\trequired_for\tresource_mode\n' > resource_manifest.tsv
emit_resource reference_fasta reference_fasta "$reference" all
if [[ "$mode" == production ]]; then
    emit_resource reference_fai reference_index "$reference_fai" all
    emit_resource reference_dict sequence_dictionary "$reference_dict" gatk
    for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
        emit_resource "bwa_mem2_index${suffix}" bwa_mem2_index "${bwa_prefix}${suffix}" alignment
    done
    if [[ -r "${bwa_prefix}.alt" && -s "${bwa_prefix}.alt" ]]; then
        emit_resource bwa_mem2_alt bwa_alt_metadata "${bwa_prefix}.alt" alignment
    fi
fi
if [[ "$run_bqsr" == true ]]; then
    for known in "${known_sites[@]}"; do
        emit_resource bqsr_known_sites known_sites_vcf "$known" bqsr
        emit_resource bqsr_known_sites_index tabix_index "${known}.tbi" bqsr
    done
fi
if [[ "$run_vep" == true ]]; then
    if [[ "$vep_mode" == cache ]]; then
        emit_resource vep_cache vep_cache "$vep_cache" annotation
    else
        emit_resource vep_gtf gtf "$vep_gtf" annotation
        emit_resource vep_gtf_index tabix_index "${vep_gtf}.tbi" annotation
    fi
fi
if [[ "$run_vcf2circos" == true ]]; then
    emit_resource vcf2circos_static vcf2circos_static "$vcf2circos_static" visualization
fi
if [[ "$gatk_somatic_profile" == production ]]; then
    emit_resource mutect2_germline_resource population_allele_frequency_vcf "$mutect2_germline_resource" mutect2
    emit_resource mutect2_germline_resource_index tabix_index "${mutect2_germline_resource}.tbi" mutect2
    emit_resource mutect2_panel_of_normals panel_of_normals_vcf "$mutect2_pon" mutect2
    emit_resource mutect2_panel_of_normals_index tabix_index "${mutect2_pon}.tbi" mutect2
    emit_resource mutect2_contamination_sites population_sites_vcf "$mutect2_contamination_sites" contamination
    emit_resource mutect2_contamination_sites_index tabix_index "${mutect2_contamination_sites}.tbi" contamination
fi

touch .resource_check.ok
echo "RESOURCE_CHECK: valid $mode resource configuration (assembly declared as $assembly)"
echo "RESOURCE_CHECK: assembly identity is checked from companion metadata where possible; the declared build remains user/resource-provider metadata"
