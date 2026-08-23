#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: generate_resource_manifest.sh RESOURCE_DIR ASSEMBLY OUTPUT_TSV" >&2
    exit 2
fi

resource_dir=$1
assembly=$2
output=$3
[[ -d "$resource_dir" ]] || { echo "Resource directory not found: $resource_dir" >&2; exit 1; }
mkdir -p "$(dirname "$output")"

tmp="${output}.tmp.$$"
trap 'rm -f "$tmp"' EXIT
printf 'resource_name\tresource_type\tassembly\tversion\tsource\tpath_or_filename\tchecksum\trequired_for\n' > "$tmp"

while IFS= read -r path; do
    relative=${path#"$resource_dir"/}
    case "$relative" in
        genome/*.fa|genome/*.fasta) name=reference_fasta; type=reference_fasta; required=all ;;
        genome/*.fai) name=reference_fai; type=reference_index; required=all ;;
        genome/*.dict) name=reference_dict; type=sequence_dictionary; required=gatk ;;
        bwa_mem2/*) name=bwa_mem2_index; type=bwa_mem2_index; required=alignment ;;
        gatk/*.vcf.gz) name=bqsr_known_sites; type=known_sites_vcf; required=bqsr ;;
        gatk/*.tbi) name=bqsr_known_sites_index; type=tabix_index; required=bqsr ;;
        somatic/germline_resource.vcf.gz) name=mutect2_germline_resource; type=population_allele_frequency_vcf; required=mutect2 ;;
        somatic/germline_resource.vcf.gz.tbi) name=mutect2_germline_resource_index; type=tabix_index; required=mutect2 ;;
        somatic/panel_of_normals.vcf.gz) name=mutect2_panel_of_normals; type=panel_of_normals_vcf; required=mutect2 ;;
        somatic/panel_of_normals.vcf.gz.tbi) name=mutect2_panel_of_normals_index; type=tabix_index; required=mutect2 ;;
        somatic/contamination_sites.vcf.gz) name=mutect2_contamination_sites; type=population_sites_vcf; required=contamination ;;
        somatic/contamination_sites.vcf.gz.tbi) name=mutect2_contamination_sites_index; type=tabix_index; required=contamination ;;
        benchmarking/truth.vcf.gz) name=benchmark_truth_vcf; type=truth_vcf; required=benchmarking ;;
        benchmarking/truth.vcf.gz.tbi) name=benchmark_truth_vcf_index; type=tabix_index; required=benchmarking ;;
        benchmarking/confident_regions.bed|benchmarking/confident_regions.bed.gz) name=benchmark_regions; type=benchmark_regions_bed; required=benchmarking ;;
        benchmarking/stratifications.tsv) name=benchmark_stratifications; type=stratification_tsv; required=benchmarking ;;
        *) continue ;;
    esac
    checksum=$(sha256sum "$path" | awk '{print $1}')
    printf '%s\t%s\t%s\t\t\t%s\tsha256:%s\t%s\n' \
        "$name" "$type" "$assembly" "$relative" "$checksum" "$required" >> "$tmp"
done < <(find "$resource_dir" -type f -print | LC_ALL=C sort)

mv "$tmp" "$output"
trap - EXIT
echo "Wrote $output. Add verified version/source metadata where known; do not guess."
