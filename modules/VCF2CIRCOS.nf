process VCF2CIRCOS_PREP {
    label 'process_low'
    tag "${meta.patient_id ?: meta.sample_id}:${caller_name}:${input_kind}"

    publishDir "${params.outdir ?: 'results'}/vcf2circos", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'

    input:
    tuple val(meta), val(caller_name), path(input_vcf), path(input_vcf_index), val(input_kind), val(analysis_type), val(calling_interval)

    output:
    tuple val(meta), val(caller_name), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.vcf2circos.input.vcf.gz"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.vcf2circos.input.vcf.gz.tbi"), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.vcf2circos_summary.tsv"), val(input_kind), val(calling_interval), emit: prepared

    script:
    def output_id = meta.patient_id ?: meta.sample_id
    """
    set -euo pipefail

    test -s '${input_vcf}'
    test -s '${input_vcf_index}'

    total_records=\$(bcftools view -H '${input_vcf}' | wc -l)
    pass_records=\$(bcftools view -f PASS -H '${input_vcf}' | wc -l)
    unfiltered_records=\$(bcftools view --include 'FILTER="."' -H '${input_vcf}' | wc -l)

    if [[ "\${pass_records}" -gt 0 ]]; then
        selection='PASS'
        bcftools view \
            --apply-filters PASS \
            --output-type z \
            --output '${output_id}.${caller_name}.vcf2circos.input.vcf.gz' \
            '${input_vcf}'
    elif [[ "\${total_records}" -gt 0 && "\${unfiltered_records}" -eq "\${total_records}" ]]; then
        selection='ALL_UNFILTERED'
        bcftools view \
            --output-type z \
            --output '${output_id}.${caller_name}.vcf2circos.input.vcf.gz' \
            '${input_vcf}'
    else
        echo 'VCF2Circos selection error: input has records but none are PASS, and not all records are unfiltered' >&2
        exit 1
    fi

    bcftools index \
        --tbi \
        --force \
        '${output_id}.${caller_name}.vcf2circos.input.vcf.gz'

    visualized_records=\$(bcftools view -H '${output_id}.${caller_name}.vcf2circos.input.vcf.gz' | wc -l)
    snvs=\$(bcftools view -v snps -H '${output_id}.${caller_name}.vcf2circos.input.vcf.gz' | wc -l)
    indels=\$(bcftools view -v indels -H '${output_id}.${caller_name}.vcf2circos.input.vcf.gz' | wc -l)

    printf 'caller\tanalysis_type\tsample_or_patient\tinput_kind\tselection\ttotal_input_records\tpass_input_records\tunfiltered_input_records\tvisualized_records\tsnvs\tindels\n' \
        > '${output_id}.${caller_name}.vcf2circos_summary.tsv'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        '${caller_name}' \
        '${analysis_type}' \
        '${output_id}' \
        '${input_kind}' \
        "\${selection}" \
        "\${total_records}" \
        "\${pass_records}" \
        "\${unfiltered_records}" \
        "\${visualized_records}" \
        "\${snvs}" \
        "\${indels}" \
        >> '${output_id}.${caller_name}.vcf2circos_summary.tsv'

    test "\${visualized_records}" -gt 0
    test -s '${output_id}.${caller_name}.vcf2circos.input.vcf.gz'
    test -s '${output_id}.${caller_name}.vcf2circos.input.vcf.gz.tbi'
    test -s '${output_id}.${caller_name}.vcf2circos_summary.tsv'
    """
}

process VCF2CIRCOS {
    label 'process_medium'
    tag "${meta.patient_id ?: meta.sample_id}:${caller_name}:${input_kind}"

    publishDir "${params.outdir ?: 'results'}/vcf2circos", mode: 'copy', overwrite: true, pattern: '*.html'
    container 'quay.io/biocontainers/vcf2circos:1.2.0--pyhdfd78af_0'

    input:
    tuple val(meta), val(caller_name), path(visualization_vcf), path(visualization_vcf_index), path(visualization_summary), val(input_kind), val(calling_interval), path(static_resources)

    output:
    tuple val(meta), val(caller_name), path("${meta.patient_id ?: meta.sample_id}.${caller_name}.vcf2circos.html"), path(visualization_summary), emit: visualizations

    script:
    def output_id = meta.patient_id ?: meta.sample_id
    def chromosome_config = calling_interval ? "\"list\": [\"${calling_interval.tokenize(':')[0]}\"],\n            \"all\": false" : '"list": [],\n            "all": true'
    def interval_refinement = calling_interval ? """
    interval='${calling_interval}'
    interval_chromosome=\"\${interval%%:*}\"
    interval_range=\"\${interval#*:}\"
    interval_start=\"\${interval_range%-*}\"
    interval_end=\"\${interval_range#*-}\"
    chromosome_length=\$(awk -v chromosome=\"\${interval_chromosome}\" '
        \$1 == chromosome { print (NF >= 3 ? \$3 : \$2); exit }
    ' '${static_resources}/Assembly/hg38/chr.hg38.sorted.txt')
    [[ \"\${chromosome_length}\" =~ ^[0-9]+\$ ]] || {
        echo \"VCF2Circos configuration error: chromosome '\${interval_chromosome}' was not found in the hg38 assembly resources\" >&2
        exit 1
    }

    mv '${output_id}.${caller_name}.vcf2circos.html' vcf2circos.unrefined.html
    refine_vcf2circos.py \\
        --input vcf2circos.unrefined.html \\
        --output '${output_id}.${caller_name}.vcf2circos.html' \\
        --sample '${output_id}' \\
        --caller '${caller_name}' \\
        --analysis-type \"\${analysis_type}\" \\
        --selection \"\${selection}\" \\
        --chromosome \"\${interval_chromosome}\" \\
        --start \"\${interval_start}\" \\
        --end \"\${interval_end}\" \\
        --chromosome-length \"\${chromosome_length}\" \\
        --records \"\${records}\" \\
        --snvs \"\${snvs}\" \\
        --indels \"\${indels}\"
    """ : ''
    """
    set -euo pipefail

    test -s '${visualization_vcf}'
    test -s '${visualization_vcf_index}'
    test -s '${visualization_summary}'

    cat > vcf2circos.options.json <<JSON
    {
        "General": {
            "title": "${output_id} variants | ${caller_name}",
            "titlefont": {
                "family": "Arial",
                "size": 25,
                "color": "#111827"
            },
            "width": 1200,
            "height": 1000,
            "plot_bgcolor": "white"
        },
        "Static": "\$PWD/${static_resources}",
        "Assembly": "hg38",
        "Chromosomes": {
            ${chromosome_config}
        },
        "Genes": {
            "only_snv_in_sv_genes": false,
            "extend": false
        },
        "Variants": {
            "annotations": {
                "fields": ["CSQ", "Gene_name", "TLOD", "NLOD", "AF", "DP", "MQ"]
            },
            "rings": {
                "position": 0.68,
                "height": 0.022,
                "space": 0.006,
                "nrings": 6
            }
        },
        "Extra": [],
        "Color": {
            "INS": "mediumorchid",
            "INV": "darkorange",
            "DEL": "crimson",
            "DUP": "royalblue",
            "CNV": "limegreen",
            "BND": "gold",
            "SNV": "#356FA8",
            "INDEL": "#D97706",
            "OTHER": "slategray",
            "MORBID_GENES": "firebrick",
            "GENES": "darkgray",
            "CYTOBAND": "gray",
            "INTERMEDIATE": "dimgrey",
            "GC": "darkblue",
            "MAPPABILITY": "indianred"
        }
    }
JSON

    vcf2circos \
        --input '${visualization_vcf}' \
        --output '${output_id}.${caller_name}.vcf2circos.html' \
        --options vcf2circos.options.json \
        --assembly hg38

    sed -i '/rel="icon"/d' \
        '${output_id}.${caller_name}.vcf2circos.html'

    analysis_type=\$(awk -F '\t' 'NR == 2 { print \$2 }' '${visualization_summary}')
    selection=\$(awk -F '\t' 'NR == 2 { print \$5 }' '${visualization_summary}')
    records=\$(awk -F '\t' 'NR == 2 { print \$9 }' '${visualization_summary}')
    snvs=\$(awk -F '\t' 'NR == 2 { print \$10 }' '${visualization_summary}')
    indels=\$(awk -F '\t' 'NR == 2 { print \$11 }' '${visualization_summary}')

    ${interval_refinement}

    test -s '${output_id}.${caller_name}.vcf2circos.html'
    grep -q 'Plotly.newPlot' '${output_id}.${caller_name}.vcf2circos.html'
    if [[ -n '${calling_interval}' ]]; then
        grep -q "\${interval_chromosome}" '${output_id}.${caller_name}.vcf2circos.html'
    fi
    grep -Eq '"name":"(SNV|INDEL|OTHER)"' '${output_id}.${caller_name}.vcf2circos.html'
    if grep -Eq '<script[^>]+src=|<link[^>]+href=' '${output_id}.${caller_name}.vcf2circos.html'; then
        echo 'VCF2Circos portability error: external HTML asset reference detected' >&2
        exit 1
    fi
    if grep -Eq '/tmp/|/work/' '${output_id}.${caller_name}.vcf2circos.html'; then
        echo 'VCF2Circos portability error: temporary or work path detected' >&2
        exit 1
    fi
    """
}
