process FIXTURE_REFERENCE_PREP {
    label 'process_low'
    tag "fixture reference preparation"

    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

    input:
    path reference_gz

    output:
    tuple path('reference.fa'), path('reference.fa.fai'), path('reference.dict'), emit: reference

    script:
    """
    set -euo pipefail
    gzip -dc '${reference_gz}' > reference.fa
    samtools faidx reference.fa
    awk 'BEGIN { print "@HD\tVN:1.6\tSO:unsorted" } { print "@SQ\tSN:" \$1 "\tLN:" \$2 }' \
        reference.fa.fai > reference.dict
    """
}

process FIXTURE_BWA_INDEX {
    label 'process_high'
    tag "fixture BWA-MEM2 indexing"

    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    input:
    tuple path(reference_fa), path(reference_fai), path(reference_dict)

    output:
    tuple path('reference.fa'), path('reference.fa.fai'), path('reference.dict'), path('bwa_index/reference.fa.*'), val('reference.fa'), emit: bundle

    script:
    """
    set -euo pipefail
    mkdir bwa_index
    cp reference.fa bwa_index/reference.fa
    bwa-mem2 index bwa_index/reference.fa
    """
}
