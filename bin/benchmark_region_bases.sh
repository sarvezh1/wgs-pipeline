#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: benchmark_region_bases.sh BENCHMARK_REGIONS.bed[.gz]" >&2
    exit 2
fi

bed=$1
if [[ "$bed" == *.gz ]]; then
    gzip -dc "$bed"
else
    command cat "$bed"
fi | awk 'BEGIN { OFS="\t" }
    /^#/ || /^track([[:space:]]|$)/ || /^browser([[:space:]]|$)/ { next }
    { print $1, $2, $3 }
' | LC_ALL=C sort -k1,1 -k2,2n -k3,3n | awk '
    NR == 1 { contig=$1; start=$2; end=$3; next }
    $1 != contig || $2 > end {
        bases += end-start
        contig=$1; start=$2; end=$3
        next
    }
    $3 > end { end=$3 }
    END {
        if (NR > 0) bases += end-start
        print bases+0
    }
'
