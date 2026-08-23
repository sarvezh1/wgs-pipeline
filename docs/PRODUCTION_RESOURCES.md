# Production resources: GRCh38 freeze

This document freezes the intended external scientific ecosystem for human short-read GRCh38 WGS. It does not bundle the biological files, authorize a production run, or replace study-specific review. The repository's compact chr20 resources remain technical regression fixtures; they are not substitutes for the resources below.

The machine-readable starting contract is [`resources/GRCh38.production-resource-manifest.template.tsv`](../resources/GRCh38.production-resource-manifest.template.tsv). Populate its blank SHA-256 fields only after obtaining or generating the exact files. Never copy an upstream MD5 value into a `sha256:` field.

## External tree contract

```text
GRCh38/
├── genome/
│   ├── GRCh38.fa
│   ├── GRCh38.fa.fai
│   └── GRCh38.dict
├── bwa_mem2/
│   ├── GRCh38.0123
│   ├── GRCh38.amb
│   ├── GRCh38.ann
│   ├── GRCh38.bwt.2bit.64
│   ├── GRCh38.pac
│   └── GRCh38.alt
├── gatk/
│   ├── Homo_sapiens_assembly38.dbsnp138.vcf.gz
│   ├── Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi
│   ├── Homo_sapiens_assembly38.known_indels.vcf.gz
│   ├── Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
│   ├── Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
│   └── Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
├── somatic/
│   ├── germline_resource.vcf.gz
│   ├── germline_resource.vcf.gz.tbi
│   ├── contamination_sites.vcf.gz
│   ├── contamination_sites.vcf.gz.tbi
│   ├── panel_of_normals.vcf.gz
│   └── panel_of_normals.vcf.gz.tbi
├── annotation/vep/homo_sapiens/116_GRCh38/
├── benchmarking/
│   ├── truth.vcf.gz
│   ├── truth.vcf.gz.tbi
│   ├── confident_regions.bed
│   └── stratifications/
│       ├── stratifications.tsv
│       └── GRCh38@all/...
├── cnv/
│   ├── intervals.interval_list
│   ├── germline/contig-ploidy-model/
│   ├── germline/gcnv-model/
│   └── somatic/
│       ├── read-count-pon.hdf5
│       ├── common-sites.vcf.gz
│       └── common-sites.vcf.gz.tbi
└── metadata/manifest.tsv
```

`truth.vcf.gz`, `confident_regions.bed`, and `stratifications.tsv` are stable logical names for one explicitly selected benchmark release. Keep immutable source filenames alongside them or use auditable links. Do not pair files from different GIAB releases.

## 1. Reference and alignment index

Freeze the Broad/GATK hg38 resource-bundle `v0` reference:

- FASTA: `Homo_sapiens_assembly38.fasta`
- FAI: `Homo_sapiens_assembly38.fasta.fai`
- dictionary: `Homo_sapiens_assembly38.dict`
- source root: `gs://gcp-public-data--broad-references/hg38/v0/`
- assembly identity: GRCh38 base assembly, associated with `GCA_000001405.15`; this Broad analysis bundle is not the plain primary-assembly FASTA
- contigs: UCSC-style `chr1`–`chr22`, `chrX`, `chrY`, `chrM`, followed by unlocalized/unplaced, ALT, decoy, and other bundle sequences; the verified dictionary contains 3,366 sequences

Keep the exact FASTA bytes unchanged. Rename or link it to `genome/GRCh38.fa`, and ensure the FAI and dictionary describe those same bytes and contig order. Build a fresh BWA-MEM2 2.2.1 index from that FASTA:

```bash
bwa-mem2 index -p "${WGS_RESOURCE_ROOT}/bwa_mem2/GRCh38" "${WGS_RESOURCE_ROOT}/genome/GRCh38.fa"
cp "${WGS_BROAD_BUNDLE}/Homo_sapiens_assembly38.fasta.64.alt" "${WGS_RESOURCE_ROOT}/bwa_mem2/GRCh38.alt"
```

The Broad FASTA includes ALT contigs. BWA-MEM2 reads `<index-prefix>.alt`; the workflow now requires and stages that sidecar whenever `_alt` contigs are present in the FAI. Confirm the alignment log reports a non-zero ALT-contig count. The legacy Broad BWA index is not a BWA-MEM2 index and must not be renamed as one.

## 2. GATK BQSR known sites

The frozen known-sites selection is the compatible Broad hg38/v0 trio:

| Resource | Frozen identity | Intended use |
|---|---|---|
| `Homo_sapiens_assembly38.dbsnp138.vcf` | dbSNP build 138 in the Broad hg38 bundle | known polymorphic sites |
| `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | Mills and 1000 Genomes gold-standard indels, Broad hg38 bundle | known indels |
| `Homo_sapiens_assembly38.known_indels.vcf.gz` | Broad hg38 known-indels resource | additional known indels |

This is the pipeline's selected production bundle, not a claim that GATK mandates exactly these three files. GATK BaseRecalibrator requires one or more trustworthy known-sites resources and documents dbSNP, gnomAD, or other polymorphism resources as possible inputs. Every configured VCF must use the exact reference's contig convention and have a tabix index. Because the Broad dbSNP file is distributed as uncompressed VCF plus Tribble index, prepare the pipeline form without changing records:

```bash
bgzip -c Homo_sapiens_assembly38.dbsnp138.vcf > Homo_sapiens_assembly38.dbsnp138.vcf.gz
tabix -p vcf Homo_sapiens_assembly38.dbsnp138.vcf.gz
```

Record SHA-256 for the distributed source and the derived bgzipped file separately. BQSR resource choice masks known variation during model fitting; it is not a truth set or a post-calling variant filter.

## 3. HG002 germline benchmarks

### Primary publication benchmark: GIAB v5.0q

Use this matched GRCh38 pair from the NIST GIAB HG002 v5.0q directory:

- `HG002_GRCh38_v5.0q_smvar.vcf.gz`
- `HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi`
- `HG002_GRCh38_v5.0q_smvar.benchmark.bed`
- release README: `NIST_HG002_v5.0q_variant-benchmarksets_README.md`

The v5.0q small-variant VCF and benchmark BED are a pair. The release is assembly-based on HG002 Q100 v1.1, marked draft/under active development by NIST, and includes challenging regions absent from earlier mapping-based sets. NIST recommends hap.py with RTG vcfeval, GIAB v3.6 stratifications, and `--gender male` for HG002. Upstream MD5 values verified from `checksum.md5` are:

| File | Upstream MD5 |
|---|---|
| benchmark BED | `3366858af85875cbb3c579412ffa4c56` |
| small-variant VCF | `c71acc71069bf7cd7f51eb8fb0c1a1ab` |
| VCF tabix index | `58cab2a06e29b74a5bdc190ead079083` |

The workflow pins hap.py 0.3.15 and RTG Tools 3.12.1. Its verified hap.py CLI accepts `--gender {male,female,auto,none}`. Configure the primary run as `--benchmark_name HG002_GIAB --benchmark_version v5.0q --benchmark_sex auto` (or explicit `male`): trusted matching v5.0q truth/region manifest metadata resolve this to `--gender male`. `none`, contradictory metadata, a female declaration, or a v4.2.1/v5.0q resource mixture fails before benchmarking. This interface support does not itself constitute full-WGS validation.

### Secondary validation: GIAB v4.2.1

Retain v4.2.1 for continuity with the compact regression baseline:

- `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`
- `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi`
- `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`

NIST now deprecates v4.2.1 for HG002 in favor of v5.0q. Label it secondary/legacy validation, not the primary current publication benchmark. The bundled chr20 subset remains a technical fixture and must not be mixed with full v5.0q regions.

## 4. GIAB/GA4GH stratifications

Freeze GIAB genome stratifications v3.6, GRCh38:

- package: `genome-stratifications-GRCh38@all.tar.gz`
- extracted root: `GRCh38@all/`
- pipeline/hap.py manifest: `GRCh38@all/GRCh38-all-stratifications.tsv`
- upstream package MD5: `f6b0666987b01aec08e9f5fa17707b10`
- upstream TSV MD5: `05f4b5d0e945e9173a9eef39b8e5dbbe`

The official TSV is already `stratum_name<TAB>relative_BED_path`, matching the pipeline stratification interface. Preserve its authoritative names and directory layout. BEDs are bgzip-compressed, may overlap, and are non-additive. Use upstream MD5 to validate acquisition, then calculate SHA-256 for the package, TSV, and every BED actually used and record those SHA-256 values in the production manifest. Do not commit the roughly 1.3 GB package.

## 5. Ensembl VEP

Freeze Ensembl VEP 116.0 with the release-116 `homo_sapiens` GRCh38 cache:

- archive: `homo_sapiens_vep_116_GRCh38.tar.gz`
- source: `https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/`
- extracted path required by the workflow: `annotation/vep/homo_sapiens/116_GRCh38/`
- cache assembly metadata: GRCh38.p14; primary chromosome coordinates are unchanged from the GRCh38 base assembly

Cache mode runs `--cache --offline --cache_version 116 --assembly GRCh38 --fasta` and needs no network. The supplied Broad FASTA must contain every input contig; VEP's chromosome-synonym support does not justify mixing assemblies. Optional GTF mode accepts a bgzip/tabix-indexed, explicitly versioned GRCh38 GTF, but it is not equivalent to the VEP cache: cache-only identifiers, colocated-variant content, transcript selection, and other annotations may be absent. Record the exact GTF source/release if used.

## 6. Mutect2 matched tumor-normal resources

Preserve the implemented path: Mutect2 → F1R2 → LearnReadOrientationModel → GetPileupSummaries → CalculateContamination → FilterMutectCalls. No thresholds change in this freeze.

Downloadable/reference inputs:

- germline population AF resource: `af-only-gnomad.hg38.vcf.gz` plus `.tbi` from `gs://gatk-best-practices/somatic-hg38/`
- pileup/common-sites resource: `small_exac_common_3.hg38.vcf.gz` plus `.tbi` from the same source

Study/platform-specific input:

- `panel_of_normals.vcf.gz` plus `.tbi`, generated from normal samples processed with the study's sequencing platform, preparation, reference, and Mutect2 version, then combined using GATK CreateSomaticPanelOfNormals

GATK describes the germline resource and PoN as recommended rather than inherently required by Mutect2. This pipeline intentionally requires both in its named production profile. A public/test PoN may be useful for workflow validation, but is not frozen as a universally defensible production PoN. Record normal cohort, platform, library preparation, reference SHA-256, GATK/container identity, generation command, date, and tree/file checksum.

## 7. Structural variants

Manta 1.6.0 requires coordinate-sorted/indexed paired-end BAM or CRAM files and the identical indexed reference FASTA used for mapping. It does not require a population database, PoN, or other external biological resource for unrestricted WGS. An optional bgzip/tabix-indexed call-regions BED is used only for deliberately restricted runs.

Publication SV benchmarking is a separate future validation activity. GIAB v5.0q provides paired GRCh38 `stvar` VCF/BED files and recommends complex-representation-aware methods such as Truvari with refine. Do not route those files into the current small-variant benchmark or add an SV benchmark engine during this freeze.

## 8. GATK CNV resources

### Germline CASE mode

- generic/reference-derived: one preprocessed GRCh38 interval list, created with a documented bin size and filtering procedure from the exact FASTA
- cohort/model-derived: DetermineGermlineContigPloidy model trained in COHORT mode with explicit contig-ploidy priors
- cohort/model-derived: matching GermlineCNVCaller model shard(s), trained on the same interval definition, reference, platform, library preparation, and broadly comparable coverage regime

No universal germline gCNV model is frozen. The model directories and interval list are a coupled set. Record cohort inclusion criteria, model/tool versions, parameters, reference and interval checksums, shard layout, and a deterministic tree SHA-256.

### Matched somatic mode

- generic/reference-derived: the exact same preprocessed interval list used for case and PoN counts
- study/platform-specific: read-count denoising PoN built from compatible normals using CreateReadCountPanelOfNormals
- population/reference input: indexed GRCh38 common-sites VCF used by CollectAllelicCounts
- sample input: matched-normal allelic counts, produced by the workflow from the matched normal

There is no universal somatic denoising PoN. The Broad public test PoN is suitable for plumbing examples only unless its cohort and assay are demonstrably compatible. Common-sites selection must be documented and consistent across tumor and normal; the current interface expects a bgzipped VCF plus `.tbi`.

## 9. Sex chromosomes and pseudoautosomal regions

Use the NCBI Genome Reference Consortium GRCh38.p14 assembly page as the authoritative coordinate source. It reports 1-based inclusive coordinates:

| Region | X | Y |
|---|---|---|
| PAR1 | 10,001–2,781,479 | 10,001–2,781,479 |
| PAR2 | 155,701,383–156,030,895 | 56,887,903–57,217,415 |

The corresponding zero-based, half-open BED intervals are `chrX 10000 2781479`, `chrX 155701382 156030895`, `chrY 10000 2781479`, and `chrY 56887902 57217415`. Store both the source coordinate convention and the conversion command/checksum if a future BED is created. These coordinates are provenance for later validation only: current X/Y/PAR caller biology and per-region ploidy routing remain unvalidated and disabled.

## Acquisition and verification policy

1. Obtain files only from the authoritative roots named above; preserve upstream README and checksum listings.
2. Verify the upstream checksum algorithm exactly as published.
3. Calculate `sha256sum` for every downloaded or derived file after final naming/indexing.
4. For model/cache directories, create a sorted relative-path plus per-file-SHA-256 inventory and hash that inventory; preserve both.
5. Populate `metadata/manifest.tsv` with exact filenames, versions, sources, SHA-256, and intended stages. Blank means not locally verified; never use a guessed checksum.
6. Confirm FASTA/FAI/dictionary identity, contig names and lengths, VCF/BED contigs, and index readability before analysis.
7. Keep raw full-GRCh38, VEP cache, GIAB stratifications, cohort models, and PoNs outside the repository.

## Attribution and release boundary

The pipeline is MIT-licensed; third-party software and data retain their own terms. Cite and attribute the Broad Institute/GATK resource bundle and tools, NIST GIAB and GA4GH benchmarking/stratification resources, Ensembl VEP, Manta, mosdepth, VCF2Circos (AGPL-3.0-or-later), and embedded Plotly.js (MIT) as applicable. Preserve upstream NIST data-use notices. Do not infer a DOI, ORCID, or publication citation that upstream has not assigned.

This resource freeze does not establish clinical validity, full-WGS accuracy or performance, SV sensitivity, CNV sensitivity, or X/Y/PAR correctness. Those require separate production-scale experiments after the v5.0q male-aware benchmark interface is resolved.
