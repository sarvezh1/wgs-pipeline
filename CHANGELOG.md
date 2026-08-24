# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0-rc1] - 2026-08-23

### Added

- Nextflow DSL2 preprocessing for paired-end FASTQ validation, FastQC, fastp, BWA-MEM2 alignment, coordinate sorting, duplicate marking, BAM indexing, and samtools QC.
- Optional GATK BQSR with standardized caller-ready BAM routing.
- Configurable calling intervals, with unrestricted supplied-reference calling when no interval is specified.
- Germline calling with GATK HaplotypeCaller and DeepVariant 1.10.0.
- Matched tumor-normal somatic calling with Mutect2, Strelka2 2.9.10, and DeepSomatic 1.10.0.
- Shared bcftools normalization, indexing, statistics, and compact variant summaries.
- Optional Ensembl VEP 116 annotation and VCF2Circos 1.2.0 visualization.
- Explicit regional (compatibility default) and production genome VCF2Circos modes, with a repository-authored offline Plotly genome ideogram, binned PASS SNV/indel densities, optional native mosdepth/GATK CNV tracks, real final-Manta breakpoint links, and machine-readable track/provenance TSVs.
- Public genome-mode Circos demonstration with canonical GRCh38 ordering, concentric density/coverage/CNV tracks, optional final PASS Manta relationships, self-contained offline HTML, and preserved regional visualization compatibility.
- Optional MultiQC 1.27.1 aggregation of native QC reports.
- Compact real HG002 and HG008 GRCh38 chr20 validation fixtures, benchmark truth subsets, known-sites data, and provenance records.
- Self-contained HG008 MultiQC and VCF2Circos example reports with browser-rendered screenshots.
- Explicit fixture and production resource modes with external GRCh38 bundle routing, reusable prebuilt BWA-MEM2 indexes, reference-companion validation, BQSR fixture safeguards, resource manifests, and run-level resource provenance.
- Optional GATK HaplotypeCaller reference-confidence mode with standard indexed GVCF artifacts, GenotypeGVCFs integration, and final-VCF-only downstream routing.
- Strict Mutect2 production profile with external germline population and Panel of Normals resources, F1R2 collection, LearnReadOrientationModel, matched tumor-normal GetPileupSummaries/CalculateContamination, and enhanced FilterMutectCalls inputs.
- GATK path provenance recording the calling mode, pinned version, contamination/orientation usage, and filtering route.
- Optional automated benchmarking with separate fixture regression, germline hap.py/vcfeval, and somatic normalized exact-allele modes.
- Mandatory benchmark-region enforcement, early truth/index/BED/contig validation, optional native hap.py stratification TSV routing, pinned hap.py 0.3.15 and RTG Tools 3.12.1 execution, and preserved native engine outputs.
- Standardized overall/SNV/indel TP, FP, FN, precision, recall, and F1 summaries plus benchmark resource/configuration provenance.
- Optional `--caller compare` orchestration that fans one caller-ready germline BAM to HaplotypeCaller and DeepVariant, or one matched tumor-normal BAM pair to Mutect2, Strelka2, and DeepSomatic.
- Descriptive normalized allele membership, exclusive overlap categories, pairwise Jaccard and agreement matrices, germline diploid genotype concordance, truth-match fractions within benchmark regions, discordant-variant tables, UpSet-ready membership, provenance, and neutral shared/caller-specific derivative VCFs.
- Optional stratified biological benchmarking through a validated `stratum_name<TAB>bed_path` manifest shared by germline hap.py, somatic/exact benchmarks, and caller-comparison derivatives.
- Native hap.py subset TP/FP/FN extraction, GA4GH `BD` record-level error tables, somatic normalized PASS-allele strata, merged effective benchmark-base accounting, overlapping-stratum semantics, stratified truth-match fractions and pairwise Jaccard, germline genotype/no-call localization, and caller-independent normalized indel-size bins.
- External production GIAB/GA4GH stratification resource contracts plus clearly labeled overlapping HG002/HG008 coordinate-window fixtures, per-BED checksums/provenance, explicit `NA` denominator behavior, and actionable negative validation for malformed/incompatible strata.
- Optional native Nextflow trace/report/timeline/DAG generation below the selected output directory, with standardized task, stage, caller, compute-component, workflow, and non-identifying performance-provenance TSVs.
- Explicit fresh-versus-cached task accounting, task-compute versus workflow-elapsed semantics, technical-fixture performance labels, and negative trace/configuration validation.
- Low/medium/high process labels, unchanged serialized local fixture requests, configurable production-scale starting requests, and institution-neutral SLURM profiles for Apptainer or Singularity.
- Documented high-speed scratch separation, reusable external-resource/container caches, full-WGS storage contributors, container-portability caveats, and deferred production calibration/GPU/retry claims.
- Immutable full Python 3.11 stratification-task image with `ps` support required by native Nextflow resource monitoring; benchmark/comparison scripts and semantics are unchanged.
- Optional pinned Manta 1.6.0 germline and matched tumor-normal SV branches from shared caller-ready BAMs, with final/candidate semantic separation and native SV class, size, PASS, and chromosome summaries.
- Production-only GATK 4.6.2.0 germline CASE gCNV and matched somatic copy-ratio/allelic-segmentation architectures with external model/PoN guards, native outputs, and conservative CNV summaries.
- Optional mosdepth 0.3.14 depth/breadth QC for germline, tumor, and normal BAMs, compact no-per-base defaults, configurable reporting thresholds, standardized coverage tables, and native MultiQC ingestion.
- Optional non-invasive XX/XY/unknown chromosome-configuration provenance, explicit autosome/non-PAR/PAR semantics, and documented deferral of per-region small-variant ploidy routing pending legitimate X/Y data and authoritative PAR resources.
- Performance stage classification for Manta, CNV, mosdepth, and integrated-report tasks without changing performance accounting semantics.
- Optional integrated, terminal WGS analysis reporting with a dependency-free static HTML renderer, embedded publication-oriented styling, inline and independent SVG figures, explicit missing/not-run/NA/valid-empty semantics, HTML escaping, machine-path scrubbing, relative MultiQC/VCF2Circos links, and preserved machine-readable source outputs.
- A permanent real-data HG008 somatic compare-mode report and 1600 × 1200 browser-rendered screenshot covering benchmark, concordance, fixture strata, tumor/normal coverage, valid-empty Manta output, performance, provenance, and graceful CNV-not-run behavior.
- Unit coverage for germline/somatic and single/compare reports, missing optional metrics, NA values, empty SV, CNV not run, malformed summaries, unsafe paths, HTML injection, local links, and deterministic rendering.
- A frozen external GRCh38 production-resource contract covering the Broad/GATK reference and BQSR bundle, GIAB HG002 v5.0q primary and v4.2.1 secondary benchmarks, GIAB/GA4GH v3.6 stratifications, VEP 116, Mutect2, Manta, GATK CNV, and GRCh38 PAR provenance without bundling large resources or universal cohort models.
- Production BWA-MEM2 ALT-sidecar validation and staging when the selected reference contains ALT contigs; compact fixture alignment and all caller parameters remain unchanged.

### Scientific scope

- Added a validated `benchmark_sex=auto|male|female|none` germline hap.py interface. Trusted production HG002 GIAB v5.0q metadata resolve `auto` to `--gender male`; legacy v4.2.1 and fixture defaults retain their prior non-sex-aware command semantics.
- Benchmark provenance records requested/resolved sex, resolution source, the hap.py gender value, and whether sex-aware handling was enabled; mixed v5.0q/v4.2.1 resources and contradictory declarations fail early.
- Full-WGS biological accuracy, compute performance, SV benchmarking, CNV biological validation, and X/Y/PAR caller behavior remain unmeasured.
