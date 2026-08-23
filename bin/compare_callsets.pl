#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my ($analysis, $sample, $expected_sample, $assembly, $bqsr, $benchmark_enabled, $benchmark_name,
    $benchmark_version, $truth_vcf, $benchmark_regions, $output_prefix);
my (@vcf_specs, @index_specs, @caller_version_specs);

GetOptions(
    'analysis=s'          => \$analysis,
    'sample=s'            => \$sample,
    'expected-sample=s'   => \$expected_sample,
    'assembly=s'          => \$assembly,
    'bqsr=s'              => \$bqsr,
    'benchmark-enabled=s' => \$benchmark_enabled,
    'benchmark-name=s'    => \$benchmark_name,
    'benchmark-version=s' => \$benchmark_version,
    'truth=s'             => \$truth_vcf,
    'regions=s'           => \$benchmark_regions,
    'output-prefix=s'     => \$output_prefix,
    'vcf=s@'              => \@vcf_specs,
    'index=s@'            => \@index_specs,
    'caller-version=s@'   => \@caller_version_specs,
) or fail('invalid command-line arguments');

for ($analysis, $sample, $assembly, $bqsr, $benchmark_enabled, $output_prefix) {
    fail('missing required comparison argument') unless defined $_ && length $_;
}
fail("unsupported analysis type '$analysis'") unless $analysis eq 'germline' || $analysis eq 'somatic';
fail('germline comparison requires expected-sample')
    if $analysis eq 'germline' && (!defined $expected_sample || !length $expected_sample);
fail("benchmark-enabled must be true or false") unless $benchmark_enabled eq 'true' || $benchmark_enabled eq 'false';

my @callers = $analysis eq 'germline'
    ? qw(haplotypecaller deepvariant)
    : qw(mutect2 strelka2 deepsomatic);
my %expected = map { $_ => 1 } @callers;
my (%vcf_path, %index_path, %caller_version);
parse_key_values(\@vcf_specs, \%vcf_path, 'VCF');
parse_key_values(\@index_specs, \%index_path, 'index');
parse_key_values(\@caller_version_specs, \%caller_version, 'caller version');

for my $caller (@callers) {
    fail("missing normalized VCF for caller '$caller'") unless defined $vcf_path{$caller};
    fail("missing normalized VCF index for caller '$caller'") unless defined $index_path{$caller};
    fail("missing or empty normalized VCF for caller '$caller': '$vcf_path{$caller}'")
        unless -r $vcf_path{$caller} && -s $vcf_path{$caller};
    fail("missing or empty normalized VCF index for caller '$caller': '$index_path{$caller}'")
        unless -r $index_path{$caller} && -s $index_path{$caller};
}
for my $caller (keys %vcf_path) {
    fail("incompatible caller '$caller' for $analysis comparison") unless $expected{$caller};
}
fail('comparison received an unexpected number of caller VCFs') unless keys(%vcf_path) == @callers;

my (%records, %contigs, %samples_by_caller);
for my $caller (@callers) {
    my ($caller_records, $caller_contigs, $caller_samples) = read_vcf(
        $vcf_path{$caller},
        $analysis eq 'somatic',
        $analysis eq 'germline' ? $expected_sample : undef,
        $caller
    );
    $records{$caller} = $caller_records;
    $contigs{$caller} = $caller_contigs;
    $samples_by_caller{$caller} = $caller_samples;
}
validate_contigs(\@callers, \%contigs, \%records);

my $benchmark_is_enabled = $benchmark_enabled eq 'true';
my (%truth, %regions);
if ($benchmark_is_enabled) {
    fail('benchmark-aware comparison requires benchmark name') unless defined $benchmark_name && length $benchmark_name;
    fail('benchmark-aware comparison requires normalized truth VCF')
        unless defined $truth_vcf && -r $truth_vcf && -s $truth_vcf;
    fail('benchmark-aware comparison requires benchmark-region BED')
        unless defined $benchmark_regions && -r $benchmark_regions && -s $benchmark_regions;
    my ($truth_records, $truth_contigs) = read_vcf($truth_vcf, 0, undef, 'truth');
    %truth = %{$truth_records};
    validate_truth_contigs($truth_contigs, \%contigs, \@callers);
    %regions = %{read_bed($benchmark_regions)};
}

my %union;
for my $caller (@callers) {
    $union{$_} = 1 for keys %{$records{$caller}};
}
my @keys = sort_keys(keys %union);
my @classes = qw(ALL SNV INDEL OTHER);
my @pair_list;
for my $i (0 .. $#callers - 1) {
    for my $j ($i + 1 .. $#callers) {
        push @pair_list, [$callers[$i], $callers[$j]];
    }
}

my (%category_for, %members_for, %class_for, %truth_for, %region_for);
for my $key (@keys) {
    my @present = grep { exists $records{$_}{$key} } @callers;
    my $pattern = join('+', @present);
    $members_for{$key} = \@present;
    $category_for{$key} = category_name($analysis, \@present);
    my (undef, undef, $ref, $alt) = split /\t/, $key;
    $class_for{$key} = variant_class($ref, $alt);
    if ($benchmark_is_enabled) {
        my ($chrom, $pos, $key_ref) = split /\t/, $key;
        my $inside = overlaps_regions($chrom, $pos, length($key_ref), \%regions);
        $region_for{$key} = $inside ? '1' : '0';
        $truth_for{$key} = $inside ? (exists $truth{$key} ? '1' : '0') : 'NA';
    } else {
        $region_for{$key} = 'NA';
        $truth_for{$key} = 'NA';
    }
}

write_membership();
write_concordance_summary();
write_pairwise_jaccard();
write_agreement_matrix();
write_discordant_variants();
write_genotype_outputs() if $analysis eq 'germline';
write_derivative_vcfs();
write_provenance();

sub fail {
    my ($message) = @_;
    print STDERR "CALLSET_COMPARISON ERROR: $message\n";
    exit 1;
}

sub parse_key_values {
    my ($specs, $destination, $label) = @_;
    for my $spec (@{$specs}) {
        my ($key, $value) = split /=/, $spec, 2;
        fail("malformed $label argument '$spec'; expected caller=value")
            unless defined $key && length $key && defined $value && length $value;
        fail("duplicate $label argument for '$key'") if exists $destination->{$key};
        $destination->{$key} = $value;
    }
}

sub open_text {
    my ($path) = @_;
    if ($path =~ /\.gz$/) {
        my $handle = IO::Uncompress::Gunzip->new($path, MultiStream => 1)
            or fail("cannot decompress '$path': $GunzipError");
        return $handle;
    }
    open my $handle, '<', $path or fail("cannot read '$path': $!");
    return $handle;
}

sub read_vcf {
    my ($path, $pass_only, $expected_sample, $label) = @_;
    my $handle = open_text($path);
    my (%vcf_records, %vcf_contigs);
    my @vcf_samples;
    my ($saw_fileformat, $saw_columns, $records_seen) = (0, 0, 0);
    while (my $line = <$handle>) {
        chomp $line;
        if ($line =~ /^##fileformat=VCF/) { $saw_fileformat = 1; next; }
        if ($line =~ /^##contig=<ID=([^,>]+)(?:,length=([0-9]+))?/) {
            $vcf_contigs{$1} = defined $2 ? $2 : '';
            next;
        }
        next if $line =~ /^##/;
        if ($line =~ /^#CHROM\t/) {
            my @header = split /\t/, $line, -1;
            fail("malformed VCF header in '$path'") if @header < 8;
            @vcf_samples = @header > 9 ? @header[9 .. $#header] : ();
            if (defined $expected_sample) {
                fail("$label VCF must contain exactly one germline sample; observed " . scalar(@vcf_samples))
                    unless @vcf_samples == 1;
                fail("mismatched germline sample: $label VCF contains '$vcf_samples[0]', expected '$expected_sample'")
                    unless $vcf_samples[0] eq $expected_sample;
            }
            $saw_columns = 1;
            next;
        }
        next if $line =~ /^#/;
        fail("VCF record encountered before #CHROM header in '$path'") unless $saw_columns;
        my @field = split /\t/, $line, -1;
        fail("malformed VCF record in '$path': expected at least 8 columns") if @field < 8;
        my ($chrom, $pos, undef, $ref, $alt, undef, $filter, undef, $format) = @field;
        fail("malformed VCF coordinate in '$path'") unless $chrom ne '' && $pos =~ /^[1-9][0-9]*$/;
        fail("malformed VCF allele in '$path'") unless $ref ne '' && $alt ne '' && $alt ne '.';
        fail("normalized VCF contract violation in '$path': multi-allelic ALT remains at $chrom:$pos") if $alt =~ /,/;
        next if $pass_only && $filter ne 'PASS';
        my $key = join("\t", $chrom, $pos, uc($ref), uc($alt));
        my $gt = 'NA';
        if (@vcf_samples) {
            my @formats = split /:/, ($format // ''), -1;
            my %format_index = map { $formats[$_] => $_ } 0 .. $#formats;
            if (exists $format_index{GT} && @field > 9) {
                my @sample_values = split /:/, $field[9], -1;
                $gt = defined $sample_values[$format_index{GT}] ? $sample_values[$format_index{GT}] : 'NA';
            }
        }
        if (exists $vcf_records{$key} && $vcf_records{$key}{gt} ne $gt) {
            fail("duplicate normalized allele with conflicting genotypes in '$path': $chrom:$pos:$ref:$alt");
        }
        $vcf_records{$key} = { gt => $gt, filter => $filter };
        $records_seen++;
    }
    close $handle;
    fail("malformed VCF '$path': missing ##fileformat") unless $saw_fileformat;
    fail("malformed VCF '$path': missing #CHROM header") unless $saw_columns;
    return (\%vcf_records, \%vcf_contigs, \@vcf_samples);
}

sub validate_contigs {
    my ($caller_order, $all_contigs, $all_records) = @_;
    my $base = $caller_order->[0];
    for my $caller (@{$caller_order}[1 .. $#{$caller_order}]) {
        my @base_names = sort keys %{$all_contigs->{$base}};
        my @other_names = sort keys %{$all_contigs->{$caller}};
        if (@base_names && @other_names && join("\t", @base_names) ne join("\t", @other_names)) {
            fail("different reference contig declarations between '$base' and '$caller'");
        }
        for my $contig (@base_names) {
            next unless exists $all_contigs->{$caller}{$contig};
            my $a = $all_contigs->{$base}{$contig};
            my $b = $all_contigs->{$caller}{$contig};
            fail("incompatible reference contig length for '$contig' between '$base' and '$caller'")
                if length($a) && length($b) && $a ne $b;
        }
    }
    for my $caller (@{$caller_order}) {
        for my $key (keys %{$all_records->{$caller}}) {
            my ($chrom) = split /\t/, $key;
            fail("allele contig '$chrom' is absent from $caller VCF contig declarations")
                if keys(%{$all_contigs->{$caller}}) && !exists $all_contigs->{$caller}{$chrom};
        }
    }
}

sub validate_truth_contigs {
    my ($truth_contigs, $all_contigs, $caller_order) = @_;
    my $reference_contigs = $all_contigs->{$caller_order->[0]};
    my $shared = 0;
    for my $contig (keys %{$truth_contigs}) {
        next unless exists $reference_contigs->{$contig};
        $shared = 1;
        my $truth_length = $truth_contigs->{$contig};
        my $query_length = $reference_contigs->{$contig};
        fail("truth/query reference length mismatch for contig '$contig'")
            if length($truth_length) && length($query_length) && $truth_length ne $query_length;
    }
    fail('benchmark truth and compared VCFs have no shared declared contigs')
        if keys(%{$truth_contigs}) && keys(%{$reference_contigs}) && !$shared;
}

sub read_bed {
    my ($path) = @_;
    my $handle = open_text($path);
    my %bed;
    my $rows = 0;
    while (my $line = <$handle>) {
        chomp $line;
        next if $line =~ /^\s*$/ || $line =~ /^#/ || $line =~ /^(?:track|browser)(?:\s|$)/;
        my @field = split /\t/, $line, -1;
        fail("malformed benchmark BED '$path'")
            if @field < 3 || $field[1] !~ /^[0-9]+$/ || $field[2] !~ /^[0-9]+$/ || $field[2] <= $field[1];
        push @{$bed{$field[0]}}, [$field[1] + 0, $field[2] + 0];
        $rows++;
    }
    close $handle;
    fail("benchmark BED '$path' has no regions") unless $rows;
    return \%bed;
}

sub overlaps_regions {
    my ($chrom, $pos, $ref_length, $all_regions) = @_;
    return 0 unless exists $all_regions->{$chrom};
    my $zero_start = $pos - 1;
    my $zero_end = $zero_start + $ref_length;
    for my $interval (@{$all_regions->{$chrom}}) {
        return 1 if $zero_start < $interval->[1] && $zero_end > $interval->[0];
    }
    return 0;
}

sub sort_keys {
    return sort {
        my @a = split /\t/, $a;
        my @b = split /\t/, $b;
        $a[0] cmp $b[0] || $a[1] <=> $b[1] || $a[2] cmp $b[2] || $a[3] cmp $b[3]
    } @_;
}

sub variant_class {
    my ($ref, $alt) = @_;
    return 'OTHER' unless $ref =~ /^[ACGTN]+$/i && $alt =~ /^[ACGTN]+$/i;
    return length($ref) == 1 && length($alt) == 1 ? 'SNV' : 'INDEL';
}

sub category_name {
    my ($kind, $present) = @_;
    my %has = map { $_ => 1 } @{$present};
    if ($kind eq 'germline') {
        return 'shared' if @{$present} == 2;
        return 'haplotypecaller_only' if $has{haplotypecaller};
        return 'deepvariant_only';
    }
    return 'all_three_shared' if @{$present} == 3;
    return join('_', @{$present}) . '_only' if @{$present} == 2;
    return $present->[0] . '_specific';
}

sub class_matches {
    my ($wanted, $observed) = @_;
    return $wanted eq 'ALL' || $wanted eq $observed;
}

sub fmt_fraction {
    my ($numerator, $denominator) = @_;
    return 'NA' unless $denominator;
    return sprintf('%.12g', $numerator / $denominator);
}

sub write_membership {
    open my $out, '>', "$output_prefix.membership.tsv" or fail("cannot write membership table: $!");
    print {$out} join("\t", qw(chrom pos ref alt variant_class), @callers,
        qw(agreement_count agreement_pattern truth_allele_match benchmark_region)), "\n";
    for my $key (@keys) {
        my ($chrom, $pos, $ref, $alt) = split /\t/, $key;
        my @presence = map { exists $records{$_}{$key} ? 1 : 0 } @callers;
        print {$out} join("\t", $chrom, $pos, $ref, $alt, $class_for{$key}, @presence,
            scalar(@{$members_for{$key}}), join('+', @{$members_for{$key}}),
            $truth_for{$key}, $region_for{$key}), "\n";
    }
    close $out;
}

sub comparison_categories {
    return $analysis eq 'germline'
        ? qw(union shared haplotypecaller_only deepvariant_only)
        : qw(union all_three_shared mutect2_strelka2_only mutect2_deepsomatic_only strelka2_deepsomatic_only mutect2_specific strelka2_specific deepsomatic_specific);
}

sub write_concordance_summary {
    open my $out, '>', "$output_prefix.concordance_summary.tsv" or fail("cannot write concordance summary: $!");
    print {$out} join("\t", qw(sample_or_patient analysis_type comparison_group variant_class count union_count fraction_of_union benchmark_region_count truth_match_count truth_match_fraction)), "\n";
    for my $class (@classes) {
        my @class_keys = grep { class_matches($class, $class_for{$_}) } @keys;
        my $union_count = scalar @class_keys;
        for my $category (comparison_categories()) {
            my @selected = $category eq 'union'
                ? @class_keys
                : grep { $category_for{$_} eq $category } @class_keys;
            my @in_region = grep { $region_for{$_} eq '1' } @selected;
            my $truth_matches = scalar grep { $truth_for{$_} eq '1' } @selected;
            print {$out} join("\t", $sample, $analysis, $category, $class,
                scalar(@selected), $union_count, fmt_fraction(scalar(@selected), $union_count),
                $benchmark_is_enabled ? scalar(@in_region) : 'NA',
                $benchmark_is_enabled ? $truth_matches : 'NA',
                $benchmark_is_enabled ? fmt_fraction($truth_matches, scalar(@in_region)) : 'NA'), "\n";
        }
    }
    close $out;
}

sub write_pairwise_jaccard {
    open my $out, '>', "$output_prefix.pairwise_jaccard.tsv" or fail("cannot write Jaccard table: $!");
    print {$out} join("\t", qw(sample_or_patient analysis_type caller_a caller_b variant_class intersection_count union_count jaccard)), "\n";
    for my $pair (@pair_list) {
        my ($a, $b) = @{$pair};
        for my $class (@classes) {
            next if $class eq 'OTHER';
            my ($intersection, $pair_union) = (0, 0);
            for my $key (@keys) {
                next unless class_matches($class, $class_for{$key});
                my $in_a = exists $records{$a}{$key};
                my $in_b = exists $records{$b}{$key};
                $pair_union++ if $in_a || $in_b;
                $intersection++ if $in_a && $in_b;
            }
            print {$out} join("\t", $sample, $analysis, $a, $b, $class,
                $intersection, $pair_union, fmt_fraction($intersection, $pair_union)), "\n";
        }
    }
    close $out;
}

sub write_agreement_matrix {
    open my $out, '>', "$output_prefix.agreement_matrix.tsv" or fail("cannot write agreement matrix: $!");
    print {$out} join("\t", 'variant_class', 'caller', @callers), "\n";
    for my $class (@classes) {
        next if $class eq 'OTHER';
        for my $row (@callers) {
            my @counts;
            for my $column (@callers) {
                my $count = scalar grep {
                    class_matches($class, $class_for{$_}) &&
                    exists $records{$row}{$_} && exists $records{$column}{$_}
                } @keys;
                push @counts, $count;
            }
            print {$out} join("\t", $class, $row, @counts), "\n";
        }
    }
    close $out;
}

sub genotype_state {
    my ($raw) = @_;
    return ('missing', 'NA') unless defined $raw && length $raw && $raw ne 'NA';
    my @alleles = split m{[/|]}, $raw, -1;
    return ('missing', 'NA') unless @alleles == 2 && $alleles[0] =~ /^[0-9]+$/ && $alleles[1] =~ /^[0-9]+$/;
    return ('callable', join('/', sort { $a <=> $b } @alleles));
}

sub write_genotype_outputs {
    open my $details, '>', "$output_prefix.genotype_details.tsv" or fail("cannot write genotype details: $!");
    print {$details} join("\t", qw(sample_or_patient chrom pos ref alt variant_class haplotypecaller_gt deepvariant_gt comparison)), "\n";
    my %genotype_result;
    for my $key (@keys) {
        next unless exists $records{haplotypecaller}{$key} && exists $records{deepvariant}{$key};
        my ($hc_state, $hc_norm) = genotype_state($records{haplotypecaller}{$key}{gt});
        my ($dv_state, $dv_norm) = genotype_state($records{deepvariant}{$key}{gt});
        my $result = $hc_state eq 'missing' || $dv_state eq 'missing'
            ? 'missing_or_nocall'
            : ($hc_norm eq $dv_norm ? 'concordant' : 'discordant');
        $genotype_result{$key} = $result;
        my ($chrom, $pos, $ref, $alt) = split /\t/, $key;
        print {$details} join("\t", $sample, $chrom, $pos, $ref, $alt, $class_for{$key},
            $records{haplotypecaller}{$key}{gt}, $records{deepvariant}{$key}{gt}, $result), "\n";
    }
    close $details;

    open my $summary, '>', "$output_prefix.genotype_concordance.tsv" or fail("cannot write genotype summary: $!");
    print {$summary} join("\t", qw(sample_or_patient variant_class shared_alleles genotype_concordant genotype_discordant missing_or_nocall callable_shared genotype_concordance_rate)), "\n";
    for my $class (@classes) {
        next if $class eq 'OTHER';
        my @shared = grep { exists $genotype_result{$_} && class_matches($class, $class_for{$_}) } @keys;
        my $concordant = scalar grep { $genotype_result{$_} eq 'concordant' } @shared;
        my $discordant = scalar grep { $genotype_result{$_} eq 'discordant' } @shared;
        my $missing = scalar grep { $genotype_result{$_} eq 'missing_or_nocall' } @shared;
        my $callable = $concordant + $discordant;
        print {$summary} join("\t", $sample, $class, scalar(@shared), $concordant,
            $discordant, $missing, $callable, fmt_fraction($concordant, $callable)), "\n";
    }
    close $summary;
    return \%genotype_result;
}

sub write_discordant_variants {
    my %genotype_result;
    if ($analysis eq 'germline') {
        for my $key (@keys) {
            next unless exists $records{haplotypecaller}{$key} && exists $records{deepvariant}{$key};
            my ($hc_state, $hc_norm) = genotype_state($records{haplotypecaller}{$key}{gt});
            my ($dv_state, $dv_norm) = genotype_state($records{deepvariant}{$key}{gt});
            $genotype_result{$key} = $hc_state eq 'missing' || $dv_state eq 'missing'
                ? 'missing_or_nocall'
                : ($hc_norm eq $dv_norm ? 'concordant' : 'discordant');
        }
    }
    open my $out, '>', "$output_prefix.discordant_variants.tsv" or fail("cannot write discordant table: $!");
    print {$out} join("\t", qw(sample_or_patient analysis_type chrom pos ref alt variant_class),
        @callers, qw(agreement_count agreement_pattern truth_allele_match benchmark_region notes)), "\n";
    for my $key (@keys) {
        my $all_present = @{$members_for{$key}} == @callers;
        my $genotype_note = $analysis eq 'germline' && exists $genotype_result{$key}
            ? $genotype_result{$key} : '';
        next if $all_present && ($analysis eq 'somatic' || $genotype_note eq 'concordant');
        my ($chrom, $pos, $ref, $alt) = split /\t/, $key;
        my @presence = map { exists $records{$_}{$key} ? 1 : 0 } @callers;
        my $notes = $analysis eq 'germline' && $all_present ? $genotype_note : 'allele_presence_disagreement';
        print {$out} join("\t", $sample, $analysis, $chrom, $pos, $ref, $alt,
            $class_for{$key}, @presence, scalar(@{$members_for{$key}}),
            join('+', @{$members_for{$key}}), $truth_for{$key}, $region_for{$key}, $notes), "\n";
    }
    close $out;
}

sub write_derivative_vcfs {
    my @groups = $analysis eq 'germline'
        ? qw(shared haplotypecaller_only deepvariant_only)
        : qw(all_three_shared mutect2_specific strelka2_specific deepsomatic_specific);
    for my $group (@groups) {
        open my $out, '>', "$output_prefix.comparison.$group.vcf" or fail("cannot write derivative VCF: $!");
        print {$out} "##fileformat=VCFv4.2\n";
        print {$out} "##source=compare_callsets.pl_1.0\n";
        print {$out} "##referenceAssembly=$assembly\n";
        print {$out} '##INFO=<ID=CALLERS,Number=.,Type=String,Description="Callers containing this normalized allele">', "\n";
        print {$out} '##INFO=<ID=AGREEMENT_PATTERN,Number=1,Type=String,Description="Neutral descriptive caller-membership category">', "\n";
        print {$out} "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        for my $key (@keys) {
            next unless $category_for{$key} eq $group;
            my ($chrom, $pos, $ref, $alt) = split /\t/, $key;
            my $caller_list = join(',', @{$members_for{$key}});
            print {$out} join("\t", $chrom, $pos, '.', $ref, $alt, '.', '.',
                "CALLERS=$caller_list;AGREEMENT_PATTERN=$group"), "\n";
        }
        close $out;
    }
}

sub write_provenance {
    open my $out, '>', "$output_prefix.comparison_provenance.tsv" or fail("cannot write comparison provenance: $!");
    print {$out} join("\t", qw(sample_or_patient analysis_type compared_callers caller_versions normalized_input_vcfs vcf_sample_columns reference_assembly bqsr_enabled benchmark_enabled benchmark_name benchmark_version comparison_engine comparison_engine_version allele_semantics genotype_semantics truth_semantics)), "\n";
    my $versions = join(';', map { $_ . '=' . ($caller_version{$_} // '') } @callers);
    my $inputs = join(';', map { $_ . '=' . $vcf_path{$_} } @callers);
    my $sample_columns = join(';', map { $_ . '=' . join(',', @{$samples_by_caller{$_}}) } @callers);
    my $allele_semantics = $analysis eq 'somatic'
        ? 'unique normalized CHROM/POS/REF/ALT presence; FILTER=PASS only'
        : 'unique normalized CHROM/POS/REF/ALT presence; raw germline callset';
    print {$out} join("\t", $sample, $analysis, join(',', @callers), $versions, $inputs,
        $sample_columns, $assembly, $bqsr, $benchmark_enabled,
        ($benchmark_name // ''), ($benchmark_version // ''), 'compare_callsets.pl+bcftools',
        '1.0;bcftools 1.20', $allele_semantics,
        $analysis eq 'germline' ? 'unordered callable diploid GT on shared alleles' : 'not applied',
        $benchmark_is_enabled ? 'allele membership within benchmark BED; distinct from official benchmark engine metrics' : 'unavailable'), "\n";
    close $out;
}
