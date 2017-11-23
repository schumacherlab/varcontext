package VariantSet;

use strict;
use warnings;
use POSIX;

use Carp;

use ensembl;
use EditTranscript;
use Bio::Seq;
use Bio::Tools::CodonTable;

my $codontable = Bio::Tools::CodonTable->new();

sub new {
	my $class = shift;
	my $self = {};
	bless $self, $class;

	my %args = @_;

	$self->{variants} = [];
	$self->{stuffadded} = 1;
	# print join(", ", @{$args{extra_field_names}} ) . "\n";
	# print scalar(@{$args{extra_field_names}}) . "\n";
	# Test length of extra_field_names
	# $self->{extra_field_names} = scalar(@{$args{extra_field_names}}) > 0 ?
	#                                     @{$args{extra_field_names}} : ();
	$self->{extra_field_names} = [ @{$args{extra_field_names}} ];
	# print scalar(@{$self->{extra_field_names}}) . "\n";
	# print 'Extra field names: ' .
	#   join(", ", @{$self->{extra_field_names}} ) . "\n";

	$self->{options}->{ensembl_build} = exists $args{ensembl_build} ?
									$args{ensembl_build} : 0;
	$self->{options}->{assembly} = exists $args{assembly} ?
								$args{assembly} : 0;
	$self->{options}->{canonical_only} = exists $args{canonical_only} ?
									$args{canonical_only} : 0;
	$self->{options}->{peptide_context} = exists $args{peptide_context} ?
									$args{peptide_context} : 0;
	$self->{options}->{protein_context} = exists $args{protein_context} ?
										$args{protein_context} : 0;
	$self->{options}->{nmd_status} = exists $args{nmd_status} ?
									$args{nmd_status} : 0;
	$self->{options}->{cdna_context} = exists $args{cdna_context} ?
									$args{cdna_context} : 0;
	$self->{options}->{cdna_contextsize} = exists $args{cdna_contextsize} ?
										$args{cdna_contextsize} : 54;
	$self->{options}->{pepcontextsize} = ceil($args{cdna_contextsize} / 3);
	$self->{options}->{print_overlap} = exists $args{print_overlap} ?
									$args{print_overlap} : 0;

	# prepare an ensembl connection wrapper
	$self->{ens} = ensembl->new(ensembl_build => $self->{options}->{ensembl_build},
								assembly      => $self->{options}->{assembly});

	return $self;
}

sub add {
	my $self = shift;
	my $toadd = shift;

	push @{$self->{variants}}, $toadd;

	$self->{stuffadded} = 1;
}

sub group_variants {
	my $self = shift;

	if ($self->{stuffadded}) {
		my %transcripts = ();
		my %transcripts_to_variants = ();
		foreach my $v (@{$self->{variants}}) {
			foreach my $t (@{$self->{ens}->get_Transcripts_for_Variant($v)}) {
				# only use translatable sequences
				next unless $t->translation;
				# is the transcript complete?
				(my $attrib) = @{$t->get_all_Attributes("CDS_start_NF")};
				next if $attrib && $attrib->value;
				($attrib) = @{$t->get_all_Attributes("CDS_end_NF")};
				next if $attrib && $attrib->value;

				# canonical only?
				next if $self->{options}->{canonical_only} && !$t->is_canonical;

				# store it
				# refetch the transcript for correct slice attachment. This can probably also be
				# achieved differently
				$transcripts{$t->stable_id} = $self->{ens}->{ta}->fetch_by_stable_id($t->stable_id)
				unless exists $transcripts{$t->stable_id};
				push @{$transcripts_to_variants{$t->stable_id}}, $v;
				# add tid to variant
				$v->add_affected_transcriptid($t->stable_id)
			}
		}

		$self->{transcripts} = \%transcripts;
		$self->{transcripts_to_variants} = \%transcripts_to_variants;
		$self->{stuffadded} = 0;
	}
}

sub apply_variants {
	my $self = shift;

	croak "Run group_variants before applying context. Variant set not clean"
		if $self->{stuffadded};

	# loop over all transcripts
	foreach my $tid (keys %{$self->{transcripts}}) {
		# the ensembl object
		my $et = $self->{transcripts}->{$tid};

		#create an editable transcript
		my $germline_transcript = EditTranscript->new($et);
		my $tumor_transcript = EditTranscript->new($et);

		# loop over all variant on this transcript
		my $nedits = 0;
		foreach my $v (@{$self->{transcripts_to_variants}->{$tid}}) {
			# apply variant as edit
			my $haseffect = $tumor_transcript->apply_variant($v);
			if($haseffect == 0) {
				$v->remove_affected_transcriptid($tid);
				next;
			}
			# apply SNP's only  to germline edit (rs identifiers are annotated by MuTect2, can contain pathogenic variants)
			$germline_transcript->apply_variant($v) if exists $v->{variant_id} && $v->{variant_id} =~ m/gs\d+/;
			$nedits += $haseffect;
		}
		if ($nedits > 0) {
			$germline_transcript->apply_edits;
			$tumor_transcript->apply_edits;
			$self->{edit_transcripts}->{tumor}->{$tid} = $tumor_transcript;
			$self->{edit_transcripts}->{germline}->{$tid} = $germline_transcript;
		} else {
			# this $tid is not affected by any variant, remove it from the lists
			delete $self->{transcripts}->{$tid};
			delete $self->{transcripts_to_variants}->{$tid};
		}
	}
}

sub print_variant_context {
	my $self = shift;

	# print header
	my @columns = qw/
		variant_id
		chromosome
		start_position
		end_position
		variant_strand
		ref_allele
		alt_allele
		gene_id
		transcript_id
		transcript_strand
		hugo_symbol
		variant_classification
		transcript_remark
		transcript_extension/;
	push @columns, qw/nmd_status nmd_remark/ if $self->{options}->{nmd_status};
	push @columns, qw/codon_ref codon_germline codon_tumor aa_ref aa_germline
					aa_tumor aa_pos_ref aa_pos_germline aa_pos_tumor_start aa_pos_tumor_stop/;
	push @columns, qw/cdna_context_ref cdna_context_germline cdna_context_tumor/ if $self->{options}->{cdna_context};
	push @columns, qw/peptide_context_ref peptide_context_germline peptide_context_tumor/ if $self->{options}->{peptide_context};
	push @columns, qw/protein_seq_ref protein_seq_germline protein_seq_tumor/ if $self->{options}->{protein_context};
	push @columns, @{$self->{extra_field_names}} if scalar(@{$self->{extra_field_names}}) > 0;
	print join("\t", @columns) . "\n";

	foreach my $v (@{$self->{variants}}) {
		# check if variant is SNP, if so, don't print a context
		# UPDATE: do print SNP info, so we can do sanety check after predictions:
		# No predicted epitopes should be found where contributing variants only contain SNP(s)
		# next if exists $v->{variant_id} && $v->{variant_id} =~ m/^[gr]s\d+$/;

		foreach my $tid (keys %{$v->{affected_transcriptids}}) {

			my $es_germline = $self->{edit_transcripts}->{germline}->{$tid};
			my $es_tumor = $self->{edit_transcripts}->{tumor}->{$tid};

			my %result = map {$_ => ""} @columns;
			$result{protein_seq_ref} = $es_tumor->{ref_protein};
			$result{protein_seq_germline} = $es_germline->{edited_protein};
			$result{protein_seq_tumor} = $es_tumor->{edited_protein};
			$result{transcript_extension} = $es_tumor->{transcript_extension} || 0;

			$result{transcript_remark} = $result{protein_seq_germline} eq $result{protein_seq_tumor} ? "identical" : "somatic_change";

			#ref rna start and end are 1-based
			my ($ref_rna_start, $ref_rna_end) = $v->map_to_transcriptid($tid);
			my $germline_rna_start = $es_germline->convert_position_to_edit($ref_rna_start);
			my $tumor_rna_start = $es_tumor->convert_position_to_edit($ref_rna_start);

			if (defined $tumor_rna_start) {
				#it is still possible $germline_rna_start is undefined

				# this translates to aa residue nr (1-based coordinates)
				my $ref_pep_start = ceil($ref_rna_start / 3);
				my $germline_pep_start = ceil($germline_rna_start/ 3) if defined $germline_rna_start;
				my $tumor_pep_start = ceil($tumor_rna_start / 3);

				if( $v->{type} eq "substitution") {
					$result{codon_ref} = $es_germline->get_codon_ref($ref_rna_start);
					$result{codon_germline} = defined $germline_rna_start ? $es_germline->get_codon_edit($germline_rna_start) : "-";
					$result{codon_tumor} = $es_tumor->get_codon_edit($tumor_rna_start);

					$result{aa_ref} = $codontable->translate($result{codon_ref});
					$result{aa_germline} = defined $germline_rna_start ? $codontable->translate($result{codon_germline}) : "-";
					$result{aa_tumor} = $codontable->translate($result{codon_tumor});

					$v->{variant_classification} = $result{aa_germline} eq $result{aa_tumor} ? "silent_mutation" : "missense_mutation";
					$v->{effect} = $v->{variant_classification};

					if ($result{aa_germline} eq "*" && $result{aa_tumor} ne "*") {
						$v->{variant_classification} = "stop_lost";
					} elsif ($result{aa_germline} ne "*" && $result{aa_tumor} eq "*") {
						$v->{variant_classification} = "stop_gained";
					}

					# check translated codon to substr in pepseq
					my $tr = substr($result{protein_seq_ref}, $ref_pep_start - 1, 1);
					my $ta = substr($result{protein_seq_tumor}, $tumor_pep_start - 1, 1);
					if ( $result{aa_ref} ne $tr || $result{aa_tumor} ne $ta) {
						croak "Codon mismatch ( $result{aa_ref} / $result{aa_tumor} ) vs. ( $tr / $ta ) for " . $v->to_string();
					}
				}

				# context sequences rna
				$result{cdna_context_ref} = $es_germline->get_rna_context_ref($ref_rna_start, $self->{options}->{cdna_contextsize});
				$result{cdna_context_germline} = defined $germline_rna_start ?
					$es_germline->get_rna_context_edit($germline_rna_start, $self->{options}->{cdna_contextsize}) : "-";
				$result{cdna_context_tumor} = $es_tumor->get_rna_context_edit($tumor_rna_start, $self->{options}->{cdna_contextsize});

				# context sequences peptide
				$result{peptide_context_ref} = $es_germline->get_protein_context_ref($ref_pep_start, $self->{options}->{pepcontextsize});
				$result{peptide_context_germline} =defined $germline_rna_start ?
					$es_germline->get_protein_context_edit($germline_pep_start, $self->{options}->{pepcontextsize}) : "-";
				$result{peptide_context_tumor} = $es_tumor->get_protein_context_edit($tumor_pep_start, $self->{options}->{pepcontextsize});

				$result{aa_pos_ref} = $ref_pep_start;
				$result{aa_pos_germline} = $germline_pep_start || "-";
				$result{aa_pos_tumor_start} = $tumor_pep_start;
				$result{aa_pos_tumor_stop} = $tumor_pep_start;

				# tumor peptides
				# if this variant induced a frame shift or mutated the stop codon clip until the end
				if ($v->{variant_classification} eq "stop_lost" || $v->{variant_classification} eq "stop_gained" || $v->{effect} eq "frameshift") {
					$result{aa_pos_tumor_start} = $tumor_pep_start;
					$result{aa_pos_tumor_stop} = length($result{protein_seq_tumor});

					# where should peptide seq start?
					my $pepseq_start = $tumor_pep_start - $self->{options}->{pepcontextsize} - 2; # -2 to correct for: 1-base $tumor_pep_start and trimming of stop codon
					# take substring from pepseq_start until end
					$result{peptide_context_tumor} = substr($result{protein_seq_tumor}, ($pepseq_start < 0 ? 0 : $pepseq_start), -1);

					# if resulting peptide seq is too short, extend
					while (length($result{peptide_context_tumor}) < (2 * $self->{options}->{pepcontextsize} + 1) && $pepseq_start > 0) {
						$pepseq_start--;
						$result{peptide_context_tumor} = substr($result{protein_seq_tumor}, ($pepseq_start < 0 ? 0 : $pepseq_start), -1);
					}
				} elsif ($v->{variant_classification} eq "insertion_inframe") {
					$result{aa_pos_tumor_start} = $tumor_pep_start + 1;
					$result{aa_pos_tumor_stop} = $tumor_pep_start + ((length($v->{alt_allele}) - length($v->{ref_allele})) / 3);
				}
				$result{transcript_remark} = "variant_after_gained_stop" if $tumor_pep_start > length($result{protein_seq_tumor});

			}

			# add remaining info
			$result{$_} = $v->{$_} // "" foreach qw/variant_id chromosome start_position end_position ref_allele alt_allele type effect variant_classification/;
			# $result{extra_fields} = $v->{extra_fields} // "";
			$result{variant_strand} = int 1;
			$result{transcript_id} = $tid;
			($result{gene_id}, $result{hugo_symbol}, $result{transcript_strand}) = $self->{ens}->transcript_info($tid);

			# get NMD status & remark
			( $result{nmd_status}, $result{nmd_remark} ) = $es_tumor->nmd_status;

			# should we print overlapping bases or trimmed bases?
			if ($self->{options}->{print_overlap}) {
				$result{start_position} = $v->{start_position_input};
				$result{ref_allele} = $v->{ref_allele_input};
				$result{alt_allele} = $v->{alt_allele_input};
			}

			# Print the results to STDOUT
			print join("\t", map {$result{$_} // ""} @columns);
			print join("\t", @{$v->{extra_fields}}) . "\n";
		}
	}
}


1;

