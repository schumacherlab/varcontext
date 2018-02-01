package EditTranscript;

use strict;
use warnings;

no if ($] >= 5.018), 'warnings' => 'experimental';
use 5.012;

use Data::Dumper;
use Carp;
use Bio::Seq;
use Bio::Tools::CodonTable;

use EditSeq;

# EditTranscript holds a reference to an ensembl transcript and a generic editseq object
# After applying edits to the mRNA sequence in the editseq sequence it has several methods that
# allow getting pre/post-edit context sequence and converted positions.
sub new {
	my $class = shift;	
	my $tr = shift;

	croak("Need transcript object as arg1") unless ref($tr) eq "Bio::EnsEMBL::Transcript";

	my $self = {};
	bless $self, $class;

	#get the transcript translatable RNA sequence
	my $rna = $tr->translateable_seq;
	my $rna_editor = EditSeq->new($rna);

	$self->{transcript} = $tr;
	$self->{rna} = $rna;
	$self->{rna_editor} = $rna_editor;

	return $self;
}


sub apply_variant {
	my $self = shift;
	my $v = shift;

	croak "Edits already applied, cannot add more variants" if $self->{applied_edits};
	croak "Need Variant object as arg1" unless ref($v) eq "Variant";

	#map variant coordinates to transcript
	my ($start, $end) = $v->map_to_Transcript($self->{transcript});
	return 0 unless defined $start && defined $end;
	#if the transcript is on the reverse strand we need to reverse the variant
	my $alt = $self->{transcript}->strand == 1 ? $v->{alt_allele} : $v->reverse_alt;
	my $ref = $self->{transcript}->strand == 1 ? $v->{ref_allele} : $v->reverse_ref;

	given($v->{type}) {
		when("substitution") {
			$self->{rna_editor}->edit_substitute($alt, $start, $ref);
		}
		when("insertion") {
			$self->{rna_editor}->edit_insert($alt, $start);
			$v->{effect} = (length($v->{alt_allele}) % 3) == 0 ? "inframe" : "frameshift";
			$v->{variant_classification} = $v->{type} . "_" . $v->{effect};
		}
		when("deletion") {
			$self->{rna_editor}->edit_delete($ref, $start);
			$v->{effect} = (length($v->{ref_allele}) % 3) == 0 ? "inframe" : "frameshift";
			$v->{variant_classification} = $v->{type} . "_" . $v->{effect};
		}
		when("complex") {
			$self->{rna_editor}->edit_complex($alt, $start, $ref);
			$v->{effect} = abs(length($v->{ref_allele}) - length($v->{alt_allele})) % 3 == 0 ? "inframe" : "frameshift";
			$v->{variant_classification} = $v->{type} . "_" . $v->{effect};
		}
	}

	return 1;
}

sub apply_edits {
	my $self = shift;

	croak "Edit already applied" if $self->{applied_edits};

	$self->{rna_editor}->apply_edits;
	my $edited_rna = $self->{rna_editor}->edited_seq;

	my $tid = $self->{transcript}->stable_id();


	# find the location of the (first) stop codon
	# sometimes the reference sequence also lacks a stopcodon. If this is the case do not extend the
	# tumor peptide sequence
	my $stop = qr/(TAG|TAA|TGA)$/; 
	if ($self->{rna} =~ m/$stop/ && !($edited_rna =~  m/$stop/ && length($edited_rna) % 3 == 0)) {
		#append genomic data until a stop is found in-frame
		my $elongated_edited_rna = $self->_get_genomic_elongation($edited_rna);
		$self->{transcript_extension} = (length($elongated_edited_rna) - length($edited_rna) > 0) ? length($elongated_edited_rna) - length($edited_rna) : 0;

		# if stop lost is combined with an earlier stop gained it is possible that the 
		# extension value becomes < 0
		$edited_rna = $elongated_edited_rna;
	}

	#translate into protein and store edited results
	my $ref_rna_bioseq = Bio::Seq->new(-seq=>$self->{rna}, -id=>$tid . "_ref");
	$self->{ref_protein} = $ref_rna_bioseq->translate->seq;
	my $edit_rna_bioseq = Bio::Seq->new(-seq=>$edited_rna, -id=>$tid . "_edit");
	$self->{edited_protein} = $edit_rna_bioseq->translate->seq;
	$self->{edited_rna} = $edited_rna;
	#in case of multiple stop codons clip after first
	if ($self->{edited_protein} =~ m/\*./) {
		$self->{edited_protein} = substr $self->{edited_protein}, 0, $-[0]+1;
		$self->{edited_rna} = substr $self->{edited_rna}, 0, ($-[0]+1)*3;
	}

	$self->{applied_edits} = 1;
}
	
sub _get_genomic_elongation {
	my $self = shift;
	my $seq = shift; #the current RNA sequence that we are going to elongate

	my $t = $self->{transcript};
	my $seq_len = length $seq;

	my $slice = $t->slice;
	my $extended = 0;

	my $stop = qr/(TAG|TAA|TGA)/;
	while(1) {
		#get the slice from the transcript
		my $start = $t->strand == -1 ? $t->coding_region_start - $extended - 100 : $t->coding_region_end + 1 + $extended;
		my $end =  $t->strand == -1 ? $t->coding_region_start - $extended - 1 : $t->coding_region_end + 100 + $extended;
		#print STDERR "adding 100bp : $start - $end\n";

		#subset slice
		my $subslice = $slice->sub_Slice($start, $end);
		my $eseq = $subslice->seq;
		#reverse complement if on reverse strand
		if( $t->strand == -1) {
			$eseq = reverse $eseq;
			$eseq =~ tr/[CGAT]/[GCTA]/;
		}
		$seq .= $eseq;

		#find a stop in this seq:
		while($seq =~ /$stop/g) {
			return substr($seq, 0, $+[0]) if ($-[0] % 3 == 0);
		}
		$extended += 100;
	}
}

sub convert_position_to_edit {
	my $self = shift;
	my $newpos = $self->{rna_editor}->convert_position_to_edit(shift);
	return $newpos <= length($self->{edited_rna}) ? $newpos : undef;
}

sub get_codon_ref {
	my $self = shift;
	my $pos = int(shift);

	#convert to 0-based
	$pos--;

	croak "Position out of range" unless $pos >= 0  && $pos < length($self->{rna});

	my $codonstart = $pos - ($pos  % 3);
	return substr($self->{rna}, $codonstart, 3);
}

sub get_codon_edit {
	my $self = shift;
	my $pos = int(shift);

	#convert to 0-based
	$pos--;

	croak "Edits not yet applied" unless $self->{applied_edits};
	croak "Position out of range" unless $pos >= 0;
	return "-" unless $pos < length($self->{edited_rna});

	my $codonstart = $pos - ($pos % 3);
	return substr($self->{edited_rna}, $codonstart, 3);
}

sub get_rna_context_ref {
	my $self = shift;
	my $pos = shift;
	my $size = shift || 54;

	#convert to 0-based
	$pos--;

	croak "Position out of range" unless $pos >= 0  && $pos < length($self->{rna});

	my $start = $pos - $size - 1;
	$start = 0 if $start < 0;
	return substr($self->{rna}, $start, $size * 2);
}

sub get_rna_context_edit {
	my $self = shift;
	my $pos = shift;
	my $size = shift || 54;

	#convert to 0-based
	$pos--;

	croak "Edits not yet applied" unless $self->{applied_edits};

	my $start = $pos - $size - 1;
  # make start fall on first position of codon
  $start = $start - ($start % 3);
	$start = 0 if $start < 0;
	return substr($self->{edited_rna}, $start, $size * 2);
}

sub get_protein_context_ref {
	my $self = shift;
	my $pos = shift;
	my $size = shift || 19;

	#convert to 0-based
	$pos--;

	croak "Edits not yet applied" unless $self->{applied_edits};

	my $start = $pos - $size;
	$start = 0 if $start < 0;
	
	my $pepseq = substr($self->{ref_protein}, $start, $size * 2 + 1);

	if (length($pepseq) < (2 * $size + 1) && substr($pepseq, -1) eq "*") {
		while (length($pepseq) < (2 * $size + 1) && $start > 0) {
			$start--;
			$pepseq = substr($self->{ref_protein}, ($start < 0 ? 0 : $start), -1);
		}
	}
	return($pepseq);
}

sub get_protein_context_edit {
	my $self = shift;
	my $pos = shift;
	my $size = shift || 19;

	#convert to 0-based
	$pos--;

	croak "Edits not yet applied" unless $self->{applied_edits};
	croak "Position out of range" unless $pos >= 0;
	return "-" unless $pos < length($self->{edited_protein});

	my $start = $pos - $size;
	$start = 0 if $start < 0;
	
	my $pepseq = substr($self->{edited_protein}, $start, $size * 2 + 1);

	if (length($pepseq) < (2 * $size + 1) && substr($pepseq, -1) eq "*") {
		while (length($pepseq) < (2 * $size + 1) && $start > 0) {
			$start--;
			$pepseq = substr($self->{edited_protein}, ($start < 0 ? 0 : $start), -1);
		}
	}
	return($pepseq);
}

sub nmd_status {
	# - Premature Termination Codon (PTC) is upstream of Normal Termination Codon (NTC)
	# - PTC is between start codon and stop codon (first exon does not necessarily contain start codon, last exon not necessarily contains stop codon)
	# - Rules below flag transcript for NMD
	#1 - PTC is > 100nt downstream of start codon
	#2 - PTC is > 50 nt upstream of the last exon junction complex (deposited on exon-exon junctions)
	#3 - PTC is not in last exon
	my $self = shift;

	croak "Edits not yet applied" unless $self->{applied_edits};

	return ( $self->{nmd_status}, $self->{nmd_remark} ) if exists $self->{nmd_status};
	$self->{nmd_status} = "FALSE";
	$self->{nmd_remark} = "";

	my $stopindex_ref = index $self->{ref_protein}, "*";
	my $stopindex_edit = index $self->{edited_protein}, "*";

	if ($stopindex_edit < $stopindex_ref) {
		my $mapper = $self->{transcript}->get_TranscriptMapper;
		my @coords = $mapper->cds2genomic($stopindex_edit * 3 + 1, $stopindex_edit * 3 + 4);
		my @exons = @{$self->{transcript}->get_all_Exons()};

		if(scalar @coords != 1) {
			# print STDERR Dumper(\@coords);
			#this is not fatal, but it means that start and end map on different features (gap+coding);
			carp "Early stop to genomic coordinate: Error during cds2genomic conversion, start & end map on different features (e.g. exon-intron boundary)";
			return;
		}

		return if(ref($coords[0]) eq "Bio::EnsEMBL::Mapper::Gap");

		my $t_strand = $self->{transcript}->strand;
		my $start_site = $self->{transcript}->coding_region_start;
		my $stop_site = $self->{transcript}->coding_region_end;
		my ($ptc_start, $ptc_end) = ($coords[0]->{start},  $coords[0]->{end});

		#print $ptc_start . "-" . $ptc_end . "\n";

		# sanity check - PTC is between start & stop
		return ( $self->{nmd_status}, $self->{nmd_remark} ) unless ( $ptc_start >= $start_site && $ptc_end <= $stop_site );

		# find out which is the last exon
		my $last_exon = 0;
		foreach my $exon (@exons) {
			$last_exon++;
			if ( $t_strand eq 1 ){
				last if ( $stop_site >= $exon->seq_region_start && $stop_site <= $exon->seq_region_end );	
			} elsif ( $t_strand eq -1 ) {
				last if ( $start_site >= $exon->seq_region_start && $start_site <= $exon->seq_region_end );	
			}
			
		}

		# rule 1 - check if PTC is within 100nt from start codon
		if ( $stopindex_edit * 3 + 1 <= 100 ) {
				$self->{nmd_status} = "FALSE";
				$self->{nmd_remark} = "PTC within 100nt of start";
				return ( $self->{nmd_status}, $self->{nmd_remark} );
		}

		my $exon_count = 0;
		foreach my $exon (@exons) {
			$exon_count++;

			if ( $t_strand eq 1 ) {
				# rule 2 - check if PTC falls within 50nt upstream of last EJC
				if ( $ptc_start >= $exon->seq_region_end - 50 && $ptc_end <= $exon->seq_region_end && $exon_count eq $last_exon - 1 ) {
					$self->{nmd_status} = "FALSE";
					$self->{nmd_remark} = "PTC <50nt upstream of last EJC";
					last;
				} elsif ( $ptc_start >= $exon->seq_region_start && $ptc_end <= $exon->seq_region_end && $exon_count ne $last_exon) {
					$self->{nmd_status} = "TRUE";
					$self->{nmd_remark} = "PTC >50nt upstream of last EJC";
					last;
				# rule 3 - check if we are in last exon; if so, stop
				} elsif ( $exon_count eq $last_exon ) {
					$self->{nmd_remark} = "PTC in last exon";
					last;
				}
			} elsif ( $t_strand eq -1 ) {
				# rule 2 - check if PTC falls within 50nt upstream of last EJC
				if ( $ptc_start >= $exon->seq_region_start && $ptc_end <= $exon->seq_region_start + 50 && $exon_count eq $last_exon - 1 ) {
					$self->{nmd_status} = "FALSE";
					$self->{nmd_remark} = "PTC <50nt upstream of last EJC";
					last;
				} elsif ( $ptc_start >= $exon->seq_region_start && $ptc_end <= $exon->seq_region_end && $exon_count ne $last_exon ) {
					$self->{nmd_status} = "TRUE";
					$self->{nmd_remark} = "PTC >50nt upstream of last EJC";
					last;
				# rule 3 - check if we are in last exon; if so, stop
				} elsif ( $exon_count eq $last_exon ) {
					$self->{nmd_remark} = "PTC in last exon";
					last;
				}
			}
		}
	}
	return ( $self->{nmd_status}, $self->{nmd_remark} );
}


1;
