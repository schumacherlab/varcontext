package Variant;

use strict;
use warnings;

no if ($] >= 5.018), 'warnings' => 'experimental';
use 5.012;

use Carp;
use Data::Dumper;
 
sub new {
	my $class = shift;	
	my $self = {};
	bless $self, $class;

	my %setters = @_;

	carp "Provide chromosome start_position ref_allele alt_allele when constructing variant"
		unless defined $setters{chromosome} && defined $setters{start_position} && defined $setters{ref_allele} && defined $setters{alt_allele};

	$self->{$_} = $setters{$_} foreach qw/chromosome start_position ref_allele alt_allele dna_ref_read_count dna_alt_read_count dna_vaf rna_expression/;

	$self->{variant_id} = $setters{variant_id} if exists $setters{variant_id};
	$self->{dna_ref_read_count} = $setters{dna_ref_read_count} if exists $setters{dna_ref_read_count};
	$self->{dna_alt_read_count} = $setters{dna_alt_read_count} if exists $setters{dna_alt_read_count};
	$self->{dna_vaf} = $setters{dna_vaf} if exists $setters{dna_vaf};
	$self->{rna_expression} = $setters{rna_expression} if exists $setters{rna_expression};
	
	#trim left identical bases for ref alt combo
	my $minlength  = length($self->{ref_allele}) < length($self->{alt_allele}) ? length($self->{ref_allele}) : length($self->{alt_allele});
	my $pos = 0;
	while($minlength > 0 && $pos < $minlength) {
		last if substr($self->{ref_allele}, $pos, 1) ne substr($self->{alt_allele}, $pos, 1);
		$pos++;
	}
	if ($pos > 0) {
		$self->{ref_allele} = substr $self->{ref_allele}, $pos;
		$self->{alt_allele} = substr $self->{alt_allele}, $pos;
		$self->{start_position} += $pos;
	}

	if($self->{ref_allele} eq $self->{alt_allele}) {
		croak $self->{variant_id} . " # Not a variant (ref_allele/alt_allele identical): '" . $self->to_string . "'";
	}

	#infer type
	if (length($self->{ref_allele}) == length($self->{alt_allele})) {
		$self->{type} = "substitution";
	} elsif (length($self->{ref_allele}) == 0 && length($self->{alt_allele}) > 0) {
		$self->{type} = "insertion";
	} elsif (length($self->{ref_allele}) > 0 && length($self->{alt_allele}) == 0) {
		$self->{type} = "deletion";
	} else {
		$self->{type} = "complex";
	}

	#calc end
	$self->{end_position} = $self->{start_position} + length($self->{ref_allele}) - 1 ;

	#prepare target maps
	$self->{affected_transcriptids} = {};
	$self->{transcript_map} = {};
	
	return $self;
}

sub reverse_alt {
	my $self = shift;

	my $rseq = $self->{alt_allele};
	$rseq =~ tr/[CGAT]/[GCTA]/;
	join("",reverse(split //, $rseq));
}

sub reverse_ref {
	my $self = shift;

	my $rseq = $self->{ref_allele};
	$rseq =~ tr/[CGAT]/[GCTA]/;
	join("",reverse(split //, $rseq));
}

sub map_to_Transcript {
	my $self = shift;
	my $tr = shift;

	if(exists $self->{transcript_map}->{$tr->stable_id}  &&
		exists $self->{transcript_map}->{$tr->stable_id}->{start_position}) {
		return ($self->{transcript_map}->{$tr->stable_id}->{start_position},
		        $self->{transcript_map}->{$tr->stable_id}->{end_position});
	}

	my $trmapper = $tr->get_TranscriptMapper;
	my @coords = $trmapper->genomic2cds( $self->{start_position}, $self->{end_position}, $tr->strand);
	
	#let do some tests and die on exceptions
	if(scalar @coords != 1) {
		# print STDERR Dumper(\@coords);
		#this is not fatal, but it means that start and end map on different features (gap+coding);
		carp $self->{variant_id} . " # Error during genomic2cds conversion, start_position & end_position map on different features (e.g. exon-intron boundary) '" . $self->to_string . "'";
		return;
	}

	#if it's a gap, we're in an intron, return undef
	return undef if(ref($coords[0]) eq "Bio::EnsEMBL::Mapper::Gap");

	if($coords[0]->{id} ne "cdna") {
		croak("Mapped coord is not cdna for '" . $self->to_string . "' on " . $tr->stable_id);
	}

	if(not defined $coords[0]->{start}) {
		croak("Start coord not on mapped for '" . $self->to_string . "' on " . $tr->stable_id);
	}

	if(not defined $coords[0]->{end}) {
		croak("End coord not on mapped for '" . $self->to_string . "' on " . $tr->stable_id);
	}

	$self->{transcript_map}->{$tr->stable_id}->{start_position} = $coords[0]->{start};
	$self->{transcript_map}->{$tr->stable_id}->{end_position} = $coords[0]->{end};
	return ($coords[0]->{start},  $coords[0]->{end});
}

sub map_to_Genome {
	my $self = shift;
	my $tr = shift;
	my $cdsstart = shift;
	my $cdsend = shift;

	my $trmapper = $tr->get_TranscriptMapper;
	my @coords = $trmapper->cds2genomic( $cdsstart, $cdsend );

	if(scalar @coords != 1) {
		# print STDERR Dumper(\@coords);
		#this is not fatal, but it means that start and end map on different features (gap+coding);
		carp $self->{variant_id} . " # Error during cds2genomic conversion, start & end map on different features (e.g. exon-intron boundary) '" . $self->to_string . "'";
		return;
	}

	#if it's a gap, we're in an intron, return undef
	return undef if(ref($coords[0]) eq "Bio::EnsEMBL::Mapper::Gap");

	if(not defined $coords[0]->{start}) {
		croak("Start coord not on mapped for '" . $self->to_string . "' on " . $tr->stable_id);
	}

	if(not defined $coords[0]->{end}) {
		croak("End coord not on mapped for '" . $self->to_string . "' on " . $tr->stable_id);
	}
	return ($coords[0]->{start},  $coords[0]->{end});
}

sub map_to_transcriptid {
	my $self = shift;
	my $tid = shift;

	if(exists $self->{transcript_map}->{$tid} && 
		exists $self->{transcript_map}->{$tid}->{start_position}) {
		return ($self->{transcript_map}->{$tid}->{start_position},
		        $self->{transcript_map}->{$tid}->{end_position});
	} else {
		croak "Not yet mapped. Need to call Variant->map_to_Transcript with a Ensembl Transcript object";
	}
}

sub add_affected_transcriptid {
	my $self = shift;
	my $tr = shift;

	$self->{affected_transcriptids}->{$tr} = 1;
}

sub remove_affected_transcriptid {
	my $self = shift;
	my $tr = shift;

	delete $self->{affected_transcriptids}->{$tr};
}

sub apply_to_Transcript {
	my $self = shift;

	my $tr = shift;
	my $es = shift;

	croak("Need transcript object as arg1") unless ref($tr) eq "Bio::EnsEMBL::Transcript";
	croak("Need editseq object as arg2") unless ref($es) eq "EditSeq";

	#map my coordinates to transcript
	my ($start, $end) = $self->map_to_Transcript($tr);
	# print STDERR "start=$start, end=$end\n";
	return 0 unless defined $start && defined $end;

	#if the transcript is on the reverse strand we need to reverse the variant
	my $alt = $tr->strand == 1 ? $self->{alt_allele} : $self->reverse_alt;
	my $ref = $tr->strand == 1 ? $self->{ref_allele} : $self->reverse_ref;

	given($self->{type}) {
		when("substitution") {
			$es->edit_substitute($alt, $start, $ref)
		}
		when("insertion") {
			$es->edit_insert($alt, $start);
			$self->{effect} = (length($self->{alt_allele}) % 3) == 0 ? "inframe" : "frameshift";
			$self->{variant_classification} = (length($self->{alt_allele}) % 3) == 0 ? $self->{type} . "_" . "inframe" : $self->{type} . "_" . "frameshift"
		}
		when("deletion") {
			$es->edit_delete($ref, $start);
			$self->{effect} = (length($self->{ref_allele}) % 3) == 0 ? "inframe" : "frameshift";
			$self->{variant_classification} = (length($self->{ref_allele}) % 3) == 0 ? $self->{type} . "_" . "inframe" : $self->{type} . "_" . "frameshift"
		}
		when("complex") {
			$es->edit_complex($alt, $start, $ref);
			$self->{effect} = abs(length($self->{ref_allele}) - length($self->{alt_allele})) % 3 != 0 ? "inframe" : "frameshift";
			$self->{variant_classification} = abs(length($self->{ref_allele}) - length($self->{alt_allele})) % 3 != 0 ? $self->{type} . "_" . "inframe" : $self->{type} . "_" . "frameshift"
		}
	}

	return 1;
	
}

sub to_string {
	my $self = shift;

	return join(" ", map {$self->{$_} // ""} qw/chromosome start_position end_position ref_allele alt_allele variant_classification/);
}


1;


