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

	carp "Provide chr start ref alt when constructing variant"
		unless defined $setters{chr} && defined $setters{start} && defined $setters{ref} && defined $setters{alt};

	$self->{$_} = $setters{$_} foreach qw/chr start ref alt/;

	$self->{id} = $setters{id} if exists $setters{id};
	
	#trim left identical bases for ref alt combo
	my $minlength  = length($self->{ref}) < length($self->{alt}) ? length($self->{ref}) : length($self->{alt});
	my $pos = 0;
	while($minlength > 0 && $pos < $minlength) {
		last if substr($self->{ref},$pos,1) ne substr($self->{alt},$pos,1);
		$pos++;
	}
	if ($pos > 0) {
		$self->{ref} = substr $self->{ref}, $pos;
		$self->{alt} = substr $self->{alt}, $pos;
		$self->{start} += $pos;
	}

	if($self->{ref} eq $self->{alt}) {
		croak $self->{id} . " # Not a variant (ref/alt identical): '" . $self->to_string . "'";
	}

	#infer type
	if (length($self->{ref}) == length($self->{alt})) {
		$self->{type} = "substitution";
	} elsif (length($self->{ref}) == 0 && length($self->{alt}) > 0) {
		$self->{type} = "insertion";
	} elsif (length($self->{ref}) > 0 && length($self->{alt}) == 0) {
		$self->{type} = "deletion";
	} else {
		$self->{type} = "complex";
	}

	#calc end
	$self->{end} = $self->{start}  + length($setters{ref}) -1 ;

	#prepare target maps
	$self->{affected_transcriptids} = {};
	$self->{transcript_map} = {};
	
	return $self;
}

sub reverse_alt {
	my $self = shift;

	my $rseq = $self->{alt};
	$rseq =~ tr/[CGAT]/[GCTA]/;
	join("",reverse(split //, $rseq));
}

sub reverse_ref {
	my $self = shift;

	my $rseq = $self->{ref};
	$rseq =~ tr/[CGAT]/[GCTA]/;
	join("",reverse(split //, $rseq));
}

sub map_to_Transcript {
	my $self = shift;
	my $tr = shift;

	if(exists $self->{transcript_map}->{$tr->stable_id}  &&
		exists $self->{transcript_map}->{$tr->stable_id}->{start}) {
		return ($self->{transcript_map}->{$tr->stable_id}->{start},
		        $self->{transcript_map}->{$tr->stable_id}->{end});
	}
	

	my $trmapper = $tr->get_TranscriptMapper;
	my @coords = $trmapper->genomic2cds( $self->{start}, $self->{end}, $tr->strand);
	
	#let do some tests and die on exceptions
	if(scalar @coords != 1) {
		# print STDERR Dumper(\@coords);
		#this is not fatal, but it means that start and end map on different features (gap+coding);
		carp $self->{id} . " # More than 1 coordinate returned for '" . $self->to_string . "'";
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

	$self->{transcript_map}->{$tr->stable_id}->{start} = $coords[0]->{start};
	$self->{transcript_map}->{$tr->stable_id}->{end} = $coords[0]->{end};
	return ($coords[0]->{start},  $coords[0]->{end});
}

sub map_to_transcriptid {
	my $self = shift;
	my $tid = shift;

	if(exists $self->{transcript_map}->{$tid} && 
		exists $self->{transcript_map}->{$tid}->{start}) {
		return ($self->{transcript_map}->{$tid}->{start},
		        $self->{transcript_map}->{$tid}->{end});
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
	my $alt = $tr->strand == 1 ? $self->{alt} : $self->reverse_alt;
	my $ref = $tr->strand == 1 ? $self->{ref} : $self->reverse_ref;

	given($self->{type}) {
		when("substitution") {
			$es->edit_substitute($alt, $start, $ref)
		}
		when("insertion") {
			$es->edit_insert($alt, $start);
			$self->{effect} = (length($self->{alt}) % 3) == 0 ? "inframe" : "frameshift";
			$self->{type_effect} = (length($self->{alt}) % 3) == 0 ? $self->{type} . "_" . "inframe" : $self->{type} . "_" . "frameshift"
		}
		when("deletion") {
			$es->edit_delete($ref, $start);
			$self->{effect} = (length($self->{ref}) % 3) == 0 ? "inframe" : "frameshift";
			$self->{type_effect} = (length($self->{ref}) % 3) == 0 ? $self->{type} . "_" . "inframe" : $self->{type} . "_" . "frameshift"
		}
		when("complex") {
			$es->edit_complex($alt, $start, $ref);
			$self->{effect} = abs(length($self->{ref}) - length($self->{alt})) % 3 != 0 ? "inframe" : "frameshift";
			$self->{type_effect} = abs(length($self->{ref}) - length($self->{alt})) % 3 != 0 ? $self->{type} . "_" . "inframe" : $self->{type} . "_" . "frameshift"
		}
	}

	return 1;
	
}

sub to_string {
	my $self = shift;

	return join(" ", map {$self->{$_} // ""} qw/chr start end ref alt type_effect/);
}


1;


