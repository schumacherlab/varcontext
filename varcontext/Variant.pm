package Variant;

use strict;
use warnings;

use 5.012;

use Carp;
use Data::Dumper;
 
sub new {
	my $class = shift;	
	my $self = {};
	bless $self, $class;

	my %setters = @_;

	#required fields:chr start ref, alt
	foreach (qw/chr start end ref alt/) {
		$self->{$_} = $setters{$_} if exists $setters{$_};
	}

	#infer type
	if (length($setters{ref}) == length($setters{alt})) {
		$self->{type} = "substitution";
	} elsif (length($setters{ref}) < length($setters{alt})) {
		$self->{type} = "insertion";
	} else {
		$self->{type} = "deletion";
	}

	#calc end if not present
	$self->{end} = $self->{start}  + length($setters{ref}) -1 ;

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
		print STDERR Dumper(\@coords);
		#this is not fatal, but it means that start and end map on different features (gap+coding);
		carp("More than 1 coordinate returned for '" . 
			$self->to_string . "' on " . $tr->stable_id);
		return;
	}

	#if it's a gap, we're in an intron, return undef
	return undef if(ref($coords[0]) eq "Bio::EnsEMBL::Mapper::Gap");

	if($coords[0]->{id} ne "cdna") {
		croak("Mapped coord is not cdna for '" . 
			$self->to_string . "' on " . $tr->stable_id);
	}

	if(not defined $coords[0]->{start}) {
		croak("Start coord not on mapped for '" . 
			$self->to_string . "' on " . $tr->stable_id);
	}

	if(not defined $coords[0]->{end}) {
		croak("End coord not on mapped for '" . 
			$self->to_string . "' on " . $tr->stable_id);
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
		croak "Not yet mapped. Need to call Variant->map_to_Transcript with a Ensmembl Transcript object";
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
	#print STDERR "start=$start, end=$end\n";
	return 0 unless defined $start && defined $end;

	#if the transcript is on the reverse strand we need to reverse the variant
	my $alt = $tr->strand == 1 ? $self->{alt} : $self->reverse_alt;
	my $ref = $tr->strand == 1 ? $self->{ref} : $self->reverse_ref;

	given($self->{type}) {
		when("substitution") {
			$es->edit_substitute($alt, $start, $ref);
		}
		when("insertion") {
			$es->edit_insert($alt, $start);
		}
		when("deletion") {
			$es->edit_delete($ref, $start);
		}
		when("complex") {
			$es->edit_complex($alt, $start, $ref);
		}
	}
	return 1;
	
}

sub to_string {
	my $self = shift;

	return join(" ", map {$self->{$_}} qw/chr start end ref alt/);
}


1;
