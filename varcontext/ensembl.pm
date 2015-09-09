package ensembl;
use strict;
use warnings;

use Carp;
 
use FindBin;
use lib ("$FindBin::Bin", map {"/net/bioinf/ensembl/ens72/". $_} qw(ensembl/modules
	ensembl-variation/modules bioperl-live));
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::TranscriptMapper;
use Bio::EnsEMBL::PredictionTranscript;
use Bio::Seq;

sub new {
	my $class = shift;

	my $self = {};
	bless $self, $class;


	#connect to ensembl database
	my $registry = 'Bio::EnsEMBL::Registry';

	$registry->load_registry_from_db(
		-host =>  'legion',
		-user =>  'ensro',
		-pass =>  'ensro',
		-verbose=>0,
		-mysql_skip_secure_auth=>1
	);

	my $sa = Bio::EnsEMBL::Registry->get_adaptor("human","core","slice");
	my $ta = Bio::EnsEMBL::Registry->get_adaptor("human","core","transcript");

	$self->{sa} = $sa;
	$self->{ta} = $ta;

	return $self;
} 

sub get_Transcripts_for_Variant {
	my $self = shift;
	my $variant = shift;

	croak "Variant must be type variant" unless ref $variant eq "Variant";

	my $slice = $self->{sa}->fetch_by_region('chromosome', $variant->{chr}, $variant->{start}, $variant->{end});

	return [] unless defined $slice;
	return $slice->get_all_Transcripts;
}

sub get_genomic_elongation_for_Transcript {
	my $self = shift;
	my $t = shift;
	my $length = shift || 1000;

	croak "Transcript must be type Bio::EnsEMBL::Transcript" unless ref $t eq "Bio::EnsEMBL::Transcript";

	#get the slice from the transcript
	my $slice = $t->slice;
	my $start = $t->strand == -1 ? $t->coding_region_start - $length : $t->coding_region_end + 1;
	my $end =  $t->strand == -1 ? $t->coding_region_start - 1 : $t->coding_region_end + $length;

	#subset slice
	my $subslice = $slice->sub_Slice($start, $end);
	my $seq = $subslice->seq;
	#reverse complement if on reverse strand
	if( $t->strand == -1) {
		$seq = reverse $seq;
		$seq =~ tr/[CGAT]/[GCTA]/;
	}
	return $seq;
}


1;

