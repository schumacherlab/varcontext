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
		-verbose=>0
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



1;

