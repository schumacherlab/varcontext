package ensembl;
use strict;
use warnings;

use Carp;
use File::Spec;
use FindBin;

BEGIN {
	warn "Please set environment variable ENSEMBLAPI to full Ensembl API path\n" && exit 1 unless defined $ENV{'ENSEMBLAPI'};
}

use lib ("$FindBin::Bin", map {File::Spec->catdir($ENV{'ENSEMBLAPI'}, $_)} qw(ensembl/modules ensembl-variation/modules bioperl-live));
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::TranscriptMapper;
use Bio::EnsEMBL::PredictionTranscript;
use Bio::Seq;

sub new {
	my $class = shift;

	my $self = {};
	bless $self, $class;

	my %args = @_;

	$self->{options}->{ensembl_build} = exists $args{ensembl_build} ?
	$args{ensembl_build} : 0;
	$self->{options}->{assembly} = exists $args{assembly} ?
	$args{assembly} : 0;

	# connect to ensembl database
	#my $registry = 'Bio::EnsEMBL::Registry';

	#$registry->load_registry_from_db(
	#	-host =>  'localhost',
	#	-user =>  'varcontext',
	#	-pass =>  'generatetranscripts',
	#	-verbose=>0,
	#	-mysql_skip_secure_auth=>1
	#);

	my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-host    => 'localhost',
	-user    => 'varcontext',
	-pass    => 'generatetranscripts',
	-port    => '3306',
	-species => 'homo_sapiens',
	-dbname  => join('_', 'homo_sapiens_core', $self->{options}->{ensembl_build}, $self->{options}->{assembly}),
	-group   => 'core',
	);

	my $sa = $dba->get_adaptor("slice");
	my $ta = $dba->get_adaptor("transcript");

	$self->{sa} = $sa;
	$self->{ta} = $ta;

	return $self;
}

sub get_Transcripts_for_Variant {
	my $self = shift;
	my $variant = shift;

	croak "Variant must be type variant" unless ref $variant eq "Variant";

	my $slice = $self->{sa}->fetch_by_region('chromosome', $variant->{chromosome}, $variant->{start_position}, $variant->{end_position});

	return [] unless defined $slice;
	return $slice->get_all_Transcripts;
}

sub get_genomic_elongation_for_Transcript {
	my $self = shift;
	my $t = shift;
	my $seq = shift; #the current RNA sequence that we are going to elongate
	my $seq_len = length $seq;

	croak "Transcript must be type Bio::EnsEMBL::Transcript" unless ref $t eq "Bio::EnsEMBL::Transcript";

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
			return substr($seq, 0, $-[0]) if ($-[0] % 3 == 0);
		}
		$extended += 100;
	}
}

sub transcript_info {
	my $self = shift;
	my $tid = shift;

	my $t = $self->{ta}->fetch_by_stable_id($tid);

	return undef unless defined $t;

	my $g = $t->get_Gene;

	return ($g->stable_id, $g->external_name, $t->strand);

}

sub exon_info {
	my $self = shift;
	my $tid = shift;

	my $t = $self->{ta}->fetch_by_stable_id($tid);

	return undef unless defined $t;

	my @exons = @{ $t->get_all_Exons() };

	return (@exons);

}

1;

