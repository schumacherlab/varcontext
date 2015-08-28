package VariantSet;

use strict;
use warnings;

use Carp;

use ensembl;
use EditSeq;
use Bio::Seq;

my $CDNACONTEXTSIZE = 54;
my $PEPCONTEXTSIZE = $CDNACONTEXTSIZE / 3 + 1;

sub new {
	my $class = shift;	
	my $self = {};
	bless $self, $class;

	$self->{variants} = [];
	$self->{stuffadded} = 1;

	#prepare an ensembl connection wrapper
	$self->{ens} = ensembl->new();

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
				#only use translatable sequences
				next unless $t->translation;
				#is the transcript complete?
				(my $attrib) = @{$t->get_all_Attributes("CDS_start_NF")};
				next if $attrib && $attrib->value;
				($attrib) = @{$t->get_all_Attributes("CDS_end_NF")};
				next if $attrib && $attrib->value;

				#store it
				#refetch the transcript for correct slice attachment. This can probably also be
				#achieved differently
				$transcripts{$t->stable_id} = $self->{ens}->{ta}->fetch_by_stable_id($t->stable_id)
				unless exists $transcripts{$t->stable_id};
				push @{$transcripts_to_variants{$t->stable_id}}, $v;
				#add tid to variant
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

	croak "Run group_variants before applying context. Variant set not cleaidn"
		if $self->{stuffadded};

	#loop over all transcripts
	foreach my $tid (keys %{$self->{transcripts}}) {
		#the ensembl object
		my $et = $self->{transcripts}->{$tid};

		my $pepseq = $et->translate->seq;
		my $cdna = $et->translateable_seq;
		#the cDNA sequence from the transcript that we are going to edit
		my $tumorcdna = EditSeq->new($cdna);
		
		#loop over all variant on this transcript
		my $nedits = 0;
		foreach my $v (@{$self->{transcripts_to_variants}->{$tid}}) {
			#apply variant as edit
			my $haseffect = $v->apply_to_Transcript($et, $tumorcdna);
			if($haseffect == 0) {
				$v->remove_affected_transcriptid($tid);
			}
			$nedits += $haseffect;
		}
		if ($nedits > 0) {
			$tumorcdna->apply_edits;
			$self->{editedtranscripts}->{$tid} = $tumorcdna;
		} else {
			#this $tid is not affected by any variant, remove it from the lists
			delete $self->{transcripts}->{$tid};
			delete $self->{transcripts_to_variants}->{$tid};
		}
	}
}

sub print_variant_context {
	my $self = shift;

	#print header
	print join("\t", qw/chr start end ref alt transcriptid cdna_context_ref cdna_context_alt peptide_context_ref peptide_context_alt codingchanges/), "\n";
	foreach my $v (@{$self->{variants}}) {
		foreach my $tid (keys %{$v->{affected_transcriptids}}) {
			my $es = $self->{editedtranscripts}->{$tid};
			my $refcdna = $es->{seq};
			my $tumorcdna = $es->edited_seq;

			#create peptides
			my $refcdnabioseq = Bio::Seq->new(-seq=>$refcdna, -id=>"${tid}_ref");
			my $refpepseq = $refcdnabioseq->translate->seq;
			my $tumorcdnabioseq = Bio::Seq->new(-seq=>$tumorcdna, -id=>"${tid}_tumor");
			my $tumorpepseq = $tumorcdnabioseq->translate->seq;

			my ($refcdnastart, $refcdnaend) = $v->map_to_transcriptid($tid);
			my $tumorcdnastart = $es->convert_position_to_edit($refcdnastart);

			#this translates to aa residue nr
			my $refpepstart = int(($refcdnastart -1)/ 3)+1;
			my $tumorpepstart = int(($tumorcdnastart -1)/ 3)+1;

			#make the context coords
			#prepare the string coordinates to clip from (substring is zero based)
			my $stringrefstart = $refcdnastart - $CDNACONTEXTSIZE - 1;
			$stringrefstart = 0 if $stringrefstart < 0;
			my $stringtumorstart = $tumorcdnastart-$CDNACONTEXTSIZE -1;
			$stringtumorstart = 0 if $stringtumorstart < 0;

			my $stringrefpepstart = $refpepstart - $PEPCONTEXTSIZE - 1;
			$stringrefpepstart = 0 if $stringrefpepstart < 0;
			my $stringtumorpepstart = $tumorpepstart - $PEPCONTEXTSIZE - 1;
			$stringtumorpepstart = 0 if $stringtumorpepstart < 0;
			
			#FIXME previous version did longer clips on frameshift effects.
		
			#prepare row result
			my @results = map {$v->{$_}} qw/chr start end ref alt/;
			push @results, $tid;
			#add substrings to result
			push @results, substr($refcdna, $stringrefstart, $CDNACONTEXTSIZE*2);
			push @results, substr($tumorcdna, $stringtumorstart, $CDNACONTEXTSIZE*2);
			push @results, substr($refpepseq, $stringrefpepstart, $PEPCONTEXTSIZE*2);
			push @results, substr($tumorpepseq, $stringtumorpepstart, $PEPCONTEXTSIZE*2);

			push @results, $results[-1] eq $results[-2] ? "identical" : "codingchanges";


			#get the context
			print join("\t", @results), "\n";
		}
	}
}


1;

