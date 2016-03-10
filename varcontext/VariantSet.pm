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

	my %args = @_;

	$self->{variants} = [];
	$self->{stuffadded} = 1;

	#prepare an ensembl connection wrapper
	$self->{ens} = ensembl->new();

	$self->{options}->{canonical} = exists $args{canonical} ? $args{canonical} : 0;
	$self->{options}->{fullpeptide} = exists $args{fullpeptide} ? $args{fullpeptide} : 0;

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

				#canonical only?
				next if $self->{options}->{canonical} && !$t->is_canonical;

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

	croak "Run group_variants before applying context. Variant set not clean"
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

# I have some issues with this sub.
# I think it should be split in a per transcript cdna->peptide part
# followed by the subsetting
# This would also avoid double work when multiple variants affect the same transcript.
sub print_variant_context {
	my $self = shift;

	#print header
	# print join("\t", qw/id chr start end ref alt transcriptid geneid externalname type cdna_context_ref cdna_context_alt peptide_pos_ref peptide_context_ref peptide_pos_alt peptide_context_alt remark effect/, $self->{options}->{fullpeptide} ? "peptide_seq_ref\tpeptide_seq_alt" : ""),  "\n";
  	
  	print join("\t", qw/id chr start end ref alt transcriptid geneid externalname type peptide_pos_ref peptide_pos_alt_start peptide_pos_alt_stop remark effect/, $self->{options}->{fullpeptide} ? "peptide_seq_ref\tpeptide_seq_alt" : ""),  "\n";

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
			
			#find the location of the (first) stop codon
			#sometimes the reference sequence also lacks a stopcodon. If this is the case do not extend the
			#tumor peptide sequence
			my $stopindex_ref = index $refpepseq, "*";
			my $stopindex = index $tumorpepseq, "*";
			my $stoplost = 0;
			if($stopindex == -1 && $stopindex_ref != -1) {
				$stoplost = 1;
				#stop is lost, append genomic data
				$tumorcdna = $self->{ens}->get_genomic_elongation_for_Transcript($self->{transcripts}->{$tid}, $tumorcdna);
				$tumorcdnabioseq = Bio::Seq->new(-seq=>$tumorcdna, -id=>"${tid}_tumor");
				$tumorpepseq = $tumorcdnabioseq->translate->seq;
				$stopindex = index $tumorpepseq, "*";
				croak $v->{id} . " # too little bases added, fix code" if $stopindex == -1;
			}
			#if we moved the stop site we need to redo the peptide generation
			if( $stoplost && $stopindex != length($tumorpepseq) -1) {
				$tumorcdna = substr $tumorcdna,0, (($stopindex+1)*3);
				$tumorcdnabioseq = Bio::Seq->new(-seq=>$tumorcdna, -id=>"${tid}_tumor");
				$tumorpepseq = $tumorcdnabioseq->translate->seq;
				$stopindex = index $tumorpepseq, "*";
				croak $v->{id} . " # too little bases added, fix code" if $stopindex == -1;
			}
			my %result;
			$result{peptide_seq_ref} = $refpepseq;
			$result{peptide_seq_alt} = $tumorpepseq;
			

			my ($refcdnastart, $refcdnaend) = $v->map_to_transcriptid($tid);
			my $tumorcdnastart = $es->convert_position_to_edit($refcdnastart);

			#this translates to aa residue nr
			my $refpepstart = int(($refcdnastart -1)/ 3)+1;
			my $tumorpepstart = int(($tumorcdnastart -1)/ 3)+1;

			#make the context coords
			#prepare the string coordinates to clip from (substring is zero based)
			
			#reference cdna
			my $stringrefstart = $refcdnastart - $CDNACONTEXTSIZE - 1;
			$stringrefstart = 0 if $stringrefstart < 0;
			$result{cdna_context_ref} = substr($refcdna, $stringrefstart, $CDNACONTEXTSIZE*2);

			#reference peptide
			my $stringrefpepstart = $refpepstart - $PEPCONTEXTSIZE - 1;
			$stringrefpepstart = 0 if $stringrefpepstart < 0;
			$result{peptide_context_ref} = substr($refpepseq, $stringrefpepstart, $PEPCONTEXTSIZE*2);
			#refpepcoord <- 
			$result{peptide_pos_ref} = $refpepstart;

			#tumor peptide	
			my $stringtumorpepstart = $tumorpepstart - $PEPCONTEXTSIZE - 1;
			$stringtumorpepstart = 0 if $stringtumorpepstart < 0;
			if($tumorpepstart <= ($stopindex+ 1)) {
				#if this variant induced a frame shift or mutated the stop codon clip until stop
				if(($tumorpepstart == length($refpepseq)) || (exists $v->{effect} && $v->{effect} eq 'frameshift')) {
					$result{peptide_context_alt} = substr($tumorpepseq, $stringtumorpepstart);
					$result{peptide_pos_alt_stop} = length($tumorpepseq);
				} else {
					$result{peptide_context_alt} = substr($tumorpepseq, $stringtumorpepstart, $PEPCONTEXTSIZE*2);
					$result{peptide_pos_alt_stop} = $tumorpepstart;
				}
				$result{peptide_pos_alt_start} = $tumorpepstart;
			} else {
				$result{peptide_pos_alt_start} = "-";
				$result{peptide_pos_alt_stop} = "-";
				$result{peptide_context_alt} = "-";
				$result{remark} = "variant after gained stop";
			}
			
			#tumor cdna
			my $stringtumorstart = $tumorcdnastart-$CDNACONTEXTSIZE -1;
			$stringtumorstart = 0 if $stringtumorstart < 0;
			if($tumorpepstart <= ($stopindex +1)) {
				$result{cdna_context_alt} = substr($tumorcdna, $stringtumorstart, $CDNACONTEXTSIZE*2);
			} else {
				$result{cdna_context_alt} = "-";
			}
		
			#add remaining info
			$result{$_} = $v->{$_} // "" foreach qw/id chr start end ref alt type effect/;
			$result{tid} = $tid;
			($result{geneid}, $result{externalname}) = $self->{ens}->transcript_info($tid);

			if(!exists $result{remark}) {
				$result{remark} = $result{peptide_context_ref} eq $result{peptide_context_alt} ? "identical" : "codingchanges";
			}

			# my @printcolumns = qw/id chr start end ref alt tid geneid externalname type cdna_context_ref cdna_context_alt peptide_pos_ref peptide_context_ref peptide_pos_alt peptide_context_alt remark effect/;

		    my @printcolumns = qw/id chr start end ref alt tid geneid externalname type remark peptide_pos_ref peptide_pos_alt_start peptide_pos_alt_stop/;

			push @printcolumns, ("peptide_seq_ref", "peptide_seq_alt") if $self->{options}->{fullpeptide};

			#get the context
			print join("\t", map {$result{$_} // ""} @printcolumns ), "\n";
		}
	}
}


1;

