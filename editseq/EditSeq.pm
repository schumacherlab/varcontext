package EditSeq;

use strict;
use warnings;

no if ($] >= 5.018), 'warnings' => 'experimental';
use 5.012;

use Data::Dumper;
use Carp;

sub edit_insert {
	my $self = shift;
	croak "Sequence has applied edits,  cannot edit anymore" if exists $self->{editedseq};

	my $ins = shift;
	my $pos = shift;
	my $inslen = length $ins;
	
	#edit seq
	#print $self->{seq},"\n";
	#substr $self->{seq}, $newpos, 0, $ins;
	#print $self->{seq},"\n";
	
	#update history
	my $edit =  {type=>'insertion', coord=>$pos, seq=>$ins, len=>$inslen};
	if( $self->insdel_overlap($edit)) {
		# warn "Edit overlaps previous edit DISCARDING!\n";
	} else {
		push @{$self->{_edit_history}}, $edit;
	}
}

sub edit_delete {
	my $self = shift;
	croak "Sequence has applied edits,  cannot edit anymore" if exists $self->{editedseq};

	my $del = shift;
	my $pos = shift;
	my $dellen = length $del;
	
	#edit seq
	#print $self->{seq},"\n";
	#substr $self->{seq}, $newpos, 0, $del;
	#print $self->{seq},"\n";
	
	#update history
	my $edit = {type=>'deletion', coord=>$pos, seq=>$del, len=>$dellen};
	if( $self->insdel_overlap($edit)) {
		# warn "Edit overlaps previous edit DISCARDING!\n";
	} else {
		push @{$self->{_edit_history}}, $edit;
	}
}

sub edit_substitute {
	my $self = shift;
	croak "Sequence has applied edits,  cannot edit anymore" if exists $self->{editedseq};

	my $sub = shift;
	my $pos = shift;
	my $original = shift;

	my $sublen = length $sub;
	
	#update history
	my $edit =  {type=>'substitution', coord=>$pos, seq=>$sub, len=>$sublen};
	$edit->{original} = $original if defined $original;
	push @{$self->{_edit_history}}, $edit;
}

sub edit_complex {
	my $self = shift;
	croak "Sequence has applied edits,  cannot edit anymore" if exists $self->{editedseq};

	my $sub = shift;
	my $pos = shift;
	my $original = shift;

	carp("Provide 3 arguments, substitution, position and original")
		unless defined $sub && defined $pos && defined $original;
	
	return $self->edit_substitute($sub, $pos, $original) if length($original) == length($sub);

	my $sublen = length $sub;
	my $edit =  {type=>'complex', coord=>$pos, seq=>$sub, len=>$sublen, original=>$original, orilen=>length($original)};

	#update history
	if( $self->insdel_overlap($edit)) {
		# warn "Edit overlaps previous edit DISCARDING!\n";
	} else {
		push @{$self->{_edit_history}}, $edit;
	}


}

sub insdel_overlap {
	my $self = shift;
	my $new = shift;

	#no history is fine
	return 0 unless exists $self->{_edit_history};

	#if the edit's an substitition it's ok.
	return 0 if $new->{type} eq "substitution";

	#the coord we need for both ins and del
	my $newmin = $new->{coord};

	#the end coord we only need if it's a del
	my $newmax = $newmin;
	$newmax = $new->{coord} + $new->{len} -1  if $new->{type} eq "deletion";
	$newmax = $new->{coord} + $new->{orilen} -1 if $new->{type} eq "complex";

	#check all previous edits
	foreach my $edit (sort { $b->{coord} <=> $a->{coord} } @{$self->{_edit_history}}) {
		given ($edit->{type}) {
			when("insertion") {
				my $min = $edit->{coord};
				return 1 if ($min > $newmin && $min <= $newmax);
			}
			when("deletion") {
				my $min = $edit->{coord};
				my $max = $edit->{coord} + $edit->{len} - 1;
				return 1 if ($newmin >= $min && $newmin <= $max)  || ($newmax >= $min && $newmax <=$max);
			}
			when("complex") {
				my $min = $edit->{coord};
				my $max = $edit->{coord} + $edit->{orilen} - 1;
				return 1 if ($newmin >= $min && $newmin <= $max)  || ($newmax >= $min && $newmax <=$max);
			}
		}
	}

	return 0;
}

sub apply_edits {
	my $self = shift;
	croak "Sequence has edits already applied." if exists $self->{editedseq};

	return ($self->{seq}, undef) unless $self->hasedits;

	#apply edits
	my $oriseq = $self->{seq};

	
	#create a map of original coordinates to string coordinates
	my %map = map {$_=>$_-1} (1 .. length $oriseq);

	#order all edits based on decreasing coords
	my @insertions = ();
	foreach my $edit (sort { $b->{coord} <=> $a->{coord} } @{$self->{_edit_history}}) {
		given ($edit->{type}) {
			when("insertion") {
				#delay insertions to end
				push @insertions, $edit;
			}
			when("deletion") {
				for my $c ($edit->{coord} .. ($edit->{coord} + $edit->{len} - 1)) {
					delete $map{$c};
				}
			}
			when("complex") {
				#delete the deletion part, put the ins in the queue
				for my $c ($edit->{coord} .. ($edit->{coord} + $edit->{orilen} - 1)) {
					delete $map{$c};
				}
				push @insertions, {type=>'insertion', coord=>$edit->{coord}-1 , seq=>$edit->{seq}, len=>$edit->{len}};
			}
			when("substitution") {
				#substitute using the original coordinates in the original string map doesn't change
				my $replaced = substr $oriseq, $edit->{coord}-1, $edit->{len}, $edit->{seq};
				if (exists $edit->{original} && $replaced ne $edit->{original}) {
					carp "Replaced seq ($replaced) doesn't match supplied original ($edit->{original})\n";
				}
			}
		}
	}

	my $newseq;
	my %inscoord = map { $_->{coord} => $_ } @insertions;
	#print Dumper(\%inscoord, \%map);
	#process the edited seq base by base. If you com to an ins coord ins that sequence first
	my @coordmap;
	foreach my $c (sort { $a <=> $b } keys %map) {
		#print this character
		$newseq .= substr $oriseq, $map{$c}, 1;
		push @coordmap, $map{$c};
		if(exists $inscoord{$c}) {
			$newseq .=  $inscoord{$c}->{seq};
			push @coordmap, (-1)x length( $inscoord{$c}->{seq});
			delete $inscoord{$c};
		}
	}

	#any inserts left?
	croak("Didn't process all insertions? overlaps?\n") if scalar(keys(%inscoord)) > 0;

	#the coordinate map now maps the new sequence to the original positions.
	#store 2 lookup table for coordinate conversion functions
	my %mapedit2original = map { $_ => $coordmap[$_]} 0 .. $#coordmap; 
	my %maporiginal2edit = map { $coordmap[$_] => $_ } 0 .. $#coordmap; 

	#print Dumper("edit2ori", \%mapedit2original, "ori2edit", \%maporiginal2edit);
	

	$self->{editedseq} = $newseq;
	$self->{mapedit2original} = \%mapedit2original;
	$self->{maporiginal2edit} = \%maporiginal2edit;

	return $newseq;
}

sub convert_position_to_original {
	my $self = shift;
	croak "Sequence edits not yet applied, call editseq->apply_edits first" unless exists $self->{editedseq};

	#go from 1 based to 0 based string coordinates
	my $coord = shift(@_) - 1;
	
	croak "No coordinate provided" unless defined $coord;
	#coordinates should be in range of edited sequence
	croak "Coordinate out of range" if $coord < 0 || $coord > length($self->{editedseq});

	my $oricoord = $self->{mapedit2original}->{$coord};
	
	return $oricoord == -1 ? -1 : $oricoord + 1;
}

#this sub will deliver the lowest coordinate in the edited string that is still available
sub convert_position_to_edit {
	my $self = shift;
	croak "Sequence edits not yet applied, call editseq->apply_edits first" unless exists $self->{editedseq};

	#go from 1 based to 0 based string coordinates
	my $coord = shift(@_) - 1;
	
	croak "No coordinate provided" unless defined $coord;
	#coordinates should be in range of original sequence
	croak "Coordinate out of range (pos=$coord)" if $coord < 0 || $coord > length($self->{seq});

	my $newcoord = 	$self->{maporiginal2edit}->{$coord};
	return $newcoord + 1 if defined $newcoord && $newcoord >=0;

	#if not found step left until coord=0
	while($coord != 0) {
		$newcoord = $self->{maporiginal2edit}->{--$coord};
		return $newcoord + 1 if defined $newcoord && $newcoord >=0;
	}
	
	return 1;
}

sub edited_seq {
	my $self = shift;
	croak "Sequence edits not yet applied, call editseq->apply_edits first" unless exists $self->{editedseq};

	return $self->{editedseq};
}

sub substring_ori {
	my $self = shift;

	#go from 1 based to 0 based string coordinates
	my $start = shift(@_) -1;
	my $length = shift;
	
	return substr $self->{seq}, $start, $length;
}

sub substring_edited {
	my $self = shift;
	croak "Sequence edits not yet applied, call editseq->apply_edits first" unless exists $self->{editedseq};

	my $oristart = shift;;
	my $length = shift;

	return substr $self->{editedseq}, $self->convert_position_to_edit($oristart) - 1, $length;
}



sub hasedits {
	my $self = shift;
	return 1 if exists $self->{_edit_history};
	return 0;
}

sub new {
	my $class = shift;
	my $seq = shift;

	my $self = {};
	bless $self, $class;

	$self->{seq} = $seq;
	return $self;
}

1;

