use strict;
use warnings;

use constant false => 0;
use constant true  => 1;

use Data::Dumper;
use Text::CSV;
use Getopt::Long;

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/../editseq" );

use Variant;
use VariantSet;
use ensembl;

my $separator = "\t";
my $canonical_flag = false;
my $peptide_flag = false;
my $nmd_flag = false;

GetOptions ("separator=s" => \$separator, "canonical" => \$canonical_flag, "peptide" => \$peptide_flag, "nmd" => \$nmd_flag );

my $vs = VariantSet->new( canonical_only => $canonical_flag, peptide_context => $peptide_flag, nmd_status => $nmd_flag);

my $csv = Text::CSV->new ( { binary => 1, sep_char => $separator } )  # should set binary attribute.
	or die "Cannot use CSV: ".Text::CSV->error_diag ();
my $file = shift;

open my $fh, "<:encoding(utf8)", $file or die "$file: $!";
$csv->column_names ($csv->getline ($fh));
while ( my $row = $csv->getline( $fh ) ) {
	my $id = $row->[2];
	(my $ref = $row->[3]) =~ s/-|\.//g;
	(my $alt = $row->[4]) =~ s/-|\.//g;
	
	# check if multiple alt's present
	if ($alt =~ m/,/) {
		warn "Pruning alt seqs to first listed alt:" . join(",", @$row) . "\n";
		$alt =~ s/,.*//;
	}

	# make new variant and add to variant set
	my $v = Variant->new(variant_id=>$id, chromosome=>$row->[0], start_position=>$row->[1], ref_allele=>$ref, alt_allele=>$alt);
	$vs->add($v);
}

$csv->eof or $csv->error_diag();
close $fh;

$vs->group_variants;
$vs->apply_variants;
$vs->print_variant_context;