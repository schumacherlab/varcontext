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

my $canonical = false;
my $fullpeptide = true;

GetOptions ("canonical" => \$canonical, "fullpeptide" => \$fullpeptide );

my $vs = VariantSet->new( canonical => $canonical, fullpeptide => $fullpeptide);

my $csv = Text::CSV->new ( { binary => 1, sep_char => "\t" } )  # should set binary attribute.
	or die "Cannot use CSV: ".Text::CSV->error_diag ();
my $file = shift;

open my $fh, "<:encoding(utf8)", $file or die "$file: $!";
$csv->column_names ($csv->getline ($fh));
while ( my $row = $csv->getline( $fh ) ) {
	my $id = $row->[2];
	(my $ref = $row->[3]) =~ s/-//g;
	(my $alt = $row->[4]) =~ s/-//g;
	
	# check if variant is SNP
	if ($id =~ m/gs.*/) {
		# warn "Discarding genomic SNP from context generation";
		next;
	}

	# check if multiple alt's present
	if ($alt =~ m/,/) {
		warn "Pruning alt seqs to first listed alt:" . join(",", @$row) . "\n";
		$alt =~ s/,.*//;
	}

	# make new variant and add to variant set
	my $v = Variant->new(id=>$id, chr=>$row->[0], start=>$row->[1], ref=>$ref, alt=>$alt);
	$vs->add($v);
}

$csv->eof or $csv->error_diag();
close $fh;

$vs->group_variants;
$vs->apply_variants;
$vs->print_variant_context;