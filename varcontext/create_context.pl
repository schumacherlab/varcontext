use strict;
use warnings;

use Data::Dumper;
use Text::CSV;
use Getopt::Long;

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/../editseq" );

use Variant;
use VariantSet;
use ensembl;

my $canonical = 1;
my $fullpeptide = 1;

GetOptions ("canonical" => \$canonical, "fullpeptide" => \$fullpeptide );

print $canonical . "\n";
print $fullpeptide . "\n";

my $vs = VariantSet->new( canonical => $canonical ? 1 : 0, fullpeptide => $fullpeptide ? 1 : 0 );

my $csv = Text::CSV->new ( { binary => 1, sep_char => "\t" } )  # should set binary attribute.
	or die "Cannot use CSV: ".Text::CSV->error_diag ();
my $file = shift;

open my $fh, "<:encoding(utf8)", $file or die "$file: $!";
$csv->column_names ($csv->getline ($fh));
while ( my $row = $csv->getline( $fh ) ) {
	my $id = $row->[2];
	(my $ref = $row->[3]) =~ s/-//g;
	(my $alt = $row->[4]) =~ s/-//g;
	if ($alt =~ m/,/) {
		warn "Discarding alt for:" . join(",", @$row) . "\n";
		$alt =~ s/,.*//;
	} elsif ($id =~ m/gs.*/) {
		warn "Discarding genomic SNP from context generation";
		next;
	}
	my $v = Variant->new(id=>$id, chr=>$row->[0], start=>$row->[1], ref=>$ref, alt=>$alt);
	$vs->add($v);
}
$csv->eof or $csv->error_diag();
close $fh;

$vs->group_variants;
$vs->apply_variants;
$vs->print_variant_context;