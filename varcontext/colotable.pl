#!/usr/bin/env perl 

use strict;
use warnings;

use Text::CSV;

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/../editseq" );

use Variant;
use VariantSet;
use ensembl;

my $vs = VariantSet->new(canonical=>1, fullpeptide=>1);

my $csv = Text::CSV->new ( { binary => 1 } )  # should set binary attribute.
	or die "Cannot use CSV: ".Text::CSV->error_diag ();
my $file = shift;

open my $fh, "<:encoding(utf8)", $file or die "$file: $!";
$csv->column_names ($csv->getline ($fh));
while ( my $row = $csv->getline( $fh ) ) {
	(my $ref = $row->[5]) =~ s/-//g;
	(my $alt = $row->[6]) =~ s/-//g;
	my $v = Variant->new(chr=>$row->[2], start=>$row->[3], end=>$row->[4], ref=>$ref, alt=>$alt);
	$vs->add($v);
}
$csv->eof or $csv->error_diag();
close $fh;

$vs->group_variants;
$vs->apply_variants;
$vs->print_variant_context;



