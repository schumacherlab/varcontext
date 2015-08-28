#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: makecontext.pl
#
#        USAGE: ./makecontext.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: YOUR NAME (), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 24-08-15 14:40:02
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

use Data::Dumper;


use lib("./");

use Variant;
use VariantSet;
use ensembl;

#use Text::CSV;
#my $vs = VariantSet->new();
#
#my $csv = Text::CSV->new ( { binary => 1 } )  # should set binary attribute.
#                     or die "Cannot use CSV: ".Text::CSV->error_diag ();
#my $file = "/net/Analysis/data/ni.v.rooij/melanomapairs/data/exome/tcgacolon/20150820\ colo\ msi\ pos/TCGA-A6-3809.csv";
#
#open my $fh, "<:encoding(utf8)", $file or die "$file: $!";
#$csv->column_names ($csv->getline ($fh));
#while ( my $row = $csv->getline( $fh ) ) {
#	my $v = Variant->new(chr=>$row->[2], start=>$row->[3], end=>$row->[4], ref=>$row->[5] =~ s/-//g, alt=>$row->[6]=~ s/-//g);
#	$vs->add($v);
#}
#$csv->eof or $csv->error_diag();
#close $fh;
#

my $vs = VariantSet->new();
my $v = Variant->new(chr=>17, start=>7577540, end=>7577541, ref=>"G", alt=>"");
$vs->add($v);
my $v2 = Variant->new(chr=>17, start=>7577552, end=>7577554, ref=>"CAT", alt=>"GCA");
$vs->add($v2);

$vs->group_variants;
$vs->apply_variants;
$vs->print_variant_context;



