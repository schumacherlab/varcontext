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


my $vs = VariantSet->new();
my $v = Variant->new(chr=>17, start=>7577540, end=>7577541, ref=>"G", alt=>"");
$vs->add($v);
my $v2 = Variant->new(chr=>17, start=>7577552, end=>7577554, ref=>"CAT", alt=>"GCA");
$vs->add($v2);

$vs->group_variants;
$vs->apply_variants;
$vs->print_variant_context;


