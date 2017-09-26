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

# Command line defaults
# Input file should be in TSV-format
my $separator = "\t";
my $canonical = true;
my $peptide_context = true;
my $protein_context = false;
my $nmd = true;
my $rna_context = false;
my $cdnacontextsize = 54;
GetOptions ("separator=s" => \$separator, 
            "canonical!" => \$canonical,
            "peptide_context!" => \$peptide_context, 
            "protein_context!" => \$protein_context, 
            "nmd!" => \$nmd,
            "rna_context!" => \$rna_context,
            "cdnacontextsize=i" => \$cdnacontextsize); 

my $vs = VariantSet->new(canonical_only => $canonical, 
                         peptide_context => $peptide_context, 
                         protein_context => $protein_context, 
                         nmd_status => $nmd,
                         rna_context => $rna_context,
                         cdnacontextsize => $cdnacontextsize); 

# should set binary attribute
my $csv = Text::CSV->new ( { binary => 1, sep_char => $separator } ) 
  or die "Cannot use CSV: ".Text::CSV->error_diag ();
my $file = shift;

open my $fh, "<:encoding(utf8)", $file or die "$file: $!";
my @cols = @{$csv->getline ($fh)};
$csv->column_names (@cols);
while ( my $row = $csv->getline_hr( $fh ) ) {
  my $chromosome = defined $row->{'chromosome'} ? $row->{'chromosome'} : '';
  my $id = $row->{'variant_id'};
  (my $ref = $row->{'ref_allele'}) =~ s/-|\.//g;
  (my $alt = $row->{'alt_allele'}) =~ s/-|\.//g;

  # Set default values for potentially missing input
  my $dna_ref_read_count = defined $row->{'dna_ref_read_count'} ? 
    $row->{'dna_ref_read_count'} : '';
  my $dna_alt_read_count = defined $row->{'dna_alt_read_count'} ? 
    $row->{'dna_alt_read_count'} : '';
  my $dna_total_read_count = defined $row->{'dna_total_read_count'} ? 
    $row->{'dna_total_read_count'} : '';
  my $dna_vaf = defined $row->{'dna_vaf'} ? $row->{'dna_vaf'} : '';

  my $rna_ref_read_count = defined $row->{'rna_ref_read_count'} ? 
    $row->{'rna_ref_read_count'} : '';
  my $rna_alt_read_count = defined $row->{'rna_alt_read_count'} ? 
    $row->{'rna_alt_read_count'} : '';
  my $rna_total_read_count = defined $row->{'rna_total_read_count'} ? 
    $row->{'rna_total_read_count'} : '';
  my $rna_alt_expression = defined $row->{'rna_alt_expression'} ? 
    $row->{'rna_alt_expression'} : '';
  my $rna_vaf = defined $row->{'rna_vaf'} ? $row->{'rna_vaf'} : '';

  # check if multiple alt's present
  if ($alt =~ m/,/) {
    warn "Pruning alt seqs to first listed alt:" . join(",", @$row) . "\n";
    $alt =~ s/,.*//;
  }
   
  # make new variant and add to variant set
  my $v = Variant->new(variant_id=>$id, chromosome=>$chromosome, 
                       start_position=>$row->{'start_position'}, 
                       ref_allele=>$ref, alt_allele=>$alt,
                       dna_ref_read_count=>$dna_ref_read_count, 
                       dna_alt_read_count=>$dna_alt_read_count, 
                       dna_total_read_count=>$dna_total_read_count, 
                       dna_vaf=>$dna_vaf,
                       rna_ref_read_count=>$rna_ref_read_count, 
                       rna_alt_read_count=>$rna_alt_read_count, 
                       rna_total_read_count=>$rna_total_read_count, 
                       rna_vaf=>$rna_vaf,
                       rna_alt_expression=>$rna_alt_expression);
  $vs->add($v);
}

$csv->eof or $csv->error_diag();
close $fh;

$vs->group_variants;
$vs->apply_variants;
$vs->print_variant_context;
