use strict;
use warnings;
use diagnostics;

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
my $ensembl_build = 90;
my $assembly = 38;
my $canonical = false;
my $peptide_context = false;
my $protein_context = true;
my $nmd = true;
my $cdna_context = false;
my $cdna_contextsize = 54;
my $trim_overlapping_bases = false;

GetOptions ("separator=s"            => \$separator,
            "ensembl=s"              => \$ensembl_build,
            "assembly=s"             => \$assembly,
            "canonical!"             => \$canonical,
            "peptide_context!"       => \$peptide_context,
            "protein_context!"       => \$protein_context,
            "nmd!"                   => \$nmd,
            "cdna_context!"          => \$cdna_context,
            "cdna_contextsize=i"     => \$cdna_contextsize,
            "trim_overlapping_bases" => \$trim_overlapping_bases);

# should set binary attribute
my $csv = Text::CSV->new ( { binary   => 1,
                             sep_char => $separator } )
  or die "Cannot use CSV: ".Text::CSV->error_diag ();
my $file = shift;

open my $fh, "<:encoding(utf8)", $file or die "$file: $!";
my @cols = @{$csv->getline ($fh)};
# Ensure column names are strictly in lowercase
@cols = map { lc } @cols;

# Determine names of non-obligatory, extra columns
# Perl does not have high-level set functionality, so we use a hash instead
my @obligatory_cols = ('variant_id', 'chromosome', 'position',
                       'start_position', 'ref_allele', 'alt_allele');
my @extra_cols = (); my %count_cols = (); my $e;
foreach $e (@cols, @obligatory_cols) {
  ## Obligatory cols will get count of 2, optional cols 1
  $count_cols{$e}++;
}
foreach $e (keys %count_cols) {
  if ($count_cols{$e} == 1 && $e !~ m/start_position/ && $e !~ m/position/) {
    push @extra_cols, $e;
  }
}
# @extra_cols = sort @extra_cols;
# print join(", ", @extra_cols) . "\n";

my $vs = VariantSet->new(ensembl_build     => $ensembl_build,
                         assembly          => $assembly,
                         canonical_only    => $canonical,
                         peptide_context   => $peptide_context,
                         protein_context   => $protein_context,
                         nmd_status        => $nmd,
                         extra_field_names => \@extra_cols,
                         cdna_context      => $cdna_context,
                         cdna_contextsize   => $cdna_contextsize);


$csv->column_names (@cols);
while ( my $row = $csv->getline_hr( $fh ) ) {
  my $chromosome = $row->{'chromosome'};
  if (length($chromosome) == 0) {
    warn "# Missing chromosome information for variant, skipping variant\n" .
         join("\t", map($row->{$_}, @cols)) . "\n";
    next;
  }
  # TODO check validity of input, not just existence
  my $position;
  if (defined $row->{'position'}) {
    $position = $row->{'position'};
  } elsif (defined $row->{'start_position'}) {
    $position = $row->{'start_position'};
  } else {
    warn "# Missing position information for variant, skipping variant\n" .
         join("\t", map($row->{$_}, @cols)) . "\n";
    next;
  }
  (my $ref = $row->{'ref_allele'}) =~ s/-|\.//g;
  (my $alt = $row->{'alt_allele'}) =~ s/-|\.//g;
  my $id = $row->{'variant_id'};

  # check if multiple alt's present
  if ($alt =~ m/,/) {
    warn "Pruning alt seqs to first listed alt:" . join("\t", map($row->{$_}, @cols)) . "\n";
    $alt =~ s/,.*//;
  }

  # Constant ordering of extra fields is ensured using this array
  my @extra_fields = map $row->{$_} // '', @extra_cols;
  # print join(", ", @extra_fields ) . "\n";

  # make new variant and add to variant set
  my $v = Variant->new(variant_id              => $id,
                       chromosome              => $chromosome,
                       start_position          => $position,
                       ref_allele              => $ref,
                       alt_allele              => $alt,
                       extra_fields            => @extra_fields,
                       trim_overlapping_bases  => $trim_overlapping_bases);
  $vs->add($v);
}

$csv->eof or $csv->error_diag();
close $fh;

$vs->group_variants;
$vs->apply_variants;
$vs->print_variant_context;
