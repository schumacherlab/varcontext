# varcontext
## Tool for variant effect prediction and transcript generation

### Introduction

Varcontext applies SNVs, indels and complex variants to transcript sequences obtained from ensembl and returns the resulting amino acid sequence, along with other variant information.  

**Germline variants are only taken into account as such, when labeled in the variant_id column by: 'gs[0-9]+'. (e.g. gs123)**

### Installation instructions

To set up a new MySQL ensembl mirror, do the following:

1. `sudo apt-get update`
2. `sudo apt-get install mysql-server`
3. `sudo mysql_secure_installation` and set up MySQL
4. to mirror the ensembl data, follow the instructions here: [Installing the Ensembl Data](http://www.ensembl.org/info/docs/webcode/mirror/install/ensembl-data.html)  
you only need the data in `homo_sapiens_core_**_**`

For the perl script to run, you also might need the following modules:

1. `perl -MCPAN -e 'install DBI'`
2. `sudo apt-get install libdbd-mysql-perl`

### Input file definition

Input files should contain at least the following columns (in no particular order):

| variant_id | chromosome | start_position | ref_allele | alt_allele |
|------------|------------|----------------|------------|------------|

Optional columns are:

| dna\_ref\_read\_count | dna\_alt\_read\_count | dna\_total\_read\_count | dna_vaf | rna\_ref\_read\_count | rna\_alt\_read\_count | rna\_total\_read\_count | rna_vaf | rna\_alt\_expression |
|--------|--------|----------|--------|--------|--------|----------|----------------|------------|

### Usage example

Varcontext can be called from a wrapper script (availabe in `neolution-prep`) or directly from the Terminal by:

`export ENSEMBLAPI=/path/to/ensembl_api/;perl /path/to/varcontext/create_context.pl --ARGUMENTS INPUT_FILE 1> OUTPUT_FILE 2> LOG_FILE`

**NOTE: environment variable `ENSEMBLAPI` must be set to full ensembl API path and include trailing slash ('/')**

### Optional arguments  

- `--separator=VALUE` - field separator for input file (default = "\t")
- `--canonical` - only fetch and apply edits to canonical transcripts (default: FALSE)
- `--peptide_context` - report peptide context (default: TRUE)
- `--nopeptide_context` - to omit peptide context
- `--protein_context` - report entire protein sequences (default: FALSE)
- `--rna_context` - report RNA sequences (default: FALSE)
- `--cdnacontextsize=VALUE` - define amount of basepairs to flank each variant by 
  (default: 54)
- `--nmd` - infer nonsense-mediated decay status (default = FALSE)




### Package definitions

- ensembl - a wrapper around the Ensembl API
- EditSeq
- EditTranscript
- Variant - describes a genomic variant using `chromosome`, `start_position`, `ref_allele` and `alt_allele`
- VariantSet - a set of variants and the ability to assign transcripts and edit them

### To do

- VCF input 
- NORMAL control samples 
- Write tests (e.g. a testset of genes with a couple of designed mutations)



