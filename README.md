# varcontext

Varcontext applies SNVs, indels and complex variants to transcript sequences obtained from ensembl and returns the resulting amino acid sequence, along with other variant information.  
*Polymorphisms can be taken into account by labeling variants in the mut_id column in the following way: 'gs[0-9]+'.*

Input files should contain at least the following columns (in identical order):

| variant_id | chromosome | start_position | ref_allele | alt_allele |
|------------|------------|----------------|------------|------------|

Optional columns are:

| dna\_ref\_read\_count | dna\_alt\_read\_count | dna\_total\_read\_count | dna_vaf | rna\_ref\_read\_count | rna\_alt\_read\_count | rna\_total\_read\_count | rna_vaf | rna\_alt\_expression |
|--------|--------|----------|--------|--------|--------|----------|----------------|------------|

**NOTE: Columns are read in order, from left to right. No header checking is done (i.e. first line is skipped), so column order is VITAL to proper processing of the input file.**

## install

To set up a new MySQL ensembl mirror, do the following:

1. `sudo apt-get update`
2. `sudo apt-get install mysql-server`
3. `sudo mysql_secure_installation`
4. follow the instructions here: [Installing the Ensembl Data](http://www.ensembl.org/info/docs/webcode/mirror/install/ensembl-data.html). You only need the data in `homo_sapiens_core_**_**`

You also might need the following perl modules:

1. `perl -MCPAN -e 'install DBI'`
2. `apt-get install libdbd-mysql-perl`

## usage

Varcontext can be called from a wrapper script (availabe in `neolution-prep`) or directly from the Terminal by:

`export ENSEMBLAPI=/path/to/ensembl_api/;perl /path/to/varcontext/create_context.pl --ARGUMENTS INPUT_FILE 1> OUTPUT_FILE 2> LOG_FILE`

**NOTE: environment variable `ENSEMBLAPI` must be set to full ensembl API path and include trailing slash ('/')**

## arguments

- separator - field separator for input file (default = "\t")
- canonical - only fetch and apply edits to canonical transcripts (default: FALSE)
- nmd - infer nonsense-mediated decay status (default = FALSE)
- peptide - report peptide context only (default: FALSE)

## packages

- ensembl - a wrapper around the Ensembl API
- EditSeq
- EditTranscript
- Variant - describes a genomic variant using `chromosome`, `start_position`, `ref_allele` and `alt_allele`
- VariantSet - a set of variants and the ability to assign transcripts and edit them

## to do

- VCF input 
- NORMAL control samples 
- Write tests (e.g. a testset of genes with a couple of designed mutations)



