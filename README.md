# varcontext

varcontext applies SNVs, indels and complex variants to transcript sequences obtained from ensembl and returns the resulting amino acid sequence, along with other variant information.  
polymorphisms can be taken into account by labeling variants in the mut_id column in the following way: 'gs[0-9]+' or 'rs[0-9]+'.

Input files should adhere to the VCF specification by using the following column order:

1. chromosome
2. start_position
3. mut_id
4. ref_allele
5. alt_allele

## usage

varcontext can be called from a wrapper script or directly from the Terminal by:

`export ENSEMBLAPI=/path/to/ensembl_api/;perl create_context.pl --ARGUMENTS INPUT_FILE 1> OUTPUT_FILE 2> LOG_FILE`

**NOTE: environment variable `ENSEMBLAPI` must be set to full ensembl API path and include trailing slash ('/')**

## arguments

- separator - field separator for input file (default = "\t")
- canonical - only fetch and apply edits to canonical transcript (default: FALSE)
- fullprotein - report full protein sequence (default: TRUE)

## packages

- ensembl - a wrapper around the Ensembl API
- EditSeq
- EditTranscript
- Variant - describes a genomic variant using chromosome start_position ref_allele and alt_allele
- VariantSet - a set of variants and the ability to assign transcripts and edit them

## to do

- VCF input 
- NORMAL control samples 
- Config module (location ensembl api etc) 
- Write tests (maybe a testset of old genes with a couple of designed mutations?)



