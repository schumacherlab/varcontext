# varcontext

varcontext applies SNVs, indels and complex variants transcript sequences obtained from ensembl and prints the resulting amino acid sequence, along with other variant information.  
polymorphisms can be taken into account by labeling variants in the ID column in the following way: 'gs[0-9]+' or 'rs[0-9]+'.

Input file should adhere to the VCF specification by using the following column order:

1. chromosome
2. position
3. id
4. ref
5. alt

## usage

varcontext can be called from a wrapper script or directly from the Terminal by:

`perl create_context.pl INPUT_FILE 1> OUTPUT_FILE 2> LOG_FILE --ARGUMENTS`

**NOTE: environment variable `ENSEMBLAPI` must be set to full ensembl API path**

## arguments

- separator - field separator for input file (default = "\t")
- canonical - only fetch and apply edits to canonical transcript (default: FALSE)
- fullprotein - report full protein sequence (default: TRUE)

## packages

- ensembl - a wrapper around the Ensembl API
- Variant - describes a genomic variant using chr start ref and alt
- VariantSet - a set of variants and the ability to assign transcripts and edit them

## to do

- VCF input 
- NORMAL control samples 
- Config module (location ensembl api etc) 
- Write tests (maybe a testset of old genes with a couple of designed mutations?)



