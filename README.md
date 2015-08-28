# varcontext

## Current modules
`editseq/` a perl module that can edit a generic string using coordinates and replacement
values. Write test first, implement code later!

```
cd editseq
perl testedit.t
```


## varcontext
Contains the basic infrastructure. Using a couple of packages
 - ensembl - a wrapper around the Ensembl API
 - Variant - Describes a genomic variant using chr start ref and alt
 - VariantSet - A set of variants and the ability to assign transcripts and edit them


The `.pl` are test programs that sort of work.

*TODO*
VCF input
NORMAL control samples
Config module (location ensembl api etc)
Write tests (maybe a testset of old genes with a couple of designed mutations?)



