# VCF-annotator

Simple tool to annotate a VCF. This tool will create a csv file with annotations for variants in the VCF, including their type, position (coding, intergenic, etc) as well as annotations from ExAC. Sample information is also include.

## Usage example

In R:

```R
source("./vcf-annotation.R")
anno<-annotateVCF("challange.vcf")
write.table(anno, file="annotated_vcf.tsv", quote=F, row.names=F, col.names=T, sep="\t")
```

Try the live REST API:

Launch your own REST API:

## Output

The resuts .tsv file contains the following tab separated columns:

ExAC_id	seqnames	start	type	REF	ALT	normal_AO	vaf5_AO	normal_RO	vaf5_RO	normal_AOP	vaf5_AOP	normal_DP	vaf5_DP	loc	exac_allele_freq	exac_consequences

| Column       | Description                |
|--------------|----------------------------|
| ExAC_id      | variant name in ExAC friendly format (chr-position-ref-alt) | 
| seqnames     | seqname (chromosome)       |
| start        | variants start position    |
| REF          | reference allele           |
| ALT          | alternate allele           |
| <sample>_AO  | alternate allele observations    |
| <sample>_RO  | reference allele observations    |
| <sample>_AOP | percentage of alternate allele observations vs reference  |
| <sample>_DP  | sequence depth             |
| loc          | location of variant relative to genes |
| exac_allele_freq  | allele frequency reported by ExAC |
| exac_consequences | variant consequences reported by ExAC separated by semilocon in the following format: *GeneID\|Symbol\|Consequence\|rsID\|PolyPhen|SIFT\|* |








## Dependencies

This has been tested in [R 4.0.3] and requires the following libraries:
* plumber
* jsonlite
* httr
* VariantAnnotation
* BSgenome.Hsapiens.UCSC.hg19
* TxDb.Hsapiens.UCSC.hg19.knownGene
* org.Hs.eg.db

In R.
```R
 install.packages("plumber")
 install.packages("jsonlite")
 install.packages("httr")
 if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
 BiocManager::install("VariantAnnotation")
 BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
 BiocManager::install('snpStats')
 BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
 BiocManager::install('org.Hs.eg.db')
```