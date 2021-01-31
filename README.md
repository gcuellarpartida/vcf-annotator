# VCF-annotator

Simple tool to annotate a VCF. This tool will create a csv file with annotations for variants in the VCF, including their type, position (coding, intergenic, etc) as well as annotations from ExAC (consequence, allele frequency, SIFT and polyphen). Sample information is also included.

## Usage example

In R:

```R
source("./vcf-annotation.R")
anno<-annotateVCF("challange.vcf")
#Sat Jan 30 16:02:47 2021: Reading VCF file challange.vcf...
#Sat Jan 30 16:02:48 2021: Read 6977 variants for 2 samples...
#Sat Jan 30 16:02:48 2021: Getting variants location relative to genes....
#Sat Jan 30 16:03:25 2021: Getting annotations for 3262 coding variants from ExAC...
#Sat Jan 30 16:03:40 2021: Generating .csv...
write.table(anno, file="annotated_vcf.tsv", quote=F, row.names=F, col.names=T, sep="\t")
```

Try the live REST API at https://genoma.io:
```bash
curl -v -F upload=@challange.vcf http://tools.genoma.io/annotate-vcf >annotated_vcf.tsv
``` 

Launch your own REST API:
```R
library(plumber)
r<-plumb("api.R")
r$run(port=8888)
```

Test it:
```bash
curl -v -F upload=@challange.vcf http://localhost:8888/annotate-vcf >annotated_vcf.tsv
```

## Output

The resuts .tsv file contains the following tab separated columns:

| Column       | Description                |
|--------------|----------------------------|
| ExAC_id      | variant name in ExAC friendly format (chr-position-ref-alt) | 
| seqnames     | seqname (chromosome)       |
| start        | variants start position    |
| REF          | reference allele           |
| ALT          | alternate allele           |
| sample_AO  | alternate allele observations    |
| sample_RO  | reference allele observations    |
| sample_AOP | percentage of alternate allele observations vs reference  |
| sample_DP  | sequence depth             |
| loc          | location of variant relative to genes |
| exac_allele_freq  | allele frequency reported by ExAC |
| exac_consequences | variant consequences reported by ExAC separated by semilocon in the following format: GeneID\|Symbol\|Consequence\|rsID\|PolyPhen|SIFT\| |


## Dependencies

This has been tested in [R 4.0.3] and requires the following libraries:
* plumber (1.0.0)
* readr (1.4.0)
* jsonlite (1.6.0)
* httr (1.4.1)
* VariantAnnotation (1.36.0)
* TxDb.Hsapiens.UCSC.hg19.knownGene (3.2.0)

In R:

```R
 install.packages("plumber")
 install.packages("jsonlite")
  install.packages("readr")
 install.packages("httr")
 if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
 BiocManager::install("VariantAnnotation")
 BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
 BiocManager::install('snpStats')
```