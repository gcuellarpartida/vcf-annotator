# VCF-annotator

Simple tool to annotate a VCF. This tool will create a csv file with annotations for variants in the VCF, including their type, position (coding, intergenic, etc) as well as annotations from ExAC. Sample information is also include.

## Example



## Dependencies

This has been tested in [R 4.0.3] and requires the following libraries:
* plumber
* VariantAnnotation
* BSgenome.Hsapiens.UCSC.hg19
* TxDb.Hsapiens.UCSC.hg19.knownGene
* org.Hs.eg.db

In R.
```R
 install.packages("plumber")
 if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
 BiocManager::install("VariantAnnotation")
 BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
 BiocManager::install('snpStats')
 BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
 BiocManager::install('org.Hs.eg.db')
```