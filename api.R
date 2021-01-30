# install.packages("plumber")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("VariantAnnotation")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install('snpStats')
# BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
# BiocManager::install('org.Hs.eg.db')
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(plumber)
library(VariantAnnotation)
library(rtracklayer)

checkPloidy <- function(vcf) {
    ploidy <- max(length(strsplit(geno(vcf)$GT[1,],"/")[[1]]),
    length(strsplit(geno(vcf)$GT[1,],"\\|")[[1]]))
    return(ploidy)
}

matchSeqLevels <- function(vcf, db) {
    if(!any(grepl("chr",seqlevels(vcf)))) 
        seqlevels(vcf) <- paste0("chr", seqlevels(vcf))
    vcf <- keepSeqlevels(vcf, intersect(seqlevels(vcf), seqlevels(db)))
    return(vcf)
}

annnotateVariants <- function(vcf, db) {
    vcf <- matchSeqLevels(vcf, db)
    loc <- locateVariants(vcf_clean, txdb, AllVariants())
    coding <- predictCoding(clean_vcf, txdb, seqSource=Hsapiens)
    coding$GENESYMBOL[coding$GENEID %in% names(as.list(org.Hs.egSYMBOL))] <- unlist(as.list(org.Hs.egSYMBOL)[coding$GENEID[coding$GENEID %in% names(as.list(org.Hs.egSYMBOL))]])
}

vcf <- VariantAnnotation::readVcf("Challenge_data_(1).vcf", "hg19")

# BiocManager::install('org.Hs.eg.db')
# > as.list(org.Hs.egSYMBOL)[["57801"]]
# [1] "HES4"
# coding_snps <- unique(names(loc[which(loc$LOCATION == "coding"),]))
# temp<-head(names(loc),500)
# mb <- toJSON(gsub("/|_|:","-",coding_snps), auto_unbox=TRUE)
# r <- POST(url = "http://exac.hms.harvard.edu/rest/bulk/variant",
#         body = mb,
#         add_headers ("Content-Type" = "application/json")
#     )


# unique(names(loc[which(loc$LOCATION == "coding"),]))

#1. Split multiallelic variants
#2. Identify the position of all variants
#3. Create dataframe with ALT/FREQ proportions after splitting multi allelic
#4. Merge with location of all variants that contain whether is coding, etc.
#5. POST request to ExAC to annotate coding variants

mungeVCF <- function(vcf) {
    #Alternate allele observation count
    AO <- apply(data.frame(geno(vcf)$AO), 2, unlist)
    #Alternate alleles
    ALT <- unlist(lapply(alt(vcf), as.character))
    ALT_AO <- data.frame(AO, ALT)
    vcf_info <- info(vcf)
    vcf_info$REF <- ref(vcf)
    vcf_info <- data.frame(vcf_info, rowRanges(vcf)[,"QUAL"])
    vcf_munged <- merge(vcf_info, ALT_AO, by=0)
    return(vcf_munged)
}

summarizeVCF <- function(vcf_munged) {

}

annotateVariants <- function(variants) {

}

annotateVariants <- function()
