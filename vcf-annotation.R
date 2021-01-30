library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(VariantAnnotation)

.suppress <- function(expression) {
    suppressMessages(suppressWarnings(expression))
}

.mergeAnnotation <- function(x) {
    paste(unique(x), collapse=",")
}

.matchSeqLevels <- function(vcf, db) {
    if(!any(grepl("chr",seqlevels(vcf)))) 
        seqlevels(vcf) <- paste0("chr", seqlevels(vcf))
    vcf <- keepSeqlevels(vcf, intersect(seqlevels(vcf), seqlevels(db)))
    return(vcf)
}

.unpackAnnotations <- function(anno) {
    anno <- anno[[1]][[1]]
    return(paste(anno[["Gene"]],
                anno[["SYMBOL"]],
                anno[["Consequence"]],
                anno[["Existing_variation"]],
                anno[["PolyPhen"]],
                anno[["SIFT"]], sep="|"))
}

checkPloidy <- function(vcf) {
    ploidy <- max(length(strsplit(geno(vcf)$GT[1,],"/")[[1]]),
    return(ploidy)
}

getVariantLocation <- function(vcf, db) {
    vcf <- .matchSeqLevels(vcf, db)
    loc <- suppressMessages(locateVariants(vcf, txdb, AllVariants()))
    loc <- tapply(loc$LOCATION, names(loc), .mergeAnnotation)
    return(data.frame(loc))
}

getVariantAnnotations <- function(variants) {
    variants <- toJSON(variants, auto_unbox=TRUE)
    r <- POST(url = "http://exac.hms.harvard.edu/rest/bulk/variant",
            body = variants,
            add_headers ("Content-Type" = "application/json")
    )
    return(content(r))
}

mungeVCF <- function(vcf) {
    #Alternate allele observation count
    AO <- apply(data.frame(geno(vcf)$AO), 2, unlist)
    #Alternate alleles
    ALT <- unlist(lapply(alt(vcf), as.character))
    ALT_AO <- data.frame(AO, ALT)
    vcf_info <- info(vcf)
    vcf_info$REF <- ref(vcf)
    vcf_info <- data.frame(vcf_info, rowRanges(vcf)[,"QUAL"])
    munged_vcf <- merge(munged_vcf, ALT_AO, by=0)
    return(munged_vcf)
}

annotateVCF <- function(vcf_file, build="hg19", db) {
    vcf <- VariantAnnotation::readVcf(vcf_file, build)
    variants_location <- .suppress(getVariantLocation(vcf, txdb))
    #Get relevant information and splitt multiallelic variants
    munged_vcf <- mungeVCF(vcf)
    annotated_vcf <- merge(munged_vcf, variants_location, by.x=1, by.y=0, all.x=T)
    #Get annotation for exonic variants from ExAC
    annotated_vcf$ExAC_id <- with(annotated_vcf,paste(gsub("chr","",seqnames),start,REF,ALT,sep="-"))
    coding_variants_ids <- which(grepl("coding",annotated_vcf$loc))
    coding_variants <- unique(annotated_vcf$ExAC_id[ex_variants_ids])
    exac_annotations <- getVariantAnnotations(coding_variants)

    exac_annotations[["9-136662928-A-G"]]$variant$allele_freq
    sapply(exac_annotations[1:2], function(x) { 
        paste(sapply(x$consequence, .unpackAnnotations ), collapse=",") 
    })
}




# BiocManager::install('org.Hs.eg.db')
# > as.list(org.Hs.egSYMBOL)[["57801"]]
# [1] "HES4"
coding_snps <- unique(names(loc[which(loc$LOCATION == "coding"),]))
mb <- toJSON(gsub("/|_|:","-",coding_snps), auto_unbox=TRUE)
r <- POST(url = "http://exac.hms.harvard.edu/rest/bulk/variant",
         body = mb,
         add_headers ("Content-Type" = "application/json")
)


# unique(names(loc[which(loc$LOCATION == "coding"),]))

#1. Split multiallelic variants
#2. Identify the position of all variants
#3. Create dataframe with ALT/FREQ proportions after splitting multi allelic
#4. Merge with location of all variants that contain whether is coding, etc.
#5. POST request to ExAC to annotate coding variants