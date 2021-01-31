library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)
library(jsonlite)
library(httr)

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

getVariantLocation <- function(vcf, db) {
    vcf <- .matchSeqLevels(vcf, db)
    loc <- suppressMessages(locateVariants(vcf, db, AllVariants()))
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
    #Splitting info for multiallelic variants into rows
    AO <- data.frame(apply(geno(vcf)$AO, 2, unlist))
    ALT <- unlist(lapply(alt(vcf), as.character))
    VARS <- sub("[0-9]$", "", rownames(AO))
    DP <- geno(vcf)$DP
    type <- unlist(info(vcf)$TYPE)
    #Alternate alleles
    RO <- geno(vcf)$RO
    colnames(RO) <- paste0(colnames(RO),"_RO")
    colnames(AO) <- paste0(colnames(AO),"_AO")
    colnames(DP) <- paste0(colnames(DP),"_DP")
    AO <- cbind(VARS, AO, ALT, type)
    multiallelic_info<-merge(AO,RO, by.x="VARS", by.y=0)
    multiallelic_info<-merge(multiallelic_info,DP, by.x="VARS", by.y=0)
    for(n in colnames(vcf))
        multiallelic_info[,paste0(n,"_AOP")] <- round(multiallelic_info[,paste0(n,"_AO")]/(multiallelic_info[,paste0(n,"_AO")]+multiallelic_info[,paste0(n,"_RO")]),2)
    vcf_info <- info(vcf)
    vcf_info$REF <- sapply(ref(vcf), as.character)
    vcf_info <- data.frame(vcf_info, rowRanges(vcf)[,"QUAL"])
    munged_vcf <- merge(multiallelic_info, vcf_info, by="VARS", by.y=0)
    #munged_vcf <- merge(munged_vcf, read_counts, by="VARS")
    return(munged_vcf)
}

annotateVCF <- function(vcf_file, build="hg19", db=TxDb.Hsapiens.UCSC.hg19.knownGene) {
    message(date(),": Reading VCF file ", vcf_file,"...")
    vcf <- VariantAnnotation::readVcf(vcf_file, build)
    message(date(),": Read ", nrow(vcf), " variants for ", length(colnames(vcf)), " samples...")
    message(date(),": Getting variants location relative to genes....")
    variants_location <- .suppress(getVariantLocation(vcf, db))
    #Get relevant information and splitt multiallelic variants
    munged_vcf <- mungeVCF(vcf)
    annotated_vcf <- merge(munged_vcf, variants_location, by.x=1, by.y=0, all.x=T)
    #Get annotation for exonic variants from ExAC
    annotated_vcf$ExAC_id <- with(annotated_vcf,paste(gsub("chr","",seqnames),start,REF,ALT,sep="-"))
    coding_variants_ids <- which(grepl("coding",annotated_vcf$loc))
    coding_variants <- unique(annotated_vcf$ExAC_id[coding_variants_ids])
    message(date(),": Getting annotations for ",length(coding_variants), " coding variants from ExAC...")
    exac_annotations <- getVariantAnnotations(coding_variants)
    #Unpack relevant annotations from ExAC API
    message(date(),": Generating .csv...")
    exac_allele_freq <- unlist(sapply(exac_annotations, function(x) {
        x$variant$allele_freq
    }))
    exac_consequences <- sapply(exac_annotations, function(x) { 
        paste(sapply(x$consequence, .unpackAnnotations ), collapse=";") 
    })
    annotated_vcf <- merge(annotated_vcf, data.frame(exac_allele_freq), by.x="ExAC_id", by.y=0, all.x=T)
    annotated_vcf <- merge(annotated_vcf, data.frame(exac_consequences), by.x="ExAC_id", by.y=0, all.x=T)
    keep_cols <- c("ExAC_id","seqnames","start","type","REF","ALT",
                    paste0(colnames(vcf),"_AO"),
                    paste0(colnames(vcf),"_RO"),
                    paste0(colnames(vcf),"_AOP"),
                    paste0(colnames(vcf),"_DP"),
                    "loc", "exac_allele_freq", "exac_consequences"
                    )
    return(annotated_vcf[,keep_cols])
}
