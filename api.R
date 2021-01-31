library(plumber)
source("./vcf-annotation.R")


#* @filter cors
cors <- function(req,res) {
    res$setHeader("Access-Control-Allow-Origin", "*")
    res$setHeader("Accept", "application/json")
  if (req$REQUEST_METHOD == "OPTIONS") {
    res$setHeader("Access-Control-Allow-Methods","*")
    res$setHeader("Access-Control-Allow-Headers", req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS)
    res$status <- 200
    return(list())
    plumber::forward()
  } else {
    plumber::forward()
  }
}

#' @post /annotate-vcf
function(req){
  multipart <- mime::parse_multipart(req)
  #in_file <- multipart$upload$name
  out_file <- multipart$upload$datapath
  annotateVCF(out_file)
}

