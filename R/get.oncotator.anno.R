#' Open Cancer Variant Explorer (CVE) Shiny app
#' @description The get.oncotator.anno retrieves annotation from the Oncotator database.
#' @param x A matrix containing the columns chromosome, start, end, reference_allele and observed_allele.
#' @examples
#' get.oncotator.anno(crcCase_input)
#' @export
get.oncotator.anno = function(x){
  urls = paste0("http://portals.broadinstitute.org/oncotator/mutation/",
                as.vector(x[,1]),"_",
                as.vector(x[,2]),"_",
                as.vector(x[,3]),"_",
                as.vector(x[,4]),"_",
                as.vector(x[,5]),"/")
  data.frame(t(sapply(1:length(urls),function(x) fromJSON(urls[x]))))
}

