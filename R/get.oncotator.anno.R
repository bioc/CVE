#' @import ConsensusClusterPlus RColorBrewer gplots plyr
#' @import ggplot2 jsonlite ape WGCNA
NULL

#' Open Cancer Variant Explorer (CVE) Shiny app
#'
#' @description The get.oncotator.anno retrieves annotation from the Oncotator database.
#' @param x A matrix containing the columns chromosome, start, end, reference_allele and observed_allele.
#' @examples
#' exCase <- data.frame(chr = rep(10, 3),
#'     start = c("100894110", "100985376", "101137905"),
#'     end = c("100894110", "100985376", "101137905"),
#'     ref_allele = c("T", "C", "G"),
#'     obs_allele = c("G", "A", "A"))
#' get.oncotator.anno(exCase)
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
