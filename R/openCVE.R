#' Open Cancer Variant Explorer (CVE) Shiny app
#' @description The openCVE function opens the CVE Shiny app.
#' The function to supplement the R package with the Shiny app was suggested by Dean Attali (http://deanattali.com). Currently, the only extension available is a melanoma co-expression network (WGCNAmelanoma).
#' @param x A dataframe (for single file) or list (for multiple oncotator output files)
#' @param sample_names A character vector with sample name(s)
#' @param extension A character vector of extention name
#' @examples
#' \donttest{
#' openCVE(oncotator_example,"case study")
#' openCVE(oncotator_example,"case study WGCNA","WGCNAmelanoma")
#' }
#' @export
openCVE <- function(x, sample_names, extension=FALSE) {
  #set app directory
  if(extension=="WGCNAmelanoma"){
    appDir <-system.file("Shiny","CVE_WGCNA_melanoma",package="CVE")
  } else {
    appDir <- system.file("Shiny", "CVE", package = "CVE")
  }
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `CVE`.", call. = FALSE)
  }

  #if data frame -> turn into list
  if(class(x)=="data.frame"){
    z = vector("list",1)
    z[[1]] = x
    names(z) = sample_names
    assign("v",z,.GlobalEnv)
  } else if (class(x)=="list"){
    z = x
    names(z) = sample_names
    assign("v",z,.GlobalEnv)
  } else {
    stop("x is neither a data frame nor list. Please change format accordingly.")
  }

  runApp(appDir, display.mode = "normal")
}

####################
#Data set documentation
####################

#' Example Oncotator output for the case study described in the paper
#'
#' A dataset containing the example Oncotator output for the case study described in the paper
#'
#'
#' @docType data
#' @keywords datasets
#' @name oncotator_example
NULL

#' UV signature gene significance (GS)
#'
#' A dataset containing the UV signature gene signficance for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name GS_UV
NULL

#' Lymphocyte score gene significance (GS)
#'
#' A dataset containing the lymphocyte score gene signficance for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name GS_lscore
NULL

#' Primary vs metastasis gene significance (GS)
#'
#' A dataset containing the primary vs metastases gene signficance for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name GS_pmet
NULL

#' Survival gene significance (GS)
#'
#' A dataset containing the survival gene signficance for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name GS_survival
NULL

#' Vemurafenib resistance gene significance (GS)
#'
#' A dataset containing the vemurafenib resistance gene signficance for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name GS_Vem
NULL

#' Gene tree of co-expression network
#'
#' A dataset containing the gene tree of co-expression network for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name METree_GO
NULL

#' Module membership
#'
#' A dataset containing the module membership for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name MM
NULL

#' UV signature module significance (MS)
#'
#' A dataset containing the UV signature module signficance for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name MS_UV
NULL

#' Lymphocyte score module significance (MS)
#'
#' A dataset containing the lymphocyte score module signficance for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name MS_lscore
NULL

#' Primary vs metastasis module significance (MS)
#'
#' A dataset containing the primary vs metastases module signficance for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name MS_pmet
NULL

#' Survival module significance (MS)
#'
#' A dataset containing the survival module signficance for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name MS_survival
NULL

#' Vemurafenib resistance module significance (MS)
#'
#' A dataset containing the vemurafenib resistance module signficance for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name MS_vem
NULL

#' UV signature module significance scaled for barplot
#'
#' A dataset containing the UV signature module signficance scaled for barplot for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name MS_UV_bar
NULL

#' Lymphocyte score module significance scaled for barplot
#'
#' A dataset containing the lymphocyte score module signficance scaled for barplot for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name MS_lscore_bar
NULL

#' Primary vs metastasis module significance scaled for barplot
#'
#' A dataset containing the primary vs metastases module signficance scaled for barplot for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name MS_pmet_bar
NULL

#' Survival module significance scaled for barplot
#'
#' A dataset containing the survival module signficance scaled for barplot for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name MS_survival_bar
NULL

#' Vemurafenib resistance module significance scaled for barplot
#'
#' A dataset containing the vemurafenib resistance module signficance scaled for barplot for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name MS_Vem_bar
NULL

#' Top 5000 most variant genes in TCGA RNAseq data
#'
#' A dataset containing the top 5000 most variant genes in TCGA RNAseq data for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name genes_WGCNA
NULL

#' Module assignment of top 5000 most variant genes in TCGA RNAseq data
#'
#' A dataset containing the module assignment top 5000 most variant genes in TCGA RNAseq data for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name modules
NULL

#' Label order of co-expression modules
#'
#' A dataset containing the label order of co-expression modules for WGCNAmelanoma extension
#'
#'
#' @docType data
#' @keywords datasets
#' @name label_order
NULL

