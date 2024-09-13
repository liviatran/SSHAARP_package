#verifyAlleleANHIG v 2.0.0 3JAN2022
#'Verifies the allele entered is present in IMGTprotalignments
#'
#'Verifies the allele entered is present in IMGTprotalignments.
#'
#'@param allele An allele name written in the IPD-IMGT/HLA Database format.
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the Solberg dataset.
#'
#'@importFrom BIGDAWG GetField
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return TRUE if allele is present in the IMGTprotalignment object. Otherwise, a vector containing FALSE and an error message is returned.
#'
#'@examples
#'#Example of an allele that exists in the IMGTprotalignment object
#'\dontrun{verifyAlleleANHIG("B*07:02:01:01", filename="your_haplotype_dataset")}
#'
#'#Example of an allele that does not exist in the IMGTprotalignment object
#'verifyAlleleANHIG("B*01:01:01", filename=solberg_dataset)

verifyAlleleANHIG<-function(allele, filename){

  locus<-getVariantInfo(allele)[[1]]

  alleleSyntaxValid<-checkAlleleSyntax(allele, filename)

  #check allele syntax
  if(grepl(TRUE, alleleSyntaxValid[[1]]) == FALSE){
    return(c(FALSE, alleleSyntaxValid[[2]]))
  }

  alleles<-gsub("*", "\\*", strsplit(allele, "~")[[1]], fixed = TRUE)

    if(all(grepl(alleles, SSHAARP::IMGTprotalignments[[locus]]$allele_name, perl=TRUE)==FALSE)){
      return(c(FALSE, "The allele you entered is not present in the current release of ANHIG/IMGTHLA alignments."))
    }

  return(TRUE)
}
