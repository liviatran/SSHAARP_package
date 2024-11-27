#verifyAlleleANHIG v 2.0.4 25NOV2024
#'Verifies the allele entered is present in the IMGT protein alignments
#'
#'Verifies the allele entered is present in IMGT protein alignments
#'
#'@param allele An allele name written in the IPD-IMGT/HLA Database format.
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the Solberg dataset.
#'@param alignments A list object of sub-lists of data frames of protein alignments for the HLA and HLA-region genes supported in the ANHIG/IMGTHLA GitHub Repository. Alignments will always be the most recent version IPD-IMGT/HLA Database version.
#'
#'@importFrom BIGDAWG GetField
#'@importFrom HLAtools buildAlignments
#'@importFrom utils data
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return TRUE if allele is present in the IMGTprotalignment object. Otherwise, a vector containing FALSE and an error message is returned.

verifyAlleleANHIG<-function(allele, filename, alignments){

  locus<-getVariantInfo(allele)[[1]]

  alleleSyntaxValid<-checkAlleleSyntax(allele, filename)

  #check allele syntax
  if(grepl(TRUE, alleleSyntaxValid[[1]]) == FALSE){
    return(c(FALSE, alleleSyntaxValid[[2]]))
  }

  if(all(grepl(allele, alignments[[locus]]$allele_name, fixed=TRUE)==FALSE)){
      return(c(FALSE, "The allele you entered is not present in the current release of ANHIG/IMGTHLA alignments."))
  }

  return(TRUE)
}
