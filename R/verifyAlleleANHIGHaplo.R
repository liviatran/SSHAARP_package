#verifyAlleleANHIGHaplo v 2.0.4 25NOV2024
#'Verifies the alleles in entered haplotype are present in IMGTprotalignments
#'
#'Verifies the alleles in entered haplotype are present in IMGTprotalignments.
#'
#'@param haplotype A haplotype where allele names are written in the IPD-IMGT/HLA Database format, and have 1-4 fields. Alleles in haplotypes may be delimited by "-" or "~".
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or pre-bundled mock haplotype dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the mock haplotype dataset bundled with the package.
#'@param alignments A list object of sub-lists of data frames of protein alignments for the HLA and HLA-region genes supported in the ANHIG/IMGTHLA GitHub Repository. Alignments will always be the most recent version IPD-IMGT/HLA Database version.
#'
#'@importFrom BIGDAWG GetField
#'@importFrom purrr flatten
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return TRUE if all alleles in a haplotype are present in the IMGTprotalignment object. Otherwise, a vector containing FALSE and an error message is returned.

verifyAlleleANHIGHaplo<-function(haplotype, filename, alignments){

  alleles <- strsplit(haplotype, "-|~")[[1]]

  verifyAlleleCheck<-sapply(alleles, function(x) NULL)

  #extract results from verifyAlleleANHIG -- make results into a list for each
  #allele
  for(i in 1:length(alleles)){
    #set variant type to haplotype for checkLocusDataset function
    verifyAlleleCheck[[i]]<-list(verifyAlleleANHIG(alleles[[i]], filename, alignments))
  }

  #if there are any FALSEs in verifyAlleleCheck, return error message
  if(any(unlist(verifyAlleleCheck)==FALSE)){

    #unnest list by one structure
    verifyAlleleCheck<-verifyAlleleCheck %>%
      flatten()

    #look at elements in list to see which alleles contain FALSE
    alleleErrors<-names(verifyAlleleCheck)[sapply(1:length(verifyAlleleCheck), function(x) "FALSE" %in% verifyAlleleCheck[[x]])]

    #paste allele names with associated errors. If there is more than one error
    #all errors will be pasted together and returned as one vector
    return(c(FALSE, paste(paste(names(verifyAlleleCheck[c(alleleErrors)]), sapply(verifyAlleleCheck[c(alleleErrors)], "[", 2), sep =" - "), collapse=" ")))
  }

  return(TRUE)
}
