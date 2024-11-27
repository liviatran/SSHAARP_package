#checkHaplotypeSyntax v2.0.3 20NOV2024
#'Check haplotype syntax
#'
#'Checks if alleles in a haplotype have correct syntax and an appropriate number of fields.
#'
#'@param haplotype A haplotype where allele names are written in the IPD-IMGT/HLA Database format, and have 1-4 fields. Alleles in haplotypes may be delimited by "-" or "~".
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or pre-bundled mock haplotype dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the mock haplotype dataset bundled with the package.
#'
#'@importFrom purrr flatten
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return TRUE if all alleles in entered haplotype have correct syntax and appropriate number of fields. Otherwise, a vector containing FALSE and an error message is returned.

checkHaplotypeSyntax <- function(haplotype, filename){

  alleles <- strsplit(haplotype, "-|~")[[1]]

  alleleSyntaxCheck<-sapply(alleles, function(x) NULL)

  #extract results from checkAlleleSyntax -- make results into a list for each
  #allele
  for(i in 1:length(alleles)){
    #set variant type to haplotype for checkLocusDataset function
    alleleSyntaxCheck[[i]]<-list(checkAlleleSyntax(alleles[[i]], filename))
  }

  #if there are any FALSEs in alleleSyntaxCheck, return error message
  if(any(unlist(alleleSyntaxCheck)==FALSE)){

    #unnest list by one structure
    alleleSyntaxCheck<-alleleSyntaxCheck %>%
      flatten()

    #look at elements in list to see which alleles contain FALSE
    alleleErrors<-names(alleleSyntaxCheck)[sapply(1:length(alleleSyntaxCheck), function(x) "FALSE" %in% alleleSyntaxCheck[[x]])]

    #paste allele names with associated errors. If there is more than one error
    #all errors will be pasted together and returned as one vector
    return(c(FALSE, paste(paste(names(alleleSyntaxCheck[c(alleleErrors)]), sapply(alleleSyntaxCheck[c(alleleErrors)], "[", 2), sep =" - "), collapse=" ")))
  }
  return(TRUE)
}
