#checkAlleleSyntax v2.0.3 18NOV2024
#'Check allele syntax
#'
#'Checks if allele syntax is valid.
#'
#'@param allele An allele name written in the IPD-IMGT/HLA Database format.
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the Solberg dataset.
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return TRUE if allele syntax is correct. Otherwise,a vector containing FALSE and an error message is returned.

checkAlleleSyntax<-function(allele, filename){

  locusANHIG<-checkLocusANHIG(allele)

  if(grepl(FALSE, locusANHIG[[1]])){
    return(c(FALSE, locusANHIG[[2]]))
  }

  #allow variantType entry to be informed by user
  locusDS<-checkLocusDataset(allele, filename)

  if(grepl(FALSE, locusDS[[1]])){
    return(c(FALSE,locusDS[[2]]))
  }

  #default isNamePresent variantType to allele
  nameValid<-isNamePresent(allele, variantType = "allele")

  if(!grepl(TRUE, nameValid[[1]])){
    return(c(FALSE, nameValid[[2]]))
  }

  #1)count if allele has between more than 4 fields
  #2)count the number of characters in a given field after splitting at ":" delimiter
  #3)after split at *, see if ":" is present in the second vector containing field information (ex DRB1*01)
  #if it is less than 2, return warning message
  #4)# of fields are the same as # of of ":" occurrences

  if((length(strsplit(allele, ":")[[1]]) > 4) | any(nchar(strsplit(strsplit(allele, "*", fixed=T)[[1]], ":")[[2]]) < 2) | !grepl(":", strsplit(allele, "*", fixed=T)[[1]][[2]]) | length(strsplit(strsplit(allele, "*", fixed=TRUE)[[1]][2], ":")[[1]]) == str_count(allele, ":")){
    return(c(FALSE, "The syntax for the allele entered is incorrect. Please ensure the allele entered follows the IPD-IMGT/HLA Database format."))
  }
  else{
    return(TRUE)
  }
}
