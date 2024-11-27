#isNamePresent v2.0.4 26NOV2024
#'Checks if name portion of entered variant is present
#'
#'Checks if the name portion of the entered variant is present. Names consist of information following the locus and asterisk of the entered variant.
#'
#'@param variant An amino acid motif or allele. Amino acid motifs must be in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Alleles must have 1-4 fields.
#'@param variantType Identifies whether the variant is an allele or motif.
#'
#'@importFrom stringr str_extract
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return TRUE if name is present. Otherwise, a vector with FALSE and an error message is returned.

isNamePresent<-function(variant, variantType){

  if(variantType=="motif"){
    if(is.na(str_extract(strsplit(variant, "*", fixed=T)[[1]][2], "[A-Z]")) | is.na(str_extract(strsplit(variant, "*", fixed=T)[[1]][2], "[0-9]"))){
      return(c(FALSE, paste0("The amino-acid motif is missing in", sep=" ", variant, ". Please use the LOCUS*MOTIF format.", sep="")))
    }
  }

  if(variantType=="allele"){
    if(is.na(str_extract(strsplit(variant, "*", fixed=T)[[1]][2], "[A-Z0-9]"))){
      return(c(FALSE, paste0("The allele name is missing in", sep=" ", variant, ". Please use the LOCUS*ALLELE-NAME format.", sep="")))
    }
  }

  return(TRUE)

}
