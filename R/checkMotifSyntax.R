#checkMotifSyntax v 2.0.3 20NOV2024
#'Check motif syntax
#'
#'Checks if motif syntax is valid.
#'
#'@param motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the Solberg dataset.
#'
#'@importFrom stringr str_extract str_count
#'@importFrom stringi stri_split_fixed
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'
#'@return TRUE if the motif syntax is valid. Otherwise, a vector containing FALSE and an error message is returned.
#'
#'@examples
#'
#'#Example with correct motif syntax where user specified dataset is the Solberg dataset
#'checkMotifSyntax("DRB1*26F~28E~30Y", filename=SSHAARP::solberg_dataset)
#'
#'#Example with incorrect motif syntax where user specified dataset is the Solberg dataset
#'checkMotifSyntax("DRB1****26F~28E", filename=SSHAARP::solberg_dataset)

checkMotifSyntax<-function(motif, filename){

  #check if locus in motif is valid
  locusANHIG<-checkLocusANHIG(motif)

  if(grepl(FALSE, locusANHIG[[1]])==TRUE){
    return(c(FALSE, locusANHIG[[2]]))
  }

  #check if locus in motif is present in user selecte dataset
  locusDS<-checkLocusDataset(motif, filename)

  if(grepl(FALSE, locusDS[[1]])==TRUE){
    return(c(FALSE,locusDS[[2]]))
  }

  nameValid<-isNamePresent(motif, variantType = "motif")

  if(grepl(FALSE, nameValid[[1]])==TRUE){
    return(c(FALSE, nameValid[[2]]))
  }

  #if conditions to catch if a motif is formatted incorrectly
  if(grepl(":", motif) | grepl("[a-z]", motif) | any((str_count(strsplit(stri_split_fixed(str = motif, pattern = "*", n = 2)[[1]], "~")[[2]], '\\*|[A-Z]')>=2)==TRUE) | any(is.na(str_extract(strsplit(strsplit(motif, "*", fixed=T)[[1]][2], "~")[[1]], "[A-Z]"))==TRUE) | ((str_count(strsplit(motif, "*", fixed=T)[[1]][2], "[A-Z]")>=2) & (grepl("~", strsplit(strsplit(motif, "*", fixed=T)[[1]][2], "*"))==FALSE)) |(length(strsplit(motif, "*", fixed=T)[[1]]) > 2) | (length(strsplit(motif, "*", fixed=T)[[1]])==1) | (  length(strsplit(strsplit(motif, "*", fixed=TRUE)[[1]][[2]], "~")[[1]]) == str_count(strsplit(motif, "*", fixed=TRUE)[[1]][[2]], "~"))){
    return(c(FALSE, "Your motif is formatted incorrectly. Please use the Locus*##$~##$~##$ format, where ## identifies a peptide position, and $ identifies an amino acid residue."))
  }
  else{
    return(TRUE)
  }
}

