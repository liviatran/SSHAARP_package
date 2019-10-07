#'Syntactic and semantic validation of HLA amino acid motifs
#'
#'Checks input motif for errors in format and amino acid positions not present in the locus alignment.
#'
#'@param motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.
#'
#'@return An error message if the input motif is formatted incorrectly, or contains an amino acid position not present in the alignment. Otherwise, a list object with extracted locus information, a correctly formatted motif, and locus specific amino acid dataframe are returned.
#'
#'@importFrom stringr str_count
#'@importFrom gtools mixedsort
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@examples
#'#Example where a motif is formatted correctly
#'checkMotif("DRB1*26F~28E~30Y")
#'
#'#Example where format is incorrect
#'checkMotif("DRB1**26F~28E~30Y")
#'
#'#Example where an amino acid position does not exist
#'checkMotif("DRB1**26F~28E~300000Y")

checkMotif<-function(motif){

  #if conditions to catch if a motif is formatted incorrectly
  ifelse(is.na(str_extract(strsplit(strsplit(motif, "*", fixed=T)[[1]][2], "~")[[1]], "[A-Z]"))==TRUE, return("Your motif is formatted incorrectly. Please use the Locus*##$~##$~##$ format, where ## identifies a peptide position, and $ identifies an amino acid residue."), "")
  if((is.na(str_extract(strsplit(motif, "*", fixed=T)[[1]][2], "[A-Z]"))==TRUE) | ((str_count(strsplit(motif, "*", fixed=T)[[1]][2], "[A-Z]")>=2) & (grepl("~", strsplit(strsplit(motif, "*", fixed=T)[[1]][2], "*"))==FALSE)) |(length(strsplit(motif, "*", fixed=T)[[1]]) > 2) | (length(strsplit(motif, "*", fixed=T)[[1]])==1)){
    return("Your motif is formatted incorrectly. Please use the Locus*##$~##$~##$ format, where ## identifies a peptide position, and $ identifies an amino acid residue.")
  }

  #extract loci information
  split_1 <- unlist(strsplit(motif,"*",fixed="true"))
  loci <- split_1[1]
  motifs <- unlist(strsplit(split_1[2],"~",fixed=TRUE))

  #if AA_segments does not exist (i.e not previously already downloaded and in the local
  #environment) then generate AA_segments)
  if(!exists("AA_segments")){
    #obtains AA_segments df
    AA_segments<-BLAASD(loci)
  }

  #if AA_segments it not a list because BLAASD() output is an error, return
  #AA_segments, which contains the error
  if(is.list(AA_segments)==FALSE){
    return(AA_segments)
  }

  #examines motifs to make sure amino acid positions are in the correct order -- sorts numerically
  #if they are not
  motifs <- mixedsort(motifs)

  #examines if amino acid positions in the motif are present in the alignment - returns error message if one or more positions is not in the alignment
  if(!all(substr(motifs,1,nchar(motifs)-1) %in% colnames(AA_segments[[loci]])[5:ncol(AA_segments[[loci]])])) {
    return("One or more of your amino acid positions is not present in the alignment. Please make sure amino acid positions of interest are present in the current release of IPD-IMGT/HLA alignments.")
  }

  #return a list object with loci and motifs information
  #if no error message is returned
  else{
    return(list(loci, motifs, AA_segments[[loci]]))}
}
