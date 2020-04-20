###Syntactic and semantic validation of HLA amino acid motifs v1 16APR20
#'Syntactic and semantic validation of HLA amino acid motifs
#'
#'Checks input motif for errors in format and amino acid positions not present in the locus alignment.
#'
#'@param motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.
#'
#'@return A warning message if the input motif is formatted incorrectly, or contains an amino acid position not present in the alignment. Otherwise, a list object with extracted locus information, a correctly formatted motif, and locus specific amino acid dataframe are returned. Note checkMotif() does not check amino acid variants in a specified motif; that is done by findMotif().
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
#'\donttest{checkMotif("DRB1*26F~28E~30Y")}
#'
#'#Example where format is incorrect
#'\donttest{checkMotif("DRB1**26F~28E~30Y")}
#'
#'#Example where an amino acid position does not exist
#'\donttest{checkMotif("DRB1**26F~28E~300000Y")}

checkMotif<-function(motif){

  #if conditions to catch if a motif is formatted incorrectly
  if(any(is.na(str_extract(strsplit(strsplit(motif, "*", fixed=T)[[1]][2], "~")[[1]], "[A-Z]"))==TRUE)){
    return(warning("Your motif is formatted incorrectly. Please use the Locus*##$~##$~##$ format, where ## identifies a peptide position, and $ identifies an amino acid residue."))}

  if((is.na(str_extract(strsplit(motif, "*", fixed=T)[[1]][2], "[A-Z]"))==TRUE) | ((str_count(strsplit(motif, "*", fixed=T)[[1]][2], "[A-Z]")>=2) & (grepl("~", strsplit(strsplit(motif, "*", fixed=T)[[1]][2], "*"))==FALSE)) |(length(strsplit(motif, "*", fixed=T)[[1]]) > 2) | (length(strsplit(motif, "*", fixed=T)[[1]])==1)){
    return(warning("Your motif is formatted incorrectly. Please use the Locus*##$~##$~##$ format, where ## identifies a peptide position, and $ identifies an amino acid residue."))
  }

  #extract loci information
  split_1 <- unlist(strsplit(motif,"*",fixed="true"))
  loci <- split_1[1]
  motifs <- unlist(strsplit(split_1[2],"~",fixed=TRUE))

  #if IMGTprotalignments does not exist (i.e not previously already downloaded and in the local
  #environment) then generate HLAalignments)
  #obtains HLAalignments df
  if(!exists("IMGTprotalignments")){
    HLAalignments<-BLAASD(loci)}

  #if IMGTprot alignments does exist, use the locus specific alignment
  #if locus is set as DRB1, use "DRB" as locus specific alignment
  if(exists("IMGTprotalignments")){
    locus<-loci
    if((locus %in% names(SSHAARP::IMGTprotalignments))==FALSE){
      return(warning(paste(locus, "is not a valid locus.")))
    }
    HLAalignments<-SSHAARP::IMGTprotalignments[[locus]]
  }


  #if HLAalignments is not a list because BLAASD() output is an error, return
  #HLAalignments, which contains the error
  if(is.list(HLAalignments)==FALSE){
    return(HLAalignments)
  }

  #examines motifs to make sure amino acid positions are in the correct order -- sorts numerically
  #if they are not
  motifs <- mixedsort(motifs)

  #examines if amino acid positions in the motif are present in the alignment - returns error message if one or more positions is not in the alignment
  if(!all(substr(motifs,1,nchar(motifs)-1) %in% colnames(HLAalignments)[5:ncol(HLAalignments)])) {
    return(warning("One or more of your amino acid positions is not present in the alignment. Please make sure amino acid positions of interest are present in the current release of ANHIG/IMGTHLA alignments."))
  }

  #return a list object with loci and motifs information
  #if no error message is returned
  else{
    return(list(loci, motifs, HLAalignments))}
}
