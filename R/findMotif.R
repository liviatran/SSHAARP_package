##Returns an alignment data frame of alleles that share a specific amino acid motif v1 16APR20
#'Returns an alignment data frame of alleles that share a specific amino acid motif
#'
#'Consumes the alignment data frame produced by BLAASD() and returns an alignment data frame of alleles that share a specific amino acid motif.
#'
#'@param motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.
#'
#'@return An amino acid alignment dataframe of alleles that share the specified motif. If the motif is not found in any alleles, or the motif has formatting errors, a warning message is returned.
#'
#'@importFrom BIGDAWG GetField
#'
#'@export
#'
#'@examples
#'
#'#example with actual motif
#' \donttest{findMotif("DRB1*26F~28E~30Y")}
#' \donttest{("DRB1*26F~28E")}
#'
#'#example with non-existent motif
#'\donttest{findMotif("DRB1*26F~28E~30Z")}
#'
#'#extracting names of alleles with user-defined motif
#'findMotif("DRB1*26F~28E~30Y")[,4]

findMotif<-function(motif){

  #check if input motif is formatted correctly or if amino acid position
  #is present in the alignment
  check_results<-suppressWarnings(checkMotif(motif))

  #if length of check_results is an error, return the error
  if(length(check_results)<2){
    return(warning(check_results))
  }

  #enters loci information from check_results
  loci<-check_results[[1]]

  #enters motifs information from check_results
  motifs<-check_results[[2]]

  #enters HLAalignments information from check_results
  HLAalignments<-check_results[[3]]

  for(x in 1:length(motifs)) {
    HLAalignments <- HLAalignments[HLAalignments[substr(motifs[x],1,nchar(motifs[x])-1)]==substr(motifs[x],nchar(motifs[x]),nchar(motifs[x])),]

    if(nrow(HLAalignments)==0)
    {
      return(warning(paste(motif, "No alleles possess this motif", sep=" : ")))
    }
  }

  #if motifs are found, HLAalignments[[loci[[i]]]] is returned
  return(HLAalignments)}
