#'Returns an alignment data frame of alleles that share a specific amino acid motif
#'
#'Consumes the alignment data frame produced by BLAASD() and returns an alignment data frame of alleles that share a specific amino acid motif.
#'
#'@param input_motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.
#'
#'@return An amino acid alignment dataframe of alleles that share the specified motif. If the motif is not found in any alleles, a data frame containing an error message is returned.
#'
#'@importFrom BIGDAWG GetField
#'
#'@export
#'
#'@examples
#'
#'#example with actual motif
#' findMotif("DRB1*26F~28E~30Y")
#' findMotif("DRB1*26F~28E")
#'
#'#example with non-existent motif
#'findMotif("DRB1*26F~28E~30Z")
#'
#'#extracting names of alleles with user-defined motif
#'findMotif("DRB1*26F~28E~30Y")[,4]

#library(gtools)
#library(BIGDAWG)
#library(stringr)

findMotif<-function(input_motif){

  #check if input motif is formatted correctly or if amino acid position
  #is present in the alignment
  check_results<-checkMotif(input_motif)

  #if length of check_results is an error, return the error
  if(length(check_results)<2){
    return(check_results)
  }

  #enters loci information from check_results
  loci<-check_results[[1]]

  #enters motifs information from check_results
  motifs<-check_results[[2]]

  #enters AA_segments information from check_results
  AA_segments<-check_results[[3]]

  #since "DRB" is used as the search criteria for the alignment (IMGTHLA/ANHIG groups all DRB loci
  #into one alignment, AA_segments consists of all DRB loci, not just DRB1)
  #if the loci is DRB1, this conditional statement subsets AA_segments to only DRB1 loci,
  #and if "NA" is present in the locus column for the alignment sequence coordinate row
  if(loci=="DRB1"){
    AA_segments<-subset(AA_segments, (loci==AA_segments$locus) | (is.na(AA_segments$locus)))
  }

  for(x in 1:length(motifs)) {
    AA_segments <- AA_segments[AA_segments[substr(motifs[x],1,nchar(motifs[x])-1)]==substr(motifs[x],nchar(motifs[x]),nchar(motifs[x])),]
    if(nrow(AA_segments)==0)
    {
      return(data.frame("Motif"=input_motif, "Error message"="No alleles possess this motif"))
    }
  }

  #if motifs are found, AA_segments[[loci[[i]]]] is returned
  return(AA_segments)}
