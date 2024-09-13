##findMotif v2.0.0 3JAN2022
#'Returns an alignment data frame of alleles that share a specific amino acid motif
#'
#'Returns an alignment data frame of alleles that share a specific amino acid motif.
#'
#'@param motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the Solberg dataset.
#'
#'@return An amino acid alignment dataframe of alleles that share the specified motif. Otherwise, a vector containing FALSE and an error message is returned.
#'
#'@export
#'
#'@examples
#'
#'#example with existing motif
#' \dontrun{findMotif("DRB1*26F~28E~30Y", filename=solberg_dataset)}
#' \dontrun{findMotif("DRB1*26F~28E", filename=solberg_dataset)}
#'
#'#example with non-existent motif
#'\dontrun{findMotif("DRB1*26F~28E~30Z", filename=solberg_dataset)}
#'
#'#extracting names of alleles with user-defined motif
#'findMotif("DRB1*26F~28E~30Y", filename=solberg_dataset)[,4]

findMotif<-function(motif, filename){

  #check if input motif is formatted correctly or if amino acid position
  #is present in the alignment
  positionCheck<-checkPosition(motif, filename)

  if(grepl(FALSE, positionCheck[[1]])){
    return(c(FALSE, positionCheck[[2]]))
  }

  locus<-getVariantInfo(motif)[[1]]
  motifs<-getVariantInfo(motif)[[2]]

  #examines motifs to make sure amino acid positions are in the correct order -- sorts numerically
  #if they are not
  motifs <- mixedsort(motifs)

  #enters HLAalignments information from check_results
  HLAalignments<-SSHAARP::IMGTprotalignments[[locus]]

  for(x in 1:length(motifs)) {
    HLAalignments <- HLAalignments[HLAalignments[substr(motifs[x],1,nchar(motifs[x])-1)]==substr(motifs[x],nchar(motifs[x]),nchar(motifs[x])),]

    if(nrow(HLAalignments)==0)
    {
      return(c(FALSE, paste(motif, "No alleles possess this motif", sep=": ")))
    }
  }

  #if motifs are found, HLAalignments[[loci[[i]]]] is returned
  return(HLAalignments)
}
