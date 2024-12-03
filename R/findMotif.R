##findMotif v 2.0.4 25NOV2024
#'Returns an alignment data frame of alleles that share a specific amino acid motif
#'
#'Returns an alignment data frame of alleles that share a specific amino acid motif.
#'
#'@param motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the Solberg dataset.
#'@param alignments A list object of sub-lists of data frames of protein alignments for the HLA and HLA-region genes supported in the ANHIG/IMGTHLA GitHub Repository. Alignments will always be the most recent version IPD-IMGT/HLA Database version.
#'
#'@note For internal SSHAARP use only.
#'
#'@return An amino acid alignment dataframe of alleles that share the specified motif. Otherwise, a vector containing FALSE and an error message is returned.
#'
#'@export

findMotif<-function(motif, filename, alignments){

  locus<-getVariantInfo(motif)[[1]]

  #check if input motif is formatted correctly or if amino acid position
  #is present in the alignment
  positionCheck<-checkPosition(motif, filename, alignments[[locus]])

  if(grepl(FALSE, positionCheck[[1]])){
    return(c(FALSE, positionCheck[[2]]))
  }

  motifs<-getVariantInfo(motif)[[2]]

  #examines motifs to make sure amino acid positions are in the correct order -- sorts numerically
  #if they are not
  motifs <- mixedsort(motifs)

  HLAalignment<-alignments[[locus]]

  for(x in 1:length(motifs)) {
    HLAalignment <- HLAalignment[HLAalignment[substr(motifs[x],1,nchar(motifs[x])-1)]==substr(motifs[x],nchar(motifs[x]),nchar(motifs[x])),]
  }

  if(nrow(HLAalignment)==0){
    return(c(FALSE, paste(motif, "No alleles possess this motif", sep=": ")))
  }
  #if motifs are found, HLAalignment[[loci[[i]]]] is returned
  return(HLAalignment)
}
