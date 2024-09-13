#checkPosition v 2.0.0 3JAN2022
#'Checks if amino acid positions in motif exist
#'
#'Checks if amino acid positions in the entered motif exist in IMGTprotalignments.
#'
#'@param motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids. This function ONLY checks if the entered amino acid positions exist in IMGTprotalignments.
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the Solberg dataset.
#'
#'@importFrom gtools mixedsort
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return TRUE if all of the amino acid positions in a motif exist. Otherwise, a vector with FALSE and an error message is returned.
#'
#'@examples
#'#Example with existent amino acid positions
#'checkPosition("DRB1*26F~28E", filename=SSHAARP::solberg_dataset)
#'
#'#Example with nonexistent amino acid positions
#'checkPosition("DRB1*199999F", filename=SSHAARP::solberg_dataset)

checkPosition<-function(motif, filename){

  motifCheck<-checkMotifSyntax(motif, filename)

  if(grepl(TRUE, motifCheck[[1]])==FALSE){
  return(c(FALSE, motifCheck[[2]]))
}

  locus<-getVariantInfo(motif)[[1]]
  motifs<-getVariantInfo(motif)[[2]]

  #examines motifs to make sure amino acid positions are in the correct order -- sorts numerically
  #if they are not
  motifs <- mixedsort(motifs)

  HLAalignments<-SSHAARP::IMGTprotalignments[[locus]]

  #examines if amino acid positions in the motif are present in the alignment
  if(!all(substr(motifs,1,nchar(motifs)-1) %in% colnames(HLAalignments)[5:ncol(HLAalignments)])) {
    return(c(FALSE, "One or more of your amino acid positions is not present in the alignment. Please make sure amino acid positions of interest are present in the current release of ANHIG/IMGTHLA alignments."))
  }

    return(TRUE)
}
