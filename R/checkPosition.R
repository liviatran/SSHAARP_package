#checkPosition v 2.0.4 25NOV2024
#'Checks if amino acid positions in motif exist
#'
#'Checks if amino acid positions in the entered motif exist in IMGTprotalignments.
#'
#'@param motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids. This function ONLY checks if the entered amino acid positions exist in IMGTprotalignments.
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the Solberg dataset.
#'@param alignments A list object of sub-lists of data frames of protein alignments for the HLA and HLA-region genes supported in the ANHIG/IMGTHLA GitHub Repository. Alignments will always be the most recent version IPD-IMGT/HLA Database version.
#'
#'@importFrom gtools mixedsort
#'@importFrom utils data
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return TRUE if all of the amino acid positions in a motif exist. Otherwise, a vector with FALSE and an error message is returned.

checkPosition<-function(motif, filename, alignments){

  motifCheck<-checkMotifSyntax(motif, filename)

  if(grepl(TRUE, motifCheck[[1]])==FALSE){
  return(c(FALSE, motifCheck[[2]]))
}

  locus<-getVariantInfo(motif)[[1]]
  motifs<-getVariantInfo(motif)[[2]]

  #examines motifs to make sure amino acid positions are in the correct order -- sorts numerically
  #if they are not
  motifs <- mixedsort(motifs)

  #examines if amino acid positions in the motif are present in the alignment
  if(!all(substr(motifs,1,nchar(motifs)-1) %in% colnames(alignments)[5:ncol(alignments)])) {
    return(c(FALSE, "One or more of your amino acid positions is not present in the alignment. Please make sure amino acid positions of interest are present in the current release of ANHIG/IMGTHLA alignments."))
  }

  return(TRUE)
}
