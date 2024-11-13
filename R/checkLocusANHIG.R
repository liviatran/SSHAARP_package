###checkLocusANHIG v2.0.1 2024OCT7
#'Check locus validity
#'
#'Checks if the locus in the entered variant is a protein-coding gene annotated by the IPD-IMGT/HLA Database
#'
#'@param variant An amino acid motif or allele with an HLA locus name followed by an asterisk. This function ONLY evaluates if the locus in the entered variant is a protein-coding gene.
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return TRUE if locus in entered variant are in the IPD-IMGT/HLA Database. Otherwise, a vector with FALSE and an error message is returned.
#'
#'@examples
#'#Example of valid locus in a motif
#'checkLocusANHIG("DRB1*26F~28E")
#'#[1] TRUE
#'
#'#Example of an invalid locus in an allele
#'checkLocusANHIG("BOO*01:01")
#'#[1] "FALSE"                     "BOO is not a valid locus."

checkLocusANHIG<-function(variant){

  locus <-getVariantInfo(variant)[[1]]
  if(any((locus %in% c('A', 'B', 'C', 'DMA', 'DMB', 'DOA', 'DOB', 'DPA1', 'DPB1', 'DQA1', 'DQA2', 'DQB1', 'DRA', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'E', 'F', 'G', 'HFE', 'MICA', 'MICB', 'TAP1', 'TAP2'))==FALSE)){
      return(c(FALSE, paste(locus, "is not a valid locus.")))
  }

  return(TRUE)
}
