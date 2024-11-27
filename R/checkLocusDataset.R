###checkLocusDataset v2.0.4 26NOV2024
#'Check locus validity and if the locus is present in the user specified dataset
#'
#'Checks if the locus in the entered variant is a protein-coding gene annotated by the IPD-IMGT/HLA Database, and if it is in the user specified dataset.
#'
#'@param variant An allele or an amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids. Haplotypes must contain alleles that follow the aforementioned format, and may be delimited by "~" or "-".
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset or mock haplotype dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the datasets bundled with the package. Allele and motif datasets should follow the Solberg dataset format, and haplotype datasets should follow the SSHAARP haplotype mock data format.
#'
#'@importFrom tools file_ext
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return TRUE if locus is a protein-coding gene and is in the specified dataset. Otherwise, a vector with FALSE and an error message is returned.

checkLocusDataset<-function(variant, filename){

  locus <-getVariantInfo(variant)[[1]]

  locusANHIGcheck<-checkLocusANHIG(variant)

  if(grepl(FALSE, locusANHIGcheck[[1]])){
    return(c(FALSE, locusANHIGcheck[[2]]))
  }

  dataset<-readFilename(filename, variant)

  #if column name "haplotype" is not present in the column name of datasets,
  #then look at the locus column to evaluate whether the locus of the entered
  #variant is in the dataset
  if(!"haplotype" %in% colnames(dataset)){
    if(grepl(TRUE, locusANHIGcheck[[1]]) & (locus %in% dataset$locus) == FALSE){
      return(c(FALSE, paste(locus, "is a valid locus, but is not in the user selected dataset")))
    }
  }
  #otherwise, search for an exact match of the locus of the entered variant in the
  #haplotype column of the dataset
  else{
    if(grepl(TRUE, locusANHIGcheck[[1]]) & (all(grepl(paste("\\<", locus, "\\>", sep=""), dataset$haplotype)==FALSE))){
      return(c(FALSE, paste(locus, "is a valid locus, but is not in the user selected dataset")))
    }
  }

  return(TRUE)
}
