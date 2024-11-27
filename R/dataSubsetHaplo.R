##Dataset manipulation for haplotypes v 2.0.4 25NOV2024
#'Dataset manipulation for haplotypes
#'
#'Returns the user input dataset that contains the selected haplotype.
#'
#'@param haplotype A haplotype where allele names are written in the IPD-IMGT/HLA Database format, and have 1-4 fields. Alleles in haplotypes may be delimited by "-" or "~".
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or pre-bundled mock haplotype dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the mock haplotype dataset bundled with the package.
#'@param AFND A logical parameter that determines whether the user specified dataset is data from AFND. This parameter is only relevant if haplotype maps are being made.
#'@param alignments A list object of sub-lists of data frames of protein alignments for the HLA and HLA-region genes supported in the ANHIG/IMGTHLA GitHub Repository. Alignments will always be the most recent version IPD-IMGT/HLA Database version.
#'
#'@importFrom utils read.delim
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return A two element list with 1) a subset data frame containing only haplotypes with alleles present in the user input haplotype, and 2) a data frame of the full dataset. Alleles with two fields will be evaluated with their three and four field allele equivalents, and alleles with three fields will be evaluated with their four field allele equivalent. Otherwise, a vector containing FALSE and an error message is returned.

dataSubsetHaplo<-function(haplotype, filename, AFND, alignments){

  dataset<-readFilename(filename)

  #if dataset is AFND dataset, divide frequency by 100
  if(AFND == TRUE){
    dataset$frequency<- dataset$frequency/100
  }

  alleleANHIGcheck<-verifyAlleleANHIGHaplo(haplotype, filename, alignments)

  #check allele is in IMGTprotalignments
  if(grepl(TRUE, alleleANHIGcheck[[1]])==FALSE){
    return(c(FALSE, alleleANHIGcheck[[2]]))
  }

  alleles<-strsplit(haplotype, "-|~")[[1]]

  #filters any haplotypes that have any of the alleles present in the user input haplotype
  #if two or three field alleles are present in the haplotype, 3 and 4, and 4 field alleles,
  #respectively, will be used to map  haplotype frequency
  filtered_ds<-dataset %>%
    filter(grepl(paste(paste("(?=.*", gsub("*", "\\*",alleles, fixed=T), ")", sep=""), collapse=""), dataset$haplotype, perl=T))

  if(nrow(filtered_ds)==0){
    return(c(FALSE, "The combination of alleles in the entered haplotype are valid and found in the latest release of ANHIG/IMGT protein alignments, but are not present in the user input dataset."))
  }

  return(list(filtered_ds, dataset))
}
