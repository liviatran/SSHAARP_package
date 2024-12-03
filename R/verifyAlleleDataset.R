#verifyAlleleDataset v 2.0.4 25NOV2024
#'Verifies the allele entered is present in specified dataset
#'
#'Verifies the allele entered is present in the specified dataset.
#'
#'@param allele An allele name written in the IPD-IMGT/HLA Database format.
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the Solberg dataset.
#'@param alignments A list object of sub-lists of data frames of protein alignments for the HLA and HLA-region genes supported in the ANHIG/IMGTHLA GitHub Repository. Alignments will always be the most recent version IPD-IMGT/HLA Database version.
#'
#'@importFrom BIGDAWG GetField
#'@importFrom dplyr %>%
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return TRUE if the allele is present in the specified data set, and the filtered allele dataset. If a user enters an allele with more than two fields and has selected the Solberg dataset as the data source, a message informing the user that the allele has been truncated is appended to the output. If an allele entered is valid, but is not present in the user provided dataset, a warning message is returned.

verifyAlleleDataset<-function(allele, filename, alignments){

  locus<-getVariantInfo(allele)[[1]]

  #check allele syntax
  alleleCheck<-checkAlleleSyntax(allele, filename)

  if(grepl(TRUE, alleleCheck[[1]]) == FALSE){
    return(c(FALSE, alleleCheck[[2]]))
  }

  alleleANHIGcheck<-verifyAlleleANHIG(allele, filename, alignments)

  #check allele is in IMGTprotalignments
  if(grepl(TRUE, alleleANHIGcheck[[1]])==FALSE){
    return(c(FALSE, alleleANHIGcheck[[2]]))
  }

  #filters dataset to entries with relevant locus
  datasetAlleles<-dataSubset(allele, filename)

  #returns FALSE and error message if locus in allele is invalid
  if(is.data.frame(datasetAlleles)==FALSE){
   return(warning(datasetAlleles))
  }

#if user specified data set is the Solberg dataset
  if(all(filename == SSHAARP::solberg_dataset)){

    endMess<-NULL

    #if allele entered has more than two fields, truncate
    if(length(strsplit(strsplit(allele, "*", fixed=T)[[1]][2], ":")[[1]]) > 2){
      endMess<-paste("The specified allele name for", allele, "was truncated to two fields, since the Solberg dataset, which only includes two field alleles, was selected as a data source.")
      allele<-GetField(allele, Res=2)
  }

  filteredAllele<-datasetAlleles[datasetAlleles$locus_allele %in% allele,]

  if(nrow(filteredAllele) == 0 & length(endMess)!= 0){
   return(list(FALSE, append(endMess, paste("The", allele, "allele entered is a valid allele and exists in the IMGT protein alignments, but is not present in the Solberg dataset."))))
  }

  if(nrow(filteredAllele) == 0){
    return(list(FALSE, paste("The", allele, "allele entered is a valid allele and exists in the IMGT protein alignments, but is not present in the Solberg dataset.")))
  }

  else{
   return(list(TRUE, endMess, data.frame(filteredAllele)))
   }
  }

#for any other user specified dataset
  else{
   filteredAllele<- datasetAlleles[datasetAlleles$locus_allele %in% allele,]

    if(nrow(filteredAllele) == 0){
     return(c(FALSE, paste("The", allele, "allele entered is a valid allele and exists in IMGTprotalignments, but is not present in the dataset."), NULL))
     }
    }
  return(list(TRUE, NA, data.frame(filteredAllele)))

}
