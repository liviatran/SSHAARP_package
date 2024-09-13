###readFilename v 2.0.0 3JAN2022
#'Designates dataset by either reading in file the user has provided, or using the Solberg dataset or mock haplotype dataset
#'
#'Returns the user specified dataset, either by reading in the file the user has provided, or by using the Solberg dataset or mock haplotype dataset.
#'
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset or mock haplotype dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the datasets bundled with the package. Allele and motif datasets should follow the Solberg dataset format, and haplotype datasets should follow the SSHAARP haplotype mock data format.
#'
#'@importFrom tools file_ext
#'@importFrom utils read.csv
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return A dataframe of the user specified dataset.

readFilename<-function(filename){

  #if filename is a data frame, then it it is the Solberg dataset, or the bundled mock haplotype dataset
  if(is.data.frame(filename)) {
    dataset <- filename
  }

  #read in file if not bundled Solberg dataset or mock haplotype dataset
  else{
    #use read.delim if .dat or .txt file
    if(file_ext(filename) %in% c("dat", "txt")){
      dataset <- data.frame(read.delim(filename, stringsAsFactors = F, sep=""), stringsAsFactors=F)
    }
    #use read.csv if .csv file
    if(file_ext(filename) %in% "csv"){
      dataset <- data.frame(read.csv(filename, stringsAsFactors = F), stringsAsFactors=F)
    }
  }
  return(dataset)
}
