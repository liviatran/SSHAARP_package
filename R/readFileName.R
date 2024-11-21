###readFilename v 2.0.3. 2024NOV20
#'Designates dataset by either reading in file the user has provided, or using the Solberg dataset or mock haplotype dataset
#'
#'Returns the user specified dataset, either by reading in the file the user has provided, or by using the Solberg dataset or mock haplotype dataset. If the user provides a dataset and the filename is not found, an error will be returned. If the user provided dataset does not have the same number of columns or the same column names as the reference datasets, an error message will be returned.
#'
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset or mock haplotype dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the datasets bundled with the package. Allele and motif datasets should follow the Solberg dataset format, and haplotype datasets should follow the SSHAARP haplotype mock data format.
#'@param variant An allele or an amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids. Haplotypes must contain alleles that follow the aforementioned format, and may be delimited by "~" or "-".
#'
#'@importFrom tools file_ext
#'@importFrom utils read.csv
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return A dataframe of the user specified dataset.

readFilename<-function(filename, variant){

  #if filename is a data frame, then it it is the Solberg dataset, or the bundled mock haplotype dataset
  if(is.data.frame(filename)) {
    dataset <- filename
  }

  #read in file if not bundled Solberg dataset or mock haplotype dataset
  else{
    #check if filename/path exists
    if(!file.exists(filename)){
      return("The input filename does not exist. Please double check the file path entered for the 'filename' parameter")
    }

    #use read.delim if .dat or .txt file
    if(file_ext(filename) %in% c("dat", "txt")){
      dataset <- data.frame(read.delim(filename, stringsAsFactors = F, sep=""), stringsAsFactors=F)
    }
    #use read.csv if .csv file
    if(file_ext(filename) %in% "csv"){
      dataset <- data.frame(read.csv(filename, stringsAsFactors = F), stringsAsFactors=F)
    }

    if(variant == 'haplotype'){
      if(ncol(dataset) != ncol(SSHAARP::mock_haplotype_dataset)){
        return('The user provided dataset does not contain the same number of columns as the mock haplotype dataset. Please make sure the input dataset follows the same format as the bundled haplotype dataset.')
      }
      if(!all(colnames(dataset) %in% colnames(SSHAARP::mock_haplotype_dataset))){
        return(sprintf('The following column names were detected in the user input file, but are not present in the mock haplotype dataset: %s. Please make sure the input dataset follows the same format as the bundled haplotype dataset.', paste(colnames(dataset)[!colnames(dataset) %in% colnames(SSHAARP::mock_haplotype_dataset)], collapse = ',')))
      }
    } else{
        if(ncol(dataset) != ncol(SSHAARP::solberg_dataset)){
          return('The user provided dataset does not contain the same number of columns as the Solberg dataset. Please make sure the input dataset follows the same format as the bundled Solberg dataset.')
        }
        if(!all(colnames(dataset) %in% colnames(SSHAARP::solberg_dataset))){
          return(sprintf('The following column names were detected in the user input file, but are not present in the Solberg dataset: %s. Please make sure the input dataset follows the same format as the bundled Solberg dataset.', paste(colnames(dataset)[!colnames(dataset) %in% colnames(SSHAARP::solberg_dataset)], collapse = ',')))        }
    }
  }
  return(dataset)
}
