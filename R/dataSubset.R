##Dataset manipulation for motifs and alleles v 2.0.3 20NOV2024
#'Dataset manipulation for motifs and alleles
#'
#'Returns a modified version of the user selected dataset that includes a column of locus*allele names, is sorted by by population name, and is reduced to the specified locus. Cardinal coordinates are converted to their Cartesian equivalents (i.e. 50S is converted to -50).
#'
#'@param variant An allele or an amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.
#'@param filename The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the Solberg dataset.
#'
#'@importFrom utils read.delim
#'@importFrom stringr str_extract
#'
#'@note For internal SSHAARP use only.
#'@note The Solberg dataset is the tab-delimited ‘1-locus-alleles.dat’ text file in the results.zip archive at http://pypop.org/popdata/.
#'@note The Solberg dataset is also prepackaged into SSHAARP as 'solberg_dataset'.
#'
#'@export
#'
#'@return A data frame containing a reformatted version of the user selected dataset, with rows ordered by population name, Cartesian coordinates in the latit and longit columns, and limited to populations with data for the specified locus. Otherwise, a vector containing FALSE and an error message is returned.

dataSubset<-function(variant, filename){

  dataset<-readFilename(filename, variant)

  locusANHIG<-checkLocusANHIG(variant)

  if(any(grepl(FALSE, locusANHIG[[1]]))==TRUE){
    return(c(FALSE, locusANHIG[[2]]))
  }

  locusDS<-checkLocusDataset(variant, filename)

  if(grepl(FALSE, locusDS[[1]])==TRUE){
    return(c(FALSE, locusDS[[2]]))
  }

  #if locus is foundin ANHIG and in the dataset, extract locus information
  if(grepl(TRUE, locusANHIG[[1]]) & grepl(TRUE, locusDS[[1]])){
    locus<-getVariantInfo(variant)[[1]]
  }

  #makes a new column with locus and trimmed allele pasted together named locus_allele
  dataset$locus_allele<-paste(dataset$locus, dataset$allele_v3, sep="*")

  #orders solberg-DS by population
  dataset<-dataset[order(dataset$popname),]

  dataset[,]<-sapply(dataset[,], as.character)

  #subsets the dataset to only the locus of interest
  dataset<-dataset[dataset$locus == locus,]

  if(any((grepl("S", dataset$latit))==TRUE)){
    dataset$latit[which((grepl("S", dataset$latit))==TRUE)]<-as.numeric(paste("-", str_extract(dataset$latit[which((grepl("S", dataset$latit))==TRUE)],"\\-*\\d+\\.*\\d*"), sep=""))}

  if(any((grepl("N", dataset$latit))==TRUE)){
    dataset$latit[which((grepl("N", dataset$latit))==TRUE)]<-as.numeric(str_extract(dataset$latit[which((grepl("N", dataset$latit))==TRUE)],"\\-*\\d+\\.*\\d*"))
  }

  #longitude conversions
  if(any((grepl("W", dataset$longit))==TRUE)){
    dataset$longit[which((grepl("W", dataset$longit))==TRUE)]<-as.numeric(paste("-", str_extract(dataset$longit[which((grepl("W", dataset$longit))==TRUE)],"\\-*\\d+\\.*\\d*"), sep=""))}

  if(any((grepl("E", dataset$longit))==TRUE)){
    dataset$longit[which((grepl("E", dataset$longit))==TRUE)]<-as.numeric(str_extract(dataset$longit[which((grepl("E", dataset$longit))==TRUE)],"\\-*\\d+\\.*\\d*"))
  }

  return(dataset)
}
