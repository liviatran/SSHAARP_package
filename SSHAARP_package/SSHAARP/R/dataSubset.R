##Solberg dataset manipulation v1 16APR20
#'Solberg dataset manipulation
#'
#'Returns a modified version of the Solberg dataset that includes a column of locus*allele names, is sorted by by population name, and is reduced to the specified locus. Cardinal coordinates are converted to their Cartesian equivalents (i.e. 50S is converted to -50).
#'
#'@param filename The filename of the local copy of the Solberg dataset - the defaulted filename is the solberg_dataset in the SSHAARP package.
#'@param motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.
#'
#'@importFrom utils read.delim
#'@importFrom stringr str_extract
#'
#'@note For internal SSHAARP use only.
#'@note The Solberg dataset is the tab-delimited ‘1-locus-alleles.dat’ text file in the results.zip archive at http://pypop.org/popdata/.
#'@note The Solberg dataset is also prepackaged into SSHAARP as 'solberg_dataset'.
#'@export
#'
#'@return A data frame containing a reformatted version of the Solberg dataset, with rows ordered by population name, Cartesian coordinates in the latit and longit columns, and limited to populations with data for the specified locus. If a motif has formatting errors, a warning message is returned.
#'
#'
dataSubset<-function(motif, filename=SSHAARP::solberg_dataset){

  if(is.data.frame(filename)) {
    solberg_DS <- filename
  } else {solberg_DS <- as.data.frame(read.delim(filename), stringsAsFactors=F)}

  #checks input motif for formatting errors
  check_results<-suppressWarnings(checkMotif(motif))

  #if length of check_results is an error, return the error
  if(length(check_results)<2){
    return(warning(check_results))
  }

  #makes a new column with locus and trimmed allele pasted together named locus_allele
  solberg_DS$locus_allele<-paste(solberg_DS$locus, solberg_DS$allele_v3, sep="*")

  #orders solberg-DS by population
  solberg_DS<-solberg_DS[order(solberg_DS$popname),]

  solberg_DS[,]<-sapply(solberg_DS[, ], as.character)

  #subsets the Solberg_DS to only the locus of interest
  solberg_DS<-subset(solberg_DS, solberg_DS$locus==check_results[[1]])

  if(any((grepl("S", solberg_DS$latit))==TRUE)){
    solberg_DS$latit[which((grepl("S", solberg_DS$latit))==TRUE)]<-as.numeric(paste("-", str_extract(solberg_DS$latit[which((grepl("S", solberg_DS$latit))==TRUE)],"\\-*\\d+\\.*\\d*"), sep=""))}

  if(any((grepl("N", solberg_DS$latit))==TRUE)){
    solberg_DS$latit[which((grepl("N", solberg_DS$latit))==TRUE)]<-as.numeric(str_extract(solberg_DS$latit[which((grepl("N", solberg_DS$latit))==TRUE)],"\\-*\\d+\\.*\\d*"))
  }

  #longitude conversions
  if(any((grepl("W", solberg_DS$longit))==TRUE)){
    solberg_DS$longit[which((grepl("W", solberg_DS$longit))==TRUE)]<-as.numeric(paste("-", str_extract(solberg_DS$longit[which((grepl("W", solberg_DS$longit))==TRUE)],"\\-*\\d+\\.*\\d*"), sep=""))}

  if(any((grepl("E", solberg_DS$longit))==TRUE)){
    solberg_DS$longit[which((grepl("E", solberg_DS$longit))==TRUE)]<-as.numeric(str_extract(solberg_DS$longit[which((grepl("E", solberg_DS$longit))==TRUE)],"\\-*\\d+\\.*\\d*"))
  }

  return(solberg_DS)}
