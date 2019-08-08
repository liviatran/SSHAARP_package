#'Solberg dataset manipulation
#'
#'Manipulates the Solberg dataset by adding in a new column containing locus*allele information, reordering the dataset based on population name, and subsetting the data to only the locus of interest.
#'
#'@param dataset The Solberg dataset, a comma-separated value (CSV) file, which is the '1-locus-alleles.dat' file in the results.zip archive at http://pypop.org/popdata/. See the vignette for more details on pertinent parts of the Solberg dataset used for this package.
#'@param motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.
#'
#'@importFrom utils read.delim
#'
#'@note For internal SSHAARP use only.
#'@export
#'
#'@return A dataframe format of the Solberg dataset, where the rows are ordered by population name, and limited to populations that have information on the locus of interest.
#'
#'
dataSubset<-function(dataset, motif){
  #reads in Solberg DS
  solberg_DS<-as.data.frame(read.delim(dataset), stringsAsFactors=F)

  #makes a new column with locus and trimmed allele pasted together named locus_allele
  solberg_DS$locus_allele<-paste(solberg_DS$locus, solberg_DS$allele_v3, sep="*")

  #orders solberg-DS by population
  solberg_DS<-solberg_DS[order(solberg_DS$popname),]

  solberg_DS[,]<-sapply(solberg_DS[, ], as.character)

  #subsets the Solberg_DS to only the locus of interest
  solberg_DS<-subset(solberg_DS, solberg_DS$locus==strsplit(motif, "\\*")[[1]][1])

  return(solberg_DS)}
