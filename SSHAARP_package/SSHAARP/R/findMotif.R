#'Reduce an alignment to alleles sharing a specific amino acid motif
#'
#'Consumes the alignment data frame produced by BLAASD() and returns an alignment data frame of alleles that share a specific amino acid motif.
#'
#'@param motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.
#'
#'@return An amino acid alignment dataframe of alleles that share the specified motif. If the motif is not found in any alleles, a data frame containing an error message is returned.
#'
#'@importFrom gtools mixedsort
#'@importFrom BIGDAWG GetField
#'@importFrom stringr str_count str_extract
#'
#'@export
#'
#'@examples
#'
#'#example with actual motif
#' findMotif("DRB1*26F~28E~30Y")
#' findMotif("DRB1*26F~28E")
#'
#'#example with non-existent motif
#'findMotif("DRB1*26F~28E~30Z")
#'
#'#extracting names of alleles with user-defined motif
#'findMotif("DRB1*26F~28E~30Y")[,4]

#library(gtools)
#library(BIGDAWG)
#library(stringr)

findMotif<-function(motif){
  #if conditions to catch if a motif is formatted incorrectly
  if(((str_count(strsplit(motif, "*", fixed=T)[[1]][2], "[A-Z]")>=2) & (grepl("~", strsplit(strsplit(motif, "*", fixed=T)[[1]][2], "*"))==FALSE)) |(length(strsplit(motif, "*", fixed=T)[[1]]) > 2) | (length(strsplit(motif, "*", fixed=T)[[1]])==1)){
    return("Your motif is formatted incorrectly. Please use the Locus*##$~##$~##$ format, where ## identifies a peptide position, and $ identifies an amino acid residue.")
  }

  #extract loci information
  split_1 <- unlist(strsplit(motif,"*",fixed="true"))
  loci <- split_1[1]
  motifs <- unlist(strsplit(split_1[2],"~",fixed=TRUE))

  #if AA_segments does not exist (i.e not previously already downloaded and in the local
  #environment) then generate AA_segments)
  if(!exists("AA_segments")){
    #obtains AA_segments df
    AA_segments<-BLAASD(loci)
  }

  #examines motifs to make sure amino acid positions are in the correct order -- sorts numerically
  #if they are not
  motifs <- mixedsort(motifs)

  #examines if amino acid positions in the motif are present in the alignment - returns error message if one or more positions is not in the alignment
  if(!all(substr(motifs,1,nchar(motifs)-1) %in% colnames(AA_segments[[loci]])[5:ncol(AA_segments[[loci]])])) {
    return("One or more of your amino acid positions is not present in the alignment. Please make sure amino acid positions of interest are present in the current release of IPD-IMGT/HLA alignments.")
  }

  #since "DRB" is used as the search criteria for the alignment (IMGTHLA/ANHIG groups all DRB loci
  #into one alignment, AA_segments consists of all DRB loci, not just DRB1)
  #if the loci is DRB1, this conditional statement subsets AA_segments to only DRB1 loci,
  #and if "NA" is present in the locus column for the alignment sequence coordinate row
  if(loci=="DRB1"){
    AA_segments[[loci]]<-subset(AA_segments[[loci]], (loci==AA_segments[[loci]]$locus) | (is.na(AA_segments[[loci]]$locus)))
  }

  #for loop for searching amino acid motifs
  #three total loops are run, where subsequent position*motifs are evaluated against which alleles are present
  #after the previous motif subset
  #each AA_segments is bound with the alignment coordinates except for the last AA_segments
  for(x in 1:length(motifs)) {
    AA_segments[[loci]] <- AA_segments[[loci]][AA_segments[[loci]][substr(motifs[x],1,nchar(motifs[x])-1)]==substr(motifs[x],nchar(motifs[x]),nchar(motifs[x])),]
    if(nrow(AA_segments[[loci]])==0)
    {
      return(data.frame("Motif"=motif, "Error message"="No alleles possess this motif"))
    }
  }

  #if motifs are found, AA_segments[[loci[[i]]]] is returned
  return(AA_segments[[loci]])
}
