#'Narrows AA_segments table to only alleles that have the amino acid motif of interest
#'
#'Isolates the individual amino acid position dataframe produced by AA_segments_maker to only alleles with the user-defined motif. If the user-defined motif does not correspond to any alleles, an error message is output.
#'
#'@param motif An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.
#'
#'@return If the motif is found in the AA_segments dataframe, a dataframe with alleles that have the user input amino acid motif is returned. If the motif is not found, a dataframe is returned, where one column has the motif, and the other column contains an error message.
#'
#'@importFrom gtools mixedsort
#'@importFrom BIGDAWG GetField
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

findMotif<-function(motif){

  #if conditions to catch if a motif is formatted incorrectly
  if((length(strsplit(motif, "*", fixed=T)[[1]]) > 2) | (length(strsplit(motif, "*", fixed=T)[[1]])==1)){
    return("Your motif is formatted incorrectly. Please use the Locus*##$~##$~##$ format, where ## identifies a peptide position, and $ identifies an amino acid residue.")
  }

  #extract loci information
  loci<-strsplit(motif, "\\*")[[1]][1]

  #if AA_segments does not exist (i.e not previously already downloaded and in the local
  #environment) then generate AA_segments)
  if(!exists("AA_segments")){
    #obtains AA_segments df
    AA_segments<-buildAAsegments(loci)}

  #since "DRB" is used as the search criteria for the alignment (IMGTHLA/ANHIG groups all DRB loci
  #into one alignment, AA_segments consists of all DRB loci, not just DRB1)
  #if the loci is DRB1, this conditional statement subsets AA_segments to only DRB1 loci,
  #and if "NA" is present in the locus column for the alignment sequence coordinate row
  if(loci=="DRB1"){
    AA_segments[[loci]]<-subset(AA_segments[[loci]], (loci==AA_segments[[loci]]$locus) | (is.na(AA_segments[[loci]]$locus)))
  }

  #examines motifs to make sure amino acid positions are in the correct order -- sorts numerically
  #if they are not
  motif<-paste(loci, sep="*",paste(mixedsort(strsplit(strsplit(motif, "*", fixed=T)[[1]][[2]], "~")[[1]]), collapse="~"))


  #for loop for searching amino acid motifs
  #three total loops are run, where subsequent position*motifs are evaluated against which alleles are present
  #after the previous motif subset
  #each AA_segments is bound with the alignment coordinates except for the last AA_segments
  for(t in 1:length(strsplit(strsplit(motif, "*", fixed=T)[[1]][[2]], "~")[[1]])){
    AA_segments[[loci]]<-AA_segments[[loci]][which((AA_segments[[loci]][5:ncol(AA_segments[[loci]])][which((str_extract(strsplit(strsplit(motif,"*",fixed=TRUE)[[1]][2],"~",fixed=TRUE)[[1]], "[0-9]+")[[t]]==names(AA_segments[[loci]][1,5:ncol(AA_segments[[loci]])]))==TRUE)]==str_extract(strsplit(strsplit(motif,"*",fixed=TRUE)[[1]][2],"~",fixed=TRUE)[[1]],"[A-Z]")[[t]])==TRUE),]
  }

  #if no motifs are found, a warning message is thrown
  if((nrow(AA_segments[[loci]])==0)){
    return(data.frame("Motif"=motif, "Error message"="No alleles match this motif"))
  }

  #if motifs are found, AA_segments[[loci[[i]]]] is returned
  if((nrow(AA_segments[[loci]])!=0)){
    return(AA_segments[[loci]])}
}

