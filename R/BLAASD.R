##BLAASD - Build Loci Amino Acid Specific Dataframe v1 16APR20
#'BLAASD - Build Loci Amino Acid Specific Dataframe
#'
#'Extracts alignment sequence information for a given locus from the ANHIG/IMGTHLA GitHub repository to produce a dataframe of individual amino acid data for each amino acid position for all alleles, for a user-defined HLA locus or loci. The first 4 columns are locus, allele, trimmed allele, and allele_name.
#'
#'@param loci A vector of un-prefixed HLA locus names
#'
#'@return A list object of data frames for each specified locus. Each list element is a data frame of allele names and the corresponding peptide sequence for each amino acid position. An error message is return if the loci input is not a locus for which petpide alignments are available in the ANHIG/IMGTHLA Github Repository.
#'
#'@importFrom stringr str_squish str_split
#'@importFrom utils head tail capture.output
#'@importFrom BIGDAWG GetField
#'@importFrom dplyr filter %>%
#'@export
#'
#'@examples
#'#BLAASD with one locus as input
#'\donttest{BLAASD("C")}
#'
#'#BLAASD with multiple loci as input
#'\donttest{BLAASD(c("A", "B", "C"))}
BLAASD<-function(loci){

  #checks if input locus is present in version 3.38.0 HLA loci
  #skip name checks for DRB1/3/4/5, as they are a part of the DRB alignment
  for(j in 1:length(loci)){
    if(loci[j]=="DRB1"|loci[j]=="DRB3"|loci[j]=="DRB4"|loci[j]=="DRB5") next
    if(loci[j] %in% names(SSHAARP::IMGTprotalignments)== FALSE){
      return(warning(paste(loci[j], "is not a valid locus.")))
    }}

  #creates empty variables for future for loops
  start<-end<-alignment<-list()

  #creates empty variables where each element is named after the used loci

  #empty variables for correspondence table
  pepsplit<-refexon<-AA_aligned<-HLAalignments<-inDels<-corr_table<-cols<-downloaded_segments<-w<-alignment_positions<-alignment_length<-alignment_start<-prot_extractions<-refblock_number<-end_char<-space_diff<-sapply(loci, function(x) NULL)


  for(i in 1:length(loci)){

    #downloads relevant locus alignment file -- readLines allows for space preservation, which is important in
    #finding where the alignment sequence starts
    #tryCatch() to ensure loci are input correctly
    alignment[[loci[i]]] <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/",paste(ifelse(loci[[i]]=="DRB1"|loci[[i]]=="DRB3"|loci[[i]]=="DRB4"|loci[[i]]=="DRB5","DRB",loci[[i]]),"_prot.txt",sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)

    #alters alignment file by cutting out non-pertinent information in beginning
    #and endind of alignment file
    alignment[[loci[i]]] <- head(alignment[[loci[i]]],-3)
    alignment[[loci[i]]] <- tail(alignment[[loci[i]]],-7)

    #see countSpaces function at beginning of script
    #Counts difference between Prot to -30 and beginning of Prot to -30 + 1 due to zero number indexing to find where
    #the alignment sequence actually starts
    space_diff[[loci[i]]]<-(nchar(strsplit(alignment[[loci[i]]][3], " ")[[1]][2])+countSpaces(alignment[[loci[i]]][3])[2]+1)-countSpaces(alignment[[loci[i]]][2])[1]

    #reduces repeated whitespace in alignment file and removes rows with empty values for proper
    #start and stop subsetting
    alignment[[loci[i]]] <-str_squish(alignment[[loci[i]]])
    alignment[[loci[i]]] <-alignment[[loci[i]]][-which(alignment[[loci[i]]] == "")]

    #determines positions of "Prot" and the end of that reference block segment
    start[[loci[i]]]<-as.numeric(grep("Prot", alignment[[loci[i]]]))
    end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,length(alignment[[loci[i]]])))

    #counts number of characters in the very last allele to add onto the last Prot enumeration block
    #to obtain end length
    end_char[[loci[i]]]<-nchar(sapply(strsplit(gsub(" ", "", sub(" ", "~", str_squish(tail(alignment[[loci[i]]], 1)))), "~"), "[", 2))-1

    #extracts rows with "Prot" and reference sequence position information
    #extracts only relevant reference sequence positions
    #NOTE: the first row containing "Prot" contains two numbers -- -30 and 1 -- where only -30, is extracted,
    #as the actual sequence start will always be 1
    for (j in 1:length(start[[loci[i]]])){

      prot_extractions[[loci[i]]][j]<-strsplit(alignment[[loci[i]]][start[[loci[i]]][j]], " ")

      refblock_number[[loci[i]]][j]<-as.numeric(sapply(prot_extractions[[loci[i]]][j], "[", 2))


      #determines the alignment start by adding -30 to the difference between white spaces found above
      alignment_start[[loci[i]]]<-refblock_number[[loci[i]]][1]+space_diff[[loci[i]]]
    }

    #closes all white space in the alignment file, except for the white space separating the allele and peptide sequence
    alignment[[loci[i]]] <-paste(substr(alignment[[loci[i]]],1,regexpr(" ",text = alignment[[loci[i]]],fixed = TRUE)), gsub(" ","",substr(alignment[[loci[i]]],regexpr(" ",text = alignment[[loci[i]]],fixed = TRUE),nchar(alignment[[loci[i]]]))),sep = "")

    #string splits at white spaces to yield allele and peptide sequences
    alignment[[loci[i]]]  <- strsplit(alignment[[loci[i]]]," ", fixed=T)

    #binds the previously split strings by row
    alignment[[loci[i]]] <- do.call(rbind,alignment[[loci[i]]])

    #if the pepseq column is equal to the allele column due to premature peptide termination,
    #insert a blank in place of the allele in the pepseq column
    alignment[[loci[i]]][which(alignment[[loci[i]]][,1]==alignment[[loci[i]]][,2]),2] <- ""

    #renames columns to "alleles" and "pepseq"
    colnames(alignment[[loci[i]]])<-c(paste(loci[[i]], "alleles", sep="_"), "pepseq")

    #due to ANHIG formatting, cases where an allele contains newly reference peptide sequences will not
    #contain the same number of rows as previous reference peptide blocks
    #this for loop is invoked to add "."for all other alleles for each character in the newly reference peptide
    #to preserve structural integrity
    #changes 10/9/19 to accommodate if there is more than one extraneous allele with an extended amino acid sequence
    for(k in 1:length(start[[loci[i]]])){
      if(nrow(alignment[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],])){
        x<-as.data.frame(alignment[[loci[i]]][,1][start[[loci[i]]][1]:end[[loci[i]]][1]][-c(1,2)], stringsAsFactors = F)
        colnames(x)<-paste(loci[[i]], "alleles", sep="_")
        x<-cbind.data.frame(x, pepseq=as.character(paste(rep(".", nchar(tail(alignment[[loci[i]]][,2], 1))), collapse = "")), stringsAsFactors=FALSE)
        y<-data.frame(tail(alignment[[loci[i]]], (nrow(alignment[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],][nrow(alignment[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), stringsAsFactors = F)
        x$pepseq[match(y[,1], x[,1])]<-y$pepseq
        alignment[[loci[i]]]<-as.matrix(rbind(head(alignment[[loci[i]]], -(nrow(alignment[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],][nrow(alignment[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), x))
        start[[loci[i]]]<-as.numeric(grep("Prot", alignment[[loci[i]]]))
        end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,nrow(alignment[[loci[i]]])))}
    }
    #if a locus has extra formatting, resulting in unqeual rows, start and end will be updated to reflect subsetting
    #if a locus has no extra formatting, start and end will remain the same, as procured by earlier code
    for(e in 1:length(start[[loci[i]]])){
      HLAalignments[[loci[i]]]<-cbind(HLAalignments[[loci[i]]], alignment[[loci[i]]][start[[loci[i]]][e]:end[[loci[i]]][e],])}

    #removes first two rows containing AA position and "Prot"
    HLAalignments[[loci[i]]] <- HLAalignments[[loci[i]]][-c(1,2),]

    #designates columns to be combined as every other so allele names are not included
    #in pasting all the amino acid sequences together
    cols<-seq(0, ncol(HLAalignments[[loci[i]]]), by=2)
    HLAalignments[[loci[i]]]<-cbind(HLAalignments[[loci[i]]][,1], apply(HLAalignments[[loci[i]]][,cols], 1 ,paste, collapse = ""))

    #creates a new matrix with the number of columns equal to the number of characters in the reference sequence
    corr_table[[loci[i]]]<-matrix(, nrow = 2, ncol = as.numeric(nchar(HLAalignments[[loci[i]]][,2][1])))

    #if the first position enumeration is negative (i.e. has a leader peptide sequence), determines alignment length based on the total number of characters plus the alignment start (which is negative)
    if(grepl("-", alignment_start[[loci[i]]][[1]])==TRUE){
      alignment_length[[loci[i]]]<-as.numeric(nchar(HLAalignments[[loci[i]]][,2][1]))+alignment_start[[loci[[i]]]]}

    #if there is no leader peptide (i.e sequence starts at 1), determines alignment length based on total number of characters
    else{
      alignment_length[[loci[i]]]<-as.numeric(nchar(HLAalignments[[loci[i]]][,2][1]))
    }

    #pastes alignment_start to alignment_length together in sequential order, with inDels accounted for
    #captures output as "w"
    w[[i]] <- capture.output(cat(alignment_start[[loci[i]]]:alignment_length[[loci[i]]]))

    #splits string formed by cat for separate character variables
    alignment_positions[[loci[i]]]<-as.character(unlist(strsplit(w[[loci[i]]], " ")))

    #eliminates "0" if present in the alignment positions, as the alignment sequence from ANHIG does not contain 0
    if("0" %in% alignment_positions[[loci[i]]]==TRUE){
      alignment_positions[[loci[i]]]<-alignment_positions[[loci[i]]][-which(alignment_positions[[loci[i]]] == 0)]
    }

    #contains alignment sequence information
    corr_table[[loci[i]]][2,]<-alignment_positions[[loci[i]]]

    #string splits to extract locus in the allele name
    #assigns to new variable "AA_aligned"
    AA_aligned[[loci[i]]]<- as.matrix(do.call(rbind,strsplit(HLAalignments[[loci[i]]][,1],"[*]")))

    #adds a new column of pasted locus and trimmed two field alleles to AA_aligned
    AA_aligned[[loci[i]]]<- cbind(AA_aligned[[loci[i]]], paste(AA_aligned[[loci[i]]][,1], apply(AA_aligned[[loci[i]]],MARGIN=c(1,2),FUN=GetField,Res=2)[,2], sep="*"))

    #binds AA_aligned and HLAalignments -- renames columns
    HLAalignments[[loci[i]]] <- cbind(AA_aligned[[loci[i]]], HLAalignments[[loci[i]]])
    colnames(HLAalignments[[loci[i]]]) <- c("locus", "full_allele", "trimmed_allele", "allele_name", "AAsequence")

    #if locus is DRB3/4/5, use DRB line 1 as reference sequence
    if(loci[[i]]=="DRB3"|loci[[i]]=="DRB4"|loci[[i]]=="DRB5"){
      #sets refexon to a reference peptide for each HLA locus based on the reference sequences in HLAalignments
      refexon[[loci[i]]] <- rbind(HLAalignments[[loci[i]]][1,])[which(rbind(HLAalignments[[loci[i]]][1,])[,"locus"]=="DRB1"),'AAsequence']
    }

    else{
      refexon[[loci[i]]] <- rbind(HLAalignments[[loci[i]]][1,])[which(rbind(HLAalignments[[loci[i]]][1,])[,"locus"]==loci[[i]]),'AAsequence']}

    #if input loci is DRB, use grep statement to match input loci to loci in HLAalignments
    if(loci[[i]]=="DRB"){
      refexon[[loci[i]]]<-rbind(HLAalignments[[loci[i]]][1,])[grepl(loci[[i]], rbind(HLAalignments[[loci[i]]][1,])[,"locus"]), "AAsequence"]}

    #splits AA_sequence column at every amino acid, resulting in a list of split amino acids for each row
    pepsplit[[loci[i]]] <- sapply(HLAalignments[[loci[i]]][,"AAsequence"],strsplit,split="*")

    #fills in space with NA for alleles with premature termination to make it the same number of characters
    #as the reference sequence
    pepsplit[[loci[i]]]<- lapply(pepsplit[[loci[i]]],function(x) c(x,rep("NA",nchar(refexon[[loci[i]]])-length(x))))

    #binds pep_split together by element in its previous list form by row
    pepsplit[[loci[i]]]<- do.call(rbind,pepsplit[[loci[i]]])

    #nullifies row names
    rownames(pepsplit[[loci[i]]]) <- NULL

    #binds all columns together to form desired ouput, as described above
    HLAalignments[[loci[i]]] <- cbind.data.frame(HLAalignments[[loci[i]]][,1:4],pepsplit[[loci[i]]], stringsAsFactors=FALSE)

    #get reference sequence from DRB alignment
    if(loci[[i]]=="DRB3"|loci[[i]]=="DRB4"|loci[[i]]=="DRB5"){
      DRBref<-HLAalignments[[loci[[i]]]][1,]}

    #if the locus is DRB1/3/4/5, subset the specific locus from DRB alignment
    if(loci[[i]]=="DRB1"|loci[[i]]=="DRB3"|loci[[i]]=="DRB4"|loci[[i]]=="DRB5"){
      HLAalignments[[loci[i]]]<-HLAalignments[[loci[i]]][HLAalignments[[loci[[i]]]]$locus==loci[[i]],]}

    #add reference sequence to DRB3/4/5, reset row names to numerical order
    if(loci[[i]]=="DRB3"|loci[[i]]=="DRB4"|loci[[i]]=="DRB5"){
      HLAalignments[[loci[[i]]]]<- rbind(DRBref, HLAalignments[[loci[[i]]]])
      rownames(HLAalignments[[loci[[i]]]])<-NULL}

    #finds positions in HLAalignments that have ".", indicating an inDel
    inDels[[loci[[i]]]]<-colnames(HLAalignments[[loci[[i]]]][1, 5:ncol(HLAalignments[[loci[[i]]]])][HLAalignments[[loci[[i]]]][1, 5:ncol(HLAalignments[[loci[[i]]]])] %in% "."])

    #inputs HLAalignments alignment sequence into the corr_table with "InDel" still present
    corr_table[[loci[[i]]]][1,]<-names(HLAalignments[[loci[[i]]]][5:ncol(HLAalignments[[loci[[i]]]])])

    #inDel inclusion if there are inDels present
    if(length(inDels[[loci[[i]]]])!=0){
      for(b in 1:length(inDels[[loci[[i]]]])){
        corr_table[[loci[[i]]]][2,][inDels[[loci[[i]]]][[b]]==corr_table[[loci[[i]]]][1,]]<-paste("INDEL", b, sep="-")
      }
    }

    #if alignment start is position 1, alignment start does not need to be accounted for
    #when determining length of corr_table in re-enumerating corr_table with InDels
    if(alignment_start[[loci[[i]]]]==1){
      #fixes enumerations following "InDel"
      corr_table[[loci[[i]]]][2,][!grepl("INDEL", corr_table[[loci[[i]]]][2,])]<-(alignment_start[[loci[[i]]]]:((length(corr_table[[loci[[i]]]][2,])-length(corr_table[[loci[[i]]]][2,][grepl("INDEL", corr_table[[loci[[i]]]][2,])]))))[!(alignment_start[[loci[[i]]]]:((length(corr_table[[loci[[i]]]][2,])-length(corr_table[[loci[[i]]]][2,][grepl("INDEL", corr_table[[loci[[i]]]][2,])]))))==0]
    }

    else{
      corr_table[[loci[[i]]]][2,][!grepl("INDEL", corr_table[[loci[[i]]]][2,])]<-(alignment_start[[loci[[i]]]]:((length(corr_table[[loci[[i]]]][2,])-length(corr_table[[loci[[i]]]][2,][grepl("INDEL", corr_table[[loci[[i]]]][2,])]))+alignment_start[[loci[[i]]]]))[!(alignment_start[[loci[[i]]]]:((length(corr_table[[loci[[i]]]][2,])-length(corr_table[[loci[[i]]]][2,][grepl("INDEL", corr_table[[loci[[i]]]][2,])]))+alignment_start[[loci[[i]]]]))==0]
    }

    #renames columns in HLAalignments
    colnames(HLAalignments[[loci[i]]]) <- c("locus","allele","trimmed_allele","allele_name", corr_table[[loci[[i]]]][2,])

    #distributes  reference sequence from row 1
    #into all other rows, if they contain a "-"
    #amino acids with changes will not be impacted
    for(k in 5:ncol(HLAalignments[[loci[i]]])) {
      HLAalignments[[loci[i]]][,k][which(HLAalignments[[loci[i]]][,k]=="-")] <- HLAalignments[[loci[i]]][,k][1]}

    if(loci[[i]]=="DQB1"){
      HLAalignments[[loci[i]]] <- specCirc(HLAalignments[[loci[i]]],loci[i])
    }
  }
  return(HLAalignments)
}

#This function was obtained from Reddit and written by Josh Bredeweg, user jbraids1421.
countSpaces <- function(x){
  counter <- 0
  coll <- numeric()
  vec <- strsplit(x," ")[[1]]
  for(i in 1:length(vec)){
    if (vec[i]==""){
      counter <- counter+1
    }
    else{
      if (counter!=0) coll <- c(coll,counter)
      counter <- 1
    }
  }
  coll
}


#This function was obtained from Stack Overflow, written by user thelatemail.
dupdiff <- function(x,y){
  x[-match(
    make.unique(as.character(y)),
    make.unique(as.character(x)),
    nomatch=0
  )]}

specCirc <- function(alObj, locus){
  AA_atlas<-SSHAARP::AA_atlas
  ## Special Custom Rules for Locus-Specific Circumstances

  if(locus=="DQB1"){

    #renames all columns from InDel 5 to rest of alignment to numerical positions
    names(alObj)[(grep(AA_atlas$DQB1$Boundary[[4]] , colnames(alObj))+1):length( alObj)]<-(AA_atlas$DQB1$Boundary[[4]]+1):((as.numeric(colnames( alObj)[ncol( alObj)])+AA_atlas$DQB1$Boundary[[5]]-AA_atlas$DQB1$Boundary[[4]])+length(names( alObj)[(grep(AA_atlas$DQB1$Boundary[[4]] , colnames( alObj))+1):length( alObj)])-length((AA_atlas$DQB1$Boundary[[4]]+1):(as.numeric(colnames( alObj)[ncol( alObj)])+AA_atlas$DQB1$Boundary[[5]]-AA_atlas$DQB1$Boundary[[4]])))

    range<-1:ncol(alObj)

    #extracts exon 5  DQB1*05:03:01:01 in lieu of using reference sequence
    exon_5<-cbind(alObj[,4, drop=F], alObj[,(range[colnames(alObj) %in% AA_atlas$DQB1$Boundary[[4]]]+1): range[colnames(alObj) %in% AA_atlas$DQB1$Boundary[[5]]]]) %>% filter(alObj$allele_name=="DQB1*05:03:01:01")
    exon_5$allele_name<-NULL

    ## checkhere for e5 size
    if(length(exon_5) > 8) { #### SM
      ##finds last InDel in exon 5
      lastInDel<-as.numeric(gsub("InDel_", "", colnames( alObj[grep("INDEL", colnames( alObj))][length(alObj[grep("INDEL", colnames( alObj))])])))

      #compares exon 5 sequence to DQB1*05:03:01:01 to see account for any future InDels
      InDeldiff<-dupdiff(str_split(paste(exon_5, collapse=""), "")[[1]],c("P","Q","G","P","P","P","A","G"))

      #conditions if more than one InDel is present
      ##      if(length(InDeldiff==".")>1){
      for(j in 1:length(InDeldiff)){ #### SM shortened this as it was extremley redundant
        names(exon_5)[grep(names(exon_5[exon_5 %in% "."])[j], names(exon_5))]<-paste("INDEL", (lastInDel+1):(lastInDel+length(InDeldiff)), sep="_")[[j]]  ## shortened this too, redundant
      }
      ##      } #### There is no point in separating 1 "." from multiple "."s, as j loops from 1 to N
      #### replacing original line 382
      lastPos <- AA_atlas$DQB1[4,]$Boundary[[1]]

      for(j in 1:ncol(exon_5)) {
        if(substr(colnames(exon_5)[j],1,1)!="I") {
          lastPos <- lastPos+1
          colnames(exon_5)[j] <- as.character(lastPos) }
      }


      names(alObj)[(range[colnames( alObj) %in% AA_atlas$DQB1$Boundary[[4]]]+1):range[colnames( alObj) %in% AA_atlas$DQB1$Boundary[[5]]]]<-names(exon_5)

      for(j in (range[colnames(alObj) %in% as.character(lastPos)]+1):(length(range))) {
        if(substr(colnames(alObj)[j],1,1)!="I") {
          lastPos <- lastPos +1
          colnames(alObj)[j] <- as.character(lastPos) }
      }
    }
  } ## end of the check on line 364
  return(alObj)
}
