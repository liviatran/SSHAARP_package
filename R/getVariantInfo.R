### Locus and allele or motif extraction v2.0.4 26NOV2024
#'Locus and allele or motif extraction for motifs, allele, or haplotype mapping
#'
#'Extracts locus and allele or motif information from variant.
#'
#'@param variant An amino acid motif or allele. Amino acid motifs must be in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Alleles must have 1-4 fields.
#'
#'@note For internal SSHAARP use only.
#'
#'@export
#'
#'@return A list object with loci and allele or motif information.

getVariantInfo<-function(variant){

    split_1 <- unlist(strsplit(variant,"*",fixed=TRUE))
    loci <- split_1[1]
    variants <- unlist(strsplit(split_1[2],"~",fixed=TRUE))

  return(list(loci, variants))
}
