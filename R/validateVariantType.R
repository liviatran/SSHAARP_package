## validateVariantType v 2.0.7 11JUNE2025
#'validateVariantType
#'
#'Validates the variantType parameter in PALM()
#'
#'@param vT Specifies whether the variant is an allele, motif, or haplotype.
#'
#'@note For internal SSHAARP use only
#'
#'@return TRUE if the entered variantType is accepted ('motif', 'allele', or 'haplotype'). FALSE for any other provided entry.
#'
#'@export
#'
validateVariantType<-function(vT){

  if(!vT %in% c('allele', 'haplotype', 'motif')){
    return(FALSE)
  }

  return(TRUE)
}
