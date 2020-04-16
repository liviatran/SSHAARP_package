## Solberg Dataset v1 16APR20
#' Solberg Dataset
#'
#' A dataframe of the original Solberg dataset, which is a global dataset of 497 population samples from 185 published datasets, representing 66,800 individuals.
#' For more information on the Solberg dataset, please see the vignette.
#' @docType data
#' @name solberg_dataset
#' @usage solberg_dataset
#' @format A dataframe with 20163 rows and 13 columns.
#' @source results.zip file from http://pypop.org/popdata/
#' @references Solberg et.al "Balancing selection and heterogeneity across the classical human leukocyte antigen loci: A meta-analytic review of 497 population studies". Human Immunology (2008) 69, 443–464
"solberg_dataset"

#' Protein alignments for all protein coding genes in the IPD-IMGT/HLA Database release v 3.39.0.
#'
#' A list object containing protein alignments for all protein coding genes in the IPD-IMGT/HLA Database release. Alignments were downloaded from the ANHIG/IMGTHLA Github respository.
#' @docType data
#' @name IMGTprotalignments
#' @usage IMGTprotalignments
#' @format A list containing protein alignments in a dataframe format for each locus.
#' @source https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments
"IMGTprotalignments"

#' Exon boundaries for all exons in protein coding genes in the IPD-IMGT/HLA Database release v 3.39.0
#'
#' A list object containing exon boundaries for all exons in protein coding genes in the IPD-IMGT/HLA Database release v 3.39.0. Exon boundaries were determined from the nucleotide alignment files, which were downloaded from the ANHIG/IMGTHLA Github respository.
#' @docType data
#' @name AA_atlas
#' @usage AA_atlas
#' @format A list containing exon boundaries in a dataframe format for each locus.
#' @source https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments
#' @note For internal use only.
"AA_atlas"
