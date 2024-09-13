## Mock Haplotype Dataset 15DEC21
#' Mock Haplotype Dataset
#'
#' A dataframe containing mock haplotype data modeled after the Allele Frequency Network Database (AFND) haplotype data structure.
#' @docType data
#' @name mock_haplotype_dataset
#' @usage mock_haplotype_dataset
"mock_haplotype_dataset"


## Solberg Dataset 18AUG2024
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

## IMGTprotalignments v 2.0.0 19AUG2024
#' Protein alignments for all protein coding genes in the IPD-IMGT/HLA Database release v 3.57.0
#'
#' A list object containing protein alignments for all protein coding genes in the IPD-IMGT/HLA Database release. Alignments were downloaded from the ANHIG/IMGTHLA Github repository.
#' @docType data
#' @name IMGTprotalignments
#' @usage IMGTprotalignments
#' @format A list containing protein alignments in a dataframe format for each locus.
#' @source https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments
#' @source Robinson, J, Barker, DJ, Georgiou, X, Cooper, MA, Flicek, P, Marsh, SGE
#'         The IPD-IMGT/HLA Database
#'         Nucleic Acids Research (2020) 43:D948-D955
"IMGTprotalignments"
