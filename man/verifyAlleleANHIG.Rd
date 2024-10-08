% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/verfiyAlleleANHIG.R
\name{verifyAlleleANHIG}
\alias{verifyAlleleANHIG}
\title{Verifies the allele entered is present in IMGTprotalignments}
\usage{
verifyAlleleANHIG(allele, filename)
}
\arguments{
\item{allele}{An allele name written in the IPD-IMGT/HLA Database format.}

\item{filename}{The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the Solberg dataset.}
}
\value{
TRUE if allele is present in the IMGTprotalignment object. Otherwise, a vector containing FALSE and an error message is returned.
}
\description{
Verifies the allele entered is present in IMGTprotalignments.
}
\note{
For internal SSHAARP use only.
}
\examples{
#Example of an allele that exists in the IMGTprotalignment object
\dontrun{verifyAlleleANHIG("B*07:02:01:01", filename="your_haplotype_dataset")}

#Example of an allele that does not exist in the IMGTprotalignment object
verifyAlleleANHIG("B*01:01:01", filename=solberg_dataset)
}
