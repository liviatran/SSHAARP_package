% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataSubset.R
\name{dataSubset}
\alias{dataSubset}
\title{Solberg dataset manipulation}
\usage{
dataSubset(motif, filename = SSHAARP::solberg_dataset)
}
\arguments{
\item{motif}{An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.}

\item{filename}{The filename of the local copy of the Solberg dataset - the defaulted filename is the solberg_dataset in the SSHAARP package.}
}
\value{
A data frame containing a reformatted version of the Solberg dataset, with rows ordered by population name, Cartesian coordinates in the latit and longit columns, and limited to populations with data for the specified locus. If a motif has formatting errors, a warning message is returned.
}
\description{
Returns a modified version of the Solberg dataset that includes a column of locus*allele names, is sorted by by population name, and is reduced to the specified locus. Cardinal coordinates are converted to their Cartesian equivalents (i.e. 50S is converted to -50).
}
\note{
For internal SSHAARP use only.

The Solberg dataset is the tab-delimited ‘1-locus-alleles.dat’ text file in the results.zip archive at http://pypop.org/popdata/.

The Solberg dataset is also prepackaged into SSHAARP as 'solberg_dataset'.
}
