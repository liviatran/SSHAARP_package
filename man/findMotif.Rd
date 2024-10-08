% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findMotif.R
\name{findMotif}
\alias{findMotif}
\title{Returns an alignment data frame of alleles that share a specific amino acid motif}
\usage{
findMotif(motif, filename)
}
\arguments{
\item{motif}{An amino acid motif in the following format: Locus*##$~##$~##$, where ## identifies a peptide position, and $ identifies an amino acid residue. Motifs can include any number of amino acids.}

\item{filename}{The full file path of the user specified dataset if the user wishes to use their own file, or the pre-bundled Solberg dataset. User provided datasets must be a .dat, .txt, or.csv file, and must conform to the structure and format of the Solberg dataset.}
}
\value{
An amino acid alignment dataframe of alleles that share the specified motif. Otherwise, a vector containing FALSE and an error message is returned.
}
\description{
Returns an alignment data frame of alleles that share a specific amino acid motif.
}
\examples{

#example with existing motif
\dontrun{findMotif("DRB1*26F~28E~30Y", filename=solberg_dataset)}
\dontrun{findMotif("DRB1*26F~28E", filename=solberg_dataset)}

#example with non-existent motif
\dontrun{findMotif("DRB1*26F~28E~30Z", filename=solberg_dataset)}

#extracting names of alleles with user-defined motif
findMotif("DRB1*26F~28E~30Y", filename=solberg_dataset)[,4]
}
