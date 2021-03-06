% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distances.R
\name{res.dist}
\alias{res.dist}
\title{Compute Distances Between Residues}
\usage{
res.dist(pdb,  rA, chainA, rB, chainB, backbone = FALSE, hatoms = FALSE)
}
\arguments{
\item{pdb}{is either a PDB id, or the path to a pdb file.}

\item{rA}{an integer indicating the position of the first residue.}

\item{chainA}{a character indicating the chain to which belong the first residue.}

\item{rB}{an integer indicating the position of the second residue.}

\item{chainB}{a character indicating the chain to which belong the second residue.}

\item{backbone}{logical, when TRUE it means that we include those atoms belonging to the main chain (CA, N, O and C) beside all the side chain atoms.}

\item{hatoms}{logical, if TRUE we include all the hydrogen atoms in the computation as long as the PDB provides their coordinates.}
}
\value{
This function returns a list of three elements, where each of these elements is, in turn, a list of three elements providing information regarding minimal, maximal and averaged distances.
}
\description{
Computes the euclidean distance between two given residues
}
\examples{
res.dist('1q8k', 51, 'A', 55, 'A', backbone = TRUE, hatoms = TRUE)
}
\seealso{
pairwise.dist(), dist2closest(), ball()
}
\author{
Juan Carlos Aledo
}
