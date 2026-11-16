\name{rings_chordless}
\alias{rings_chordless}
\alias{rings_chordless.dseams_frame}
\alias{seams_rings_chordless}
\title{Chordless cycles of one size}
\usage{
rings_chordless(frame, size = 6L, cutoff = 3.5, type = 0L, graph = "cutoff")

seams_rings_chordless(sys, cutoff = 3.5, size = 6L, type = 0L, graph = 0L)
}
\arguments{
\item{frame}{A \code{dseams_frame}.}
\item{sys}{A handle from \code{seams_system_new}.}
\item{size}{Cycle length, from 3 through 12.}
\item{cutoff}{Neighbour cutoff in angstrom.}
\item{type}{LAMMPS type. 0 bonds every atom.}
\item{graph}{\code{"cutoff"} or \code{"mutual"} for \code{rings_chordless};
0 (cutoff) or 1 (mutual) for \code{seams_rings_chordless}.}
}
\value{
\code{rings_chordless} returns an integer matrix with one row per cycle and
\code{size} columns of 1-based atom rows. The \code{per_atom} attribute
counts the cycles through each atom.

\code{seams_rings_chordless} returns a list with \code{rings}, the same
matrix with 0-based indices as \code{seams_rings_chordless} in C writes
them, and \code{per_atom}.
}
\description{
A cycle is chordless when no bond joins two of its atoms that are not next
to each other on it. Every primitive ring is chordless; a chordless cycle
need not be primitive. Each cycle starts at its smallest atom and runs
towards the smaller of that atom's two cycle neighbours, and the rows are
sorted.
}
