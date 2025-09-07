\name{seqi1epi}
\alias{seqi1epi}

\title{
At least one episode in the states
}
\description{
  Returns whether each sequence comprises at least one episode in the states.
  }
  
\usage{
seqi1epi(seqdata)
}

\arguments{
  \item{seqdata}{a sequence object (see \code{\link[TraMineR]{seqdef}} function).}
}

\references{
Gabadinho, A., G. Ritschard, N. S. Müller and M. Studer (2011). Analyzing and Visualizing State Sequences in R with TraMineR. \emph{Journal of Statistical Software} \bold{40}(4), 1-37.
}

\author{Nicolas Robette}

\seealso{
\code{\link[TraMineR]{seqistatd}}, \code{\link{seqinepi}}, \code{\link{seqifpos}}
}

\examples{
data(trajact)
seqact <- seqdef(trajact)
stat <- seqi1epi(seqact)
head(stat)
}

\keyword{Longitudinal characteristics}
\keyword{State sequences} 
