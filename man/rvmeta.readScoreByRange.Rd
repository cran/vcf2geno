\name{rvmeta.readScoreByRange}
\alias{rvmeta.readScoreByRange}
\title{Read score test statistics}
\usage{
  rvmeta.readScoreByRange(scoreFile, tabixRange)
}
\arguments{
  \item{scoreFile}{score test output file (rvtests
  outputs)}

  \item{tabixRange}{a text indicating which range in the
  VCF file to extract. e.g. 1:100-200}
}
\value{
  score test statistics in given range
}
\description{
  Read score test statistics
}

