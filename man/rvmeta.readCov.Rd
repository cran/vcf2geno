\name{rvmeta.readCov}
\alias{rvmeta.readCov}
\title{Read covariance}
\usage{
  rvmeta.readCov(covFile, tabixRange)
}
\arguments{
  \item{covFile}{covaraite file (vcf2ld_window outputs)}

  \item{tabixRange}{a text indicating which range in the
  VCF file to extract. e.g. 1:100-200}
}
\value{
  a matrix of covariance within gene
}
\description{
  Read covariance
}

