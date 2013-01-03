\name{rvmeta.readDataByRange}
\alias{rvmeta.readDataByRange}
\title{Read a range from VCF file and return a genotypes matrix}
\usage{
  rvmeta.readDataByRange(scoreTestFiles, covFiles, ranges)
}
\arguments{
  \item{scoreTestFiles}{score test output file (rvtests
  outputs)}

  \item{covFiles}{covaraite files (vcf2ld_gene outputs)}

  \item{ranges}{speicify which range to extract (e.g.
  1:100-200)}
}
\value{
  a list of variants, and each elements contain genotype
  covariance within gene and associated score test
  statsitics.
}
\description{
  Read a range from VCF file and return a genotypes matrix
}

