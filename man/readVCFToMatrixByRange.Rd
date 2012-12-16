\name{readVCFToMatrixByRange}
\alias{readVCFToMatrixByRange}
\title{Read a gene from VCF file and return a genotypes matrix}
\usage{
  readVCFToMatrixByRange(fileName, range, annoType)
}
\arguments{
  \item{fileName}{input VCF file (Bgzipped, with Tabix
  index)}

  \item{range}{a text indicating which range in the VCF
  file to extract. e.g. 1:100-200}

  \item{annoType}{which annotation you would like to
  extract}
}
\value{
  genotype matrix
}
\description{
  Read a gene from VCF file and return a genotypes matrix
}

