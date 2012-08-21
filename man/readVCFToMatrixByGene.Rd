\name{readVCFToMatrixByGene}
\alias{readVCFToMatrixByGene}
\title{Read a gene from VCF file and return a genotypes matrix}
\usage{
  readVCFToMatrixByGene(fileName, geneFile, geneName,
    annoType)
}
\arguments{
  \item{fileName}{input VCF file}

  \item{geneFile}{a text file listing all genes in refFlat
  format}

  \item{geneName}{which gene(s) to be extracted}

  \item{annoType}{which annotation you would like to
  extract}
}
\value{
  genotype matrix
}
\description{
  Read a gene from VCF file and return a genotypes matrix
}

