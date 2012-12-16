\name{rvmeta.readData}
\alias{rvmeta.readData}
\title{Read a gene from VCF file and return a genotypes matrix}
\usage{
  rvmeta.readData(scoreTestFiles, covFiles, geneFile, gene)
}
\arguments{
  \item{scoreTestFiles}{score test output file (rvtests
  outputs)}

  \item{covFiles}{covaraite files (vcf2ld_gene outputs)}

  \item{geneFile}{speicify which gene file to extract (need
  to be refFlat format)}

  \item{gene}{speicify which gene to extract}
}
\value{
  a list of genes, and each elements contain genotype
  covariance within gene and associated score test
  statsitics.
}
\description{
  Read a gene from VCF file and return a genotypes matrix
}

