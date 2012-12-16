\name{readVCFToListByGene}
\alias{readVCFToListByGene}
\title{Read a gene from VCF file and return a genotypes matrix}
\usage{
  readVCFToListByGene(fileName, geneFile, geneName,
    annoType, vcfColumn, vcfInfo, vcfIndv)
}
\arguments{
  \item{fileName}{input VCF file (Bgzipped, with Tabix
  index)}

  \item{geneFile}{a text file listing all genes in refFlat
  format}

  \item{geneName}{which gene(s) to be extracted}

  \item{annoType}{which annotation you would like to
  extract}

  \item{vcfColumn}{which vcf columns to extract, e.g.
  CHROM, POS, FILTER...}

  \item{vcfInfo}{which tag in the INFO columns to extarct.
  e.g. DP, AC...}

  \item{vcfIndv}{which individual tag to extract, e.g. GT,
  GQ...}
}
\value{
  a list of genes, and each elements has specified
  vcfColumn, vcfinfo, vcfIndv
}
\description{
  Read a gene from VCF file and return a genotypes matrix
}

