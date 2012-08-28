\name{readVCFToListByRange}
\alias{readVCFToListByRange}
\title{Read a gene from VCF file and return a genotypes matrix}
\usage{
  readVCFToListByRange(fileName, range, annoType,
    vcfColumn, vcfInfo, vcfIndv)
}
\arguments{
  \item{fileName}{input VCF file}

  \item{range}{a text indicating which range in the VCF
  file to extract. e.g. 1:100-200}

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

