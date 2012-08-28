\name{vcf2geno-package}
\alias{vcf2geno-package}
\alias{vcf2geno}
\docType{package}
\title{
  Vcf2geno make reading VCF easy
}
\description{
  Help read VCF files into R
}
\details{
\tabular{ll}{
Package: \tab vcf2geno\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-05-24\cr
License: \tab BSD \cr
}

~~ An overview of how to use the package, including the most important ~~
You will need tabix indexed vcf files. Then use readVCF series functions to read.

~~ functions ~~
readVCFToMatrixByGene
readVCFToMatrixByRange
readVCFToListByGene
readVCFToListByRange
rvmeta.readData
rvmeta.readCov
rvmeta.readScoreByRange

}
\author{
Xiaowei Zhan <zhanxw@gmail.com>
}
\keyword{ vcf }
\examples{
##readVCFToListByGene("input.vcf.gz", "refFlatGeneFile.txt.gz", "CFH", "Synonymous", c("CHROM","POS"), c("DP"), c("GT")

}
