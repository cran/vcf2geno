#' Read a gene from VCF file and return a genotypes matrix
#'
#' @param fileName input VCF file (Bgzipped, with Tabix index)
#' @param geneFile a text file listing all genes in refFlat format
#' @param geneName which gene(s) to be extracted
#' @param annoType which annotation you would like to extract
#' @return genotype matrix
readVCFToMatrixByGene <- function(fileName, geneFile, geneName, annoType) {
  stopifnot(file.exists(fileName))

  storage.mode(fileName) <- "character"
  storage.mode(geneFile) <- "character"
  storage.mode(geneName) <- "character"
  storage.mode(annoType) <- "character"
  .Call("readVCFToMatrixByGene", fileName, geneFile, geneName, annoType, PACKAGE="vcf2geno");
};

#' Read a gene from VCF file and return a genotypes matrix
#'
#' @param fileName input VCF file (Bgzipped, with Tabix index)
#' @param range a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @param annoType which annotation you would like to extract
#' @return genotype matrix
readVCFToMatrixByRange <- function(fileName, range, annoType) {
  stopifnot(file.exists(fileName))
  storage.mode(fileName) <- "character"
  storage.mode(range)    <- "character"
  storage.mode(annoType) <- "character"
  .Call("readVCFToMatrixByRange", fileName, range, annoType, PACKAGE="vcf2geno");
};

#' Read a gene from VCF file and return a genotypes matrix
#'
#' @param fileName input VCF file (Bgzipped, with Tabix index)
#' @param geneFile a text file listing all genes in refFlat format
#' @param geneName which gene(s) to be extracted
#' @param annoType which annotation you would like to extract
#' @param vcfColumn which vcf columns to extract, e.g. CHROM, POS, FILTER...
#' @param vcfInfo which tag in the INFO columns to extarct. e.g. DP, AC...
#' @param vcfIndv which individual tag to extract, e.g. GT, GQ...
#' @return a list of genes, and each elements has specified vcfColumn, vcfinfo, vcfIndv
readVCFToListByGene <- function(fileName, geneFile, geneName, annoType, vcfColumn, vcfInfo, vcfIndv) {
  stopifnot(file.exists(fileName))
  storage.mode(fileName) <- "character"
  storage.mode(geneFile) <- "character"
  storage.mode(geneName) <- "character"
  storage.mode(annoType) <- "character"
  storage.mode(vcfColumn)<- "character"
  storage.mode(vcfInfo)  <- "character"
  storage.mode(vcfIndv)  <- "character"
  .Call("readVCFToListByGene", fileName, geneFile, geneName, annoType, vcfColumn, vcfInfo, vcfIndv, PACKAGE="vcf2geno");
};

#' Read a gene from VCF file and return a genotypes matrix
#'
#' @param fileName input VCF file (Bgzipped, with Tabix index)
#' @param range a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @param annoType which annotation you would like to extract
#' @param vcfColumn which vcf columns to extract, e.g. CHROM, POS, FILTER...
#' @param vcfInfo which tag in the INFO columns to extarct. e.g. DP, AC...
#' @param vcfIndv which individual tag to extract, e.g. GT, GQ...
#' @return a list of genes, and each elements has specified vcfColumn, vcfinfo, vcfIndv
readVCFToListByRange <- function(fileName, range, annoType, vcfColumn, vcfInfo, vcfIndv) {
  stopifnot(file.exists(fileName))
  storage.mode(fileName) <- "character"
  storage.mode(range)    <- "character"
  storage.mode(annoType) <- "character"
  storage.mode(vcfColumn)<- "character"
  storage.mode(vcfInfo)  <- "character"
  storage.mode(vcfIndv)  <- "character"
  .Call("readVCFToListByRange", fileName, range, annoType, vcfColumn, vcfInfo, vcfIndv, PACKAGE="vcf2geno");
};

#' Read a gene from VCF file and return a genotypes matrix
#'
#' @param scoreTestFiles score test output file (rvtests outputs)
#' @param covFiles covaraite files (vcf2ld_gene outputs)
#' @param geneFile speicify which gene file to extract (need to be refFlat format)
#' @param gene speicify which gene to extract
#' @return a list of genes, and each elements contain genotype covariance within gene and associated score test statsitics. 
rvmeta.readData <- function(scoreTestFiles, covFiles, geneFile, gene) {
  stopifnot(file.exists(scoreTestFiles))
  stopifnot(file.exists(covFiles))  
  storage.mode(scoreTestFiles) <- "character"
  storage.mode(covFiles) <- "character"
  storage.mode(geneFile) <- "character"
  storage.mode(gene) <- "character"
  .Call("rvMetaReadData", scoreTestFiles, covFiles, geneFile, gene, PACKAGE="vcf2geno");
};

## rvmeta.CMH <- function(score.stat.vec.list, maf.vec.list, cov.mat.list, var.Y.list, N.list, alternative=c('two.sided','greater','less'), no.boot, rv.test, extra.pars=list())  {
  
## };

#' Read covariance
#'
#' @param covFile covaraite file (vcf2ld_window outputs)
#' @param tabixRange a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @return a matrix of covariance within gene
rvmeta.readCovByRange <- function(covFile, tabixRange) {
  stopifnot(file.exists(covFile))    
  storage.mode(covFile) <- "character"
  storage.mode(tabixRange) <- "character"
  .Call("readCovByRange", covFile, tabixRange, PACKAGE="vcf2geno");
};

#' Read score test statistics
#'
#' @param scoreFile score test output file (rvtests outputs)
#' @param tabixRange a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @return score test statistics in given range
rvmeta.readScoreByRange <- function(scoreFile, tabixRange) {
  stopifnot(file.exists(scoreFile))
  storage.mode(scoreFile) <- "character"
  storage.mode(tabixRange) <- "character"
  .Call("readScoreByRange", scoreFile, tabixRange, PACKAGE="vcf2geno");
};
