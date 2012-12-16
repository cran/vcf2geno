#include <iostream> // have to include this to avoid R screw up "length" in iostream/sstream...
#include <R.h>
#include <Rinternals.h>

extern "C" SEXP impl_readVCFToMatrixByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType);
extern "C" SEXP impl_readVCFToMatrixByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType);

extern "C" SEXP impl_readVCFToListByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag);
extern "C" SEXP impl_readVCFToListByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag);

extern "C" SEXP impl_rvMetaReadData(SEXP arg_pvalFile, SEXP arg_covFile, SEXP arg_geneFile, SEXP arg_gene);
extern "C" SEXP impl_readCovByRange(SEXP arg_covFile, SEXP arg_range);
extern "C" SEXP impl_readScoreByRange(SEXP arg_scoreFile, SEXP arg_range);
