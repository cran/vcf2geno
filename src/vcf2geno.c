#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>

#if 0
SEXP out(SEXP x, SEXP y) {
  R_len_t i, j, nx = length(x), ny = length(y);
  double tmp, *rx = REAL(x), *ry = REAL(y), *rans;
  SEXP ans, dim, dimnames;
  PROTECT(ans = allocMatrix(REALSXP, nx, ny));
  rans = REAL(ans);
  for(i = 0; i < nx; i++) {
    tmp = rx[i];
    for(j = 0; j <ny ; j++ )
      rans[i + nx *j ] = tmp * ry[j];
  }

  PROTECT(dim = allocVector(INTSXP, 2));
  INTEGER(dim)[0] = nx; INTEGER(dim)[1] = ny;
  setAttrib(ans, R_DimSymbol, dim);
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, getAttrib(x, R_NamesSymbol));
  SET_VECTOR_ELT(dimnames, 1, getAttrib(y, R_NamesSymbol));
  setAttrib(ans, R_DimNamesSymbol, dimnames);
  UNPROTECT(3);
  return(ans);
};
#endif

extern SEXP impl_readVCFToMatrixByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType);
SEXP readVCFToMatrixByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType) {
  return impl_readVCFToMatrixByGene(arg_fileName, arg_geneFile, arg_geneName, arg_annoType);
}

extern SEXP impl_readVCFToMatrixByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType);
SEXP readVCFToMatrixByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType) {
  return impl_readVCFToMatrixByRange(arg_fileName, arg_range, arg_annoType);
}


extern SEXP impl_readVCFToListByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag);
SEXP readVCFToListByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag) {
  return impl_readVCFToListByGene(arg_fileName, arg_geneFile, arg_geneName, arg_annoType, arg_columns, arg_infoTag, arg_indvTag);
}

extern SEXP impl_readVCFToListByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag);
SEXP readVCFToListByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag) {
  return impl_readVCFToListByRange(arg_fileName, arg_range, arg_annoType, arg_columns, arg_infoTag, arg_indvTag);
}

extern SEXP impl_rvMetaReadData(SEXP arg_pvalFile, SEXP arg_covFile, SEXP arg_gene);
SEXP rvMetaReadData(SEXP arg_pvalFile, SEXP arg_covFile, SEXP arg_gene) {
  return impl_rvMetaReadData(arg_pvalFile, arg_covFile, arg_gene);
}

extern SEXP impl_readCovByRange(SEXP arg_covFile, SEXP arg_range);
SEXP readCovByRange(SEXP arg_covFile, SEXP arg_range) {
  return impl_readCovByRange(arg_covFile, arg_range);
}

extern SEXP impl_readScoreByRange(SEXP arg_covFile, SEXP arg_range);
SEXP readScoreByRange(SEXP arg_covFile, SEXP arg_range) {
  return impl_readScoreByRange(arg_covFile, arg_range);
}

