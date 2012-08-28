#include "vcf2genoLoader.h"

#include <string>
#include <map>
#include <vector>
#include <set>

#include "VCFUtil.h"
#include "tabix.h"
#include "R_CPP_interface.h"

#include "R.h"

void loadGeneFile(const std::string& geneFile, const std::string& geneName, std::map< std::string, std::string>* geneRangePtr){
  std::map< std::string, std::string>& geneRange = *geneRangePtr;

  // load gene ranges
  if (geneName.size() ){
    if (geneFile.size() == 0) {
      //REprintf("Have to provide --geneFile to extract by gene.\n");
      //REprintf("Critical error happening!"); //abort();
      error("Have to provide --geneFile to extract by gene.\n");
    }
    LineReader lr(geneFile);
    std::vector< std::string > fd;
    while (lr.readLineBySep(&fd, "\t ")) {
      if (geneName != fd[0]) continue;
      fd[2] = chopChr(fd[2]); // chop "chr1" to "1"
      if (geneRange.find(fd[0])  == geneRange.end()) {
        geneRange[fd[0]] = fd[2] + ":" + fd[4] + "-" + fd[5];
      } else {
        geneRange[fd[0]] += "," + fd[2] + ":" + fd[4] + "-" + fd[5];
      }
    };
  }
}

/**
 * Read from @param vin and return a matrix of marker by people
 */
SEXP readVCF2Matrix(VCFExtractor* vin) {
  std::vector<double> genoVec;
  std::vector<std::string> posVec;
  std::vector<std::string> idVec;
  std::string posString;

  // print header
  std::vector<std::string>& names = idVec;
  vin->getVCFHeader()->getPeopleName(&names);

  while (vin->readRecord()){
    // REprintf("read a record\n");
    VCFRecord& r = vin->getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;

    // store all results here
    posString = r.getChrom();
    posString += ':';
    posString += r.getPosStr();
    posVec.push_back(posString);

    for (size_t i = 0; i < people.size(); i++) {
      indv = people[i];
      int g = indv->justGet(0).getGenotype();
      //fprintf(stdout, "\t%d", g);
      genoVec.push_back(g);
    }
    //fprintf(stdout, "\n");
  }; // end while

  //  REprintf("posVec = %zu, idVec = %zu, genoVec = %zu\n", posVec.size(), idVec.size(), genoVec.size());

  // pass value back to R (see Manual Chapter 5)


  int nx = (int) posVec.size();
  int ny = (int) idVec.size();

  SEXP ans = R_NilValue;

  PROTECT(ans = allocMatrix(REALSXP, nx, ny));
  double* rans = REAL(ans);
  int idx = 0;
  for(int i = 0; i < nx; i++) {
    for(int j = 0; j <ny ; j++ ) {
      // Rprintf("idx = %d, i = %d, j=%d, geno = %g\n", idx, i, j, genoVec[idx]);
      rans[i + nx * j] = genoVec[idx];
      ++idx;
    }
  }

  // set row and col names
  SEXP dim;
  PROTECT(dim = allocVector(INTSXP, 2));
  INTEGER(dim)[0] = nx; INTEGER(dim)[1] = ny;
  setAttrib(ans, R_DimSymbol, dim);

  SEXP rowName;
  PROTECT(rowName=allocVector(STRSXP, nx));
  for (int i = 0; i < nx; i++ )
    SET_STRING_ELT(rowName, i, mkChar(posVec[i].c_str()));
  SEXP colName;
  PROTECT(colName=allocVector(STRSXP, ny));
  for (int i = 0; i < ny; i++ )
    SET_STRING_ELT(colName, i, mkChar(idVec[i].c_str()));

  SEXP dimnames;
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, rowName);
  SET_VECTOR_ELT(dimnames, 1, colName);
  setAttrib(ans, R_DimNamesSymbol, dimnames);

  // finish up
  UNPROTECT(5);
  return(ans);

}


/**
 * @param arg_fileName: a string character
 * @param arg_geneFile: which gene file to use
 * @param arg_geneName: which gene we are interested. (just allow One gene name).
 * @param arg_annoType: allow annotation type, can be regular expression. (e.g. Synonymous|Nonsynonymous)
 */
SEXP impl_readVCFToMatrixByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType) {
  SEXP ans = R_NilValue;

  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName,0));
  std::string FLAG_range = CHAR(STRING_ELT(arg_range,0));
  std::string FLAG_annoType = CHAR(STRING_ELT(arg_annoType,0));

  if (FLAG_fileName.size() == 0) {
    error("Please provide VCF file name");
    return ans;
  }
  if (FLAG_range.size() == 0) {
    error("Please provide a given range, e.g. '1:100-200'");
    return ans;
  }

  // REprintf("range = %s\n", FLAG_range.c_str());
  VCFExtractor vin(FLAG_fileName.c_str());
  vin.setRangeList(FLAG_range.c_str());

  if (FLAG_annoType.size()) {
    vin.setAnnoType(FLAG_annoType.c_str());
  }

  // real working part
  return readVCF2Matrix(&vin);
}

/**
 * @param arg_fileName: a string character
 * @param arg_geneFile: which gene file to use
 * @param arg_geneName: which gene we are interested. (just allow One gene name).
 * @param arg_annoType: allow annotation type, can be regular expression. (e.g. Synonymous|Nonsynonymous)
 */
SEXP impl_readVCFToMatrixByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType) {
  SEXP ans = R_NilValue;

  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName,0));
  std::string FLAG_geneFile = CHAR(STRING_ELT(arg_geneFile,0));
  std::string FLAG_geneName = CHAR(STRING_ELT(arg_geneName,0));
  std::string FLAG_annoType = CHAR(STRING_ELT(arg_annoType,0));

  if (FLAG_fileName.size() == 0) {
    error("Please provide VCF file name");
  }
  if (FLAG_geneName.size() && FLAG_geneFile.size() == 0) {
    error("Please provide gene file name when extract genotype by gene");
  }
#if 0
  // default values to easy testing
  if (!FLAG_fileName.size()) {
    FLAG_fileName = "/net/fantasia/home/zhanxw/rvtests/executable/all.anno.filtered.vcf.gz";
    warning("Use example indexed annotated VCF file: /net/fantasia/home/zhanxw/rvtests/executable/all.anno.filtered.vcf.gz");
  }
  if (FLAG_geneName.size()) {
    if (!FLAG_geneFile.size()) {
      FLAG_geneFile = "/net/fantasia/home/zhanxw/rvtests/executable/refFlat_hg19_uniq_gene.txt.gz";
      warning(" gene files: /net/fantasia/home/zhanxw/rvtests/executable/refFlat_hg19_uniq_gene.txt.gz");
    }
  } else {
    FLAG_geneName = "CFH";
    warning("We set gene name to be CFH, otherwise you may dump the whole VCF file (huge memory!)");
  }
#endif

  std::map< std::string, std::string> geneRange;
  loadGeneFile(FLAG_geneFile, FLAG_geneName, &geneRange);
  std::string range;
  for (std::map< std::string, std::string>::iterator it = geneRange.begin();
       it != geneRange.end();
       it++) {
    if (range.size() > 0) {
      range += ",";
    }
    range += it->second;
  };

  //fprintf(stdout, "range = %s\n", range.c_str());
  VCFExtractor vin(FLAG_fileName.c_str());
  if (range.size())
    vin.setRangeList(range.c_str());
  else {
    warning("Gene name [ %s ] does not exists in provided gene file", FLAG_geneName.c_str());
    return (ans);
  };

  if (FLAG_annoType.size()) {
    vin.setAnnoType(FLAG_annoType.c_str());
  }

  // real working part
  return readVCF2Matrix(&vin);
}


SEXP readVCF2List(VCFInputFile* vin,
                  const std::set<std::string>& FLAG_vcfColumn,
                  const std::vector<std::string>& FLAG_infoTag,
                  const std::vector<std::string>& FLAG_indvTag) {


  // Rprintf("vcfColumn.size() = %u\n", FLAG_vcfColumn.size());
  // Rprintf("vcfInfo.size() = %u\n", FLAG_infoTag.size());
  // Rprintf("vcfIndv.size() = %u\n", FLAG_indvTag.size());
  // also append sample names at the end
  int retListLen = FLAG_vcfColumn.size() + FLAG_infoTag.size() + FLAG_indvTag.size() + 1;
  if (retListLen == 0) {
    return R_NilValue;
  }

  int numAllocated = 0; // record how many times we allocate (using PROTECT in R);
  SEXP ret;
  PROTECT(ret = allocVector(VECSXP, retListLen));
  numAllocated ++;

  //  store results
  std::vector<std::string> idVec;
  std::vector<std::string> chrom;
  std::vector<int> pos;
  std::vector<std::string> rsId;
  std::vector<std::string> ref;
  std::vector<std::string> alt;
  std::vector<std::string> qual;
  std::vector<std::string> filt;
  std::vector<std::string> info;
  std::vector<std::string> format;

  std::map<std::string, std::vector<std::string> > infoMap;

  // std::vector<int> gtVec;
  // std::vector<int> gdVec;
  // std::vector<int> gqVec;

  std::map<std::string, std::vector<std::string> > indvMap;
  int nRow = 0; // # of positions that will be outputed

  // print header
  std::vector<std::string>& names = idVec;
  vin->getVCFHeader()->getPeopleName(&names);


  bool FLAG_variantOnly = false;
  // real working part
  bool nonVariantSite;
  while (vin->readRecord()){
    // REprintf("read a record\n");
    VCFRecord& r = vin->getVCFRecord();
    VCFPeople& people = r.getPeople();
    VCFIndividual* indv;
    if (FLAG_variantOnly) {
      // REprintf("filter by var\n");
      bool hasVariant = false;
      int geno;
      int GTidx = r.getFormatIndex("GT");
      for (size_t i = 0; i < people.size() ;i ++) {
        indv = people[i];
        geno = indv->justGet(GTidx).getGenotype();
        if (geno != 0 && geno != MISSING_GENOTYPE)
          hasVariant = true;
      }
      if (!hasVariant) {
        nonVariantSite++;
        continue;
      }
    }

    // store results here
    nRow++;
    if (FLAG_vcfColumn.count("CHROM")){
      chrom.push_back(r.getChrom());
    }
    if (FLAG_vcfColumn.count("POS")){
      pos.push_back(r.getPos());
    }
    if (FLAG_vcfColumn.count("ID")){
      rsId.push_back(r.getID());
    }
    if (FLAG_vcfColumn.count("REF")){
      ref.push_back(r.getRef());
    }
    if (FLAG_vcfColumn.count("ALT")){
      alt.push_back(r.getAlt());
    }
    if (FLAG_vcfColumn.count("QUAL")){
      qual.push_back(r.getQual());
    }
    if (FLAG_vcfColumn.count("FILTER")){
      filt.push_back(r.getFilt());
    }
    if (FLAG_vcfColumn.count("INFO")){
      info.push_back(r.getInfo());
    }
    if (FLAG_vcfColumn.count("FORMAT")){
      format.push_back(r.getFormat());
    }

    // store INFO field
    for (std::vector<std::string>::const_iterator it = FLAG_infoTag.begin(); it != FLAG_infoTag.end(); ++it) {
      bool missing;
      VCFValue v = r.getInfoTag(it->c_str(), &missing);
      if (missing) {
        infoMap[ *it ].push_back("");
      } else {
        infoMap[ *it ].push_back(v.toStr());
        // Rprintf("add info field [ %s ] = %s\n", it->c_str(), v.toStr());
      }
    };
    // Rprintf("Done add info\n");

    // store indv values
    for (size_t i = 0; i < people.size(); i++) {
      indv = people[i];

      for (std::vector<std::string>::const_iterator it = FLAG_indvTag.begin(); it != FLAG_indvTag.end(); ++it) {
        int idx = r.getFormatIndex(it->c_str());
        if (idx < 0) {
          indvMap[ *it ].push_back("");
        } else {
          bool missing;
          VCFValue v = indv->get(idx, &missing);
          if (missing) {
            indvMap[ *it ].push_back("");
          } else{
            indvMap[ *it ].push_back(v.toStr());
            // Rprintf("add indv field [ %s ] = %s\n", it->c_str(), v.toStr());
          }
        }
      };
    }
    // Rprintf("Done add indv\n");
  }; // end while
//   Rprintf("indvMap.size() = %zu\n", indvMap.size());
  // REprintf("posVec = %zu, idVec = %zu, genoVec = %zu\n", posVec.size(), idVec.size(), genoVec.size());

  // pass value back to R (see Manual Chapter 5)
  std::vector<std::string> listNames;
  int retListIdx = 0;
  if (FLAG_vcfColumn.count("CHROM")) {
    numAllocated += storeResult(chrom, ret, retListIdx++);
    listNames.push_back("CHROM");
  }
  if (FLAG_vcfColumn.count("POS")) {
    numAllocated += storeResult(pos, ret, retListIdx++);
    listNames.push_back("POS");
  }
  if (FLAG_vcfColumn.count("ID"))  {
    numAllocated += storeResult(rsId, ret, retListIdx++);
    listNames.push_back("ID");
  }
  if (FLAG_vcfColumn.count("REF")) {
    numAllocated += storeResult(ref, ret, retListIdx++);
    listNames.push_back("REF");
  }
  if (FLAG_vcfColumn.count("ALT")) {
    numAllocated += storeResult(alt, ret, retListIdx++);
    listNames.push_back("ALT");
  }
  if (FLAG_vcfColumn.count("QUAL")) {
    numAllocated += storeResult(qual, ret, retListIdx++);
    listNames.push_back("QUAL");
  }
  if (FLAG_vcfColumn.count("FILTER")) {
    numAllocated += storeResult(filt, ret, retListIdx++);
    listNames.push_back("FILTER");
  }
  if (FLAG_vcfColumn.count("INFO")) {
    numAllocated += storeResult(info, ret, retListIdx++);
    listNames.push_back("INFO");
  }
  if (FLAG_vcfColumn.count("FORMAT")) {
    numAllocated += storeResult(format, ret, retListIdx++);
    listNames.push_back("FORMAT");
  }
  // pass info values to R
  for ( std::map<std::string, std::vector<std::string> >::iterator it = infoMap.begin();
        it != infoMap.end();
        ++it) {
    numAllocated += storeResult(it->first, it->second, ret, retListIdx++);
    listNames.push_back(it->first);
  }
  // pass indv tags to R
  // Rprintf("pass idnv tags\n");
  for ( std::map<std::string, std::vector<std::string> >::iterator it = indvMap.begin();
        it != indvMap.end();
        ++it) {

    dump(it->second);
    numAllocated += storeResult(it->first, it->second, ret, retListIdx);
    // Rprintf("results done\n");
    // NOTE: R internally store values into matrix by column first!
    // thus the matrix is people by marker
    numAllocated += setDim(idVec.size(), nRow, ret, retListIdx);
    retListIdx ++;
    listNames.push_back(it->first);
  }
  // Rprintf("pass idnv tags done.\n");

  // store sample ids
  // Rprintf("set sample id");
  listNames.push_back("sampleId");
  numAllocated += storeResult(idVec, ret, retListIdx++);

  // Rprintf("set list names\n");
  SEXP sListNames;
  PROTECT(sListNames = allocVector(STRSXP, listNames.size()));
  numAllocated ++;
  for (unsigned int i = 0; i != listNames.size(); ++i){
    SET_STRING_ELT(sListNames, i, mkChar(listNames[i].c_str()));
  }
  setAttrib(ret, R_NamesSymbol, sListNames);

  // finish up
  UNPROTECT(numAllocated);
  // Rprintf("Unprotected: %d\n", (retListLen + 1));
  return(ret);

}

/**
 * @param arg_fileName: a string character
 * @param arg_range: which range to extract
 * @param arg_annoType: allow annotation type, can be regular expression. (e.g. Synonymous|Nonsynonymous)
 * @param arg_columns: a list of which columns to extract (e.g. CHROM, POS ...)
 * @param arg_infoTag: a list of which tag under INFO tag will be extracted (e.g. ANNO, ANNO_FULL, AC ...)
 * @param arg_indvTag: a list of which tag given in individual's column (e.g. GT, GD, GQ ...)
 */
SEXP impl_readVCFToListByRange(SEXP arg_fileName, SEXP arg_range, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag){
  // begin
  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName,0));
  std::string FLAG_range = CHAR(STRING_ELT(arg_range,0));
  std::string FLAG_annoType = CHAR(STRING_ELT(arg_annoType,0));

  std::set<std::string> FLAG_vcfColumn;
  std::vector<std::string> FLAG_infoTag, FLAG_indvTag;

  extractStringSet(arg_columns, &FLAG_vcfColumn);
  extractStringArray(arg_infoTag, &FLAG_infoTag);
  extractStringArray(arg_indvTag, &FLAG_indvTag);

  //fprintf(stdout, "range = %s\n", range.c_str());
  VCFExtractor vin(FLAG_fileName.c_str());
  if (FLAG_range.size())
    vin.setRangeList(FLAG_range.c_str());
  else {
    error("Please provide a range before we can continue\n");
  };

  if (FLAG_annoType.size()) {
    vin.setAnnoType(FLAG_annoType.c_str());
  }
  return readVCF2List(&vin, FLAG_vcfColumn, FLAG_infoTag, FLAG_indvTag);
}

/**
 * @param arg_fileName: a string character
 * @param arg_geneFile: which gene file to use
 * @param arg_geneName: which gene we are interested. (just allow One gene name).
 * @param arg_annoType: allow annotation type, can be regular expression. (e.g. Synonymous|Nonsynonymous)
 * @param arg_columns: a list of which columns to extract (e.g. CHROM, POS ...)
 * @param arg_infoTag: a list of which tag under INFO tag will be extracted (e.g. ANNO, ANNO_FULL, AC ...)
 * @param arg_indvTag: a list of which tag given in individual's column (e.g. GT, GD, GQ ...)
 */
SEXP impl_readVCFToListByGene(SEXP arg_fileName, SEXP arg_geneFile, SEXP arg_geneName, SEXP arg_annoType, SEXP arg_columns, SEXP arg_infoTag, SEXP arg_indvTag){

  // begin
  std::string FLAG_fileName = CHAR(STRING_ELT(arg_fileName,0));
  std::string FLAG_geneFile = CHAR(STRING_ELT(arg_geneFile,0));
  std::string FLAG_geneName = CHAR(STRING_ELT(arg_geneName,0));
  std::string FLAG_annoType = CHAR(STRING_ELT(arg_annoType,0));

  std::set<std::string> FLAG_vcfColumn;
  std::vector<std::string> FLAG_infoTag, FLAG_indvTag;

  extractStringSet(arg_columns, &FLAG_vcfColumn);
  extractStringArray(arg_infoTag, &FLAG_infoTag);
  extractStringArray(arg_indvTag, &FLAG_indvTag);

  // Rprintf("vcfColumn.size() = %u\n", FLAG_vcfColumn.size());
  // Rprintf("vcfInfo.size() = %u\n", FLAG_infoTag.size());
  // Rprintf("vcfIndv.size() = %u\n", FLAG_indvTag.size());
  // also append sample names at the end
  int retListLen = FLAG_vcfColumn.size() + FLAG_infoTag.size() + FLAG_indvTag.size() + 1;
  if (retListLen == 0) {
    return R_NilValue;
  }

  std::map< std::string, std::string> geneRange;
  loadGeneFile(FLAG_geneFile, FLAG_geneName, &geneRange);
  std::string range;
  for (std::map< std::string, std::string>::iterator it = geneRange.begin();
       it != geneRange.end();
       it++) {
    if (range.size() > 0) {
      range += ",";
    }
    range += it->second;
  };

  //fprintf(stdout, "range = %s\n", range.c_str());
  VCFExtractor vin(FLAG_fileName.c_str());
  if (range.size())
    vin.setRangeList(range.c_str());
  else {
    error("Please provide a range before we can continue\n");
  };

  if (FLAG_annoType.size()) {
    vin.setAnnoType(FLAG_annoType.c_str());
  }

  return readVCF2List(&vin, FLAG_vcfColumn, FLAG_infoTag, FLAG_indvTag);
} // end readVCF2List

/**
 * @param out will be a concatenated @param in separated by @param sep
 */
void set2string(const std::set<std::string>& in,
                std::string* out,
                const char sep) {
  out->clear();
  std::set<std::string>::const_iterator iter;
  for (iter = in.begin();
       iter != in.end();
       ++iter) {
    if (!out->empty()) {
      out->push_back(sep);
    };
    out->append(*iter);
  };
}

#define PVAL_FILE_CHROM_COL 0
#define PVAL_FILE_POS_COL 1
#define PVAL_FILE_REF_COL 2
#define PVAL_FILE_ALT_COL 3
#define PVAL_FILE_NSAMPLE_COL 4
#define PVAL_FILE_MAF_COL 5
#define PVAL_FILE_STAT_COL 6
#define PVAL_FILE_DIRECTION_COL 7
#define PVAL_FILE_PVAL_COL 8
#define PVAL_FILE_ANNO_COL 9
#define PVAL_FILE_ANNOFULL_COL 10
#define PVAL_FILE_FIELD_LEN 9  // required at least 9 columns (or it's invalid file)

#define COV_FILE_CHROM_COL 0
#define COV_FILE_BEGIN_POS_COL 1
#define COV_FILE_END_POS_COL 2
#define COV_FILE_GENE_COL 3
#define COV_FILE_POS_COL 4
#define COV_FILE_COV_COL 5

#define RET_REF_INDEX 0
#define RET_ALT_INDEX 1
#define RET_NSAMPLE_INDEX 2
#define RET_MAF_INDEX 3
#define RET_STAT_INDEX 4
#define RET_DIRECTION_INDEX 5
#define RET_PVAL_INDEX 6
#define RET_COV_INDEX 7
#define RET_POS_INDEX 8
#define RET_ANNO_INDEX 9

/**
 * Read serveral @param arg_pvalFile and equal number of @param arg_covFile,
 * return reesult for @arg_gene.
 * @return a list, each index is a gene.
 * then under gene, we have a list of values: numSample, maf, ...
 * then under gene/values level, we have a list (with length equals to numStudy) of double/int
 * NOTE:
 * we don't handle variants at duplicated positions.
 */
SEXP impl_rvMetaReadData(SEXP arg_pvalFile, SEXP arg_covFile, SEXP arg_gene) {
  int numAllocated = 0;
  SEXP ret = R_NilValue;

  // load by position
  std::set<std::string> FLAG_gene;
  std::vector<std::string> FLAG_pvalFile, FLAG_covFile;
  extractStringSet(arg_gene, &FLAG_gene);
  extractStringArray(arg_pvalFile, &FLAG_pvalFile);
  extractStringArray(arg_covFile, &FLAG_covFile);

  if (FLAG_pvalFile.size() != FLAG_covFile.size()){
    Rprintf("Unequal size!\n");
    error("Quitting");
  }

  int nStudy = FLAG_pvalFile.size();

  // get locations for the same gene from all provided cov files
  typedef std::map<std::string, std::map<std::string, int> > GeneLocationMap;
  GeneLocationMap geneLocationMap; // gene -> location -> location index
  std::set<std::string> locations; // only these locations will be considered in analysis.
  std::map<std::string, std::string> locationGeneMap; // location -> gene name
  std::map<std::string, int> geneIndex; // gene -> its index in returned list

  for (int i = 0; i < nStudy; ++i){
    LineReader lr(FLAG_covFile[i].c_str());
    std::vector<std::string> fd;
    std::vector<std::string> pos;
    std::string buf;
    while(lr.readLineBySep(&fd, " \t")){
      std::string& gene = fd[COV_FILE_GENE_COL];
      // if (gene != "CFH") continue;
      if (FLAG_gene.count(gene) == 0) continue;

      stringNaturalTokenize(fd[COV_FILE_POS_COL], ',', &pos);
      for (size_t i = 0; i < pos.size(); i++){
        // sprintf(buf, "%s:%s", fd[0].c_str(), pos[i].c_str());
        buf = fd[COV_FILE_CHROM_COL];
        buf += ':';
        buf += pos[i];
        if (geneLocationMap[gene].count(buf) == 0) {
          int s = geneLocationMap[gene].size();
          geneLocationMap[gene] [buf] = s;
          locations.insert(buf);
          // Rprintf("Gene %s location %s has index %d\n", gene.c_str(), buf.c_str(), geneLocationMap[gene][buf]);
        }
        if (locationGeneMap.count(buf) == 0) {
          locationGeneMap[buf] = gene;
        } else {
          // NOTE: need to take care of it later!
          // Rprintf("Same location [ %s ] have more than one gene [ %s ]\n", buf.c_str(), gene.c_str());
        };
        // Rprintf("add gene [ %s ] position %s\n", gene.c_str(), buf.c_str());
      }
    }
  }
  int nGene = geneLocationMap.size();
  Rprintf("%d gene loaded from covariate file.\n", nGene);

  PROTECT(ret = allocVector(VECSXP, nGene));
  numAllocated ++;
  // set list names
  SEXP geneNames;
  PROTECT(geneNames = allocVector(STRSXP, nGene));
  numAllocated ++;
  {
    int i = 0;
    for ( GeneLocationMap::iterator iter = geneLocationMap.begin();
          iter != geneLocationMap.end() ; ++iter){
      // Rprintf("assign gene name: %s\n", iter->first.c_str());
      SET_STRING_ELT(geneNames, i, mkChar(iter->first.c_str()));
      geneIndex[iter->first.c_str()] = i;
      ++i;
    }
  }
  setAttrib(ret, R_NamesSymbol, geneNames);

  // create n, maf, p, cov
  std::vector<std::string> names;
  names.push_back("ref");
  names.push_back("alt");
  names.push_back("nSample");
  names.push_back("maf");
  names.push_back("stat");
  names.push_back("direction");
  names.push_back("pVal");
  names.push_back("cov");
  names.push_back("pos");
  names.push_back("anno");

  GeneLocationMap::iterator iter = geneLocationMap.begin();
  for ( int i = 0;
        iter != geneLocationMap.end() ; ++iter, ++i){
    SEXP s = VECTOR_ELT(ret, i);
    numAllocated += createList(names.size(), &s); // a list with 10 elements: ref, alt, n, maf, stat, direction, p, cov, pos, anno
    numAllocated += setListNames(names, &s);

    SEXP ref, alt, n, maf, stat, direction, p, cov, pos, anno;
    numAllocated += createList(nStudy, &ref);
    numAllocated += createList(nStudy, &alt);
    numAllocated += createList(nStudy, &n);
    numAllocated += createList(nStudy, &maf);
    numAllocated += createList(nStudy, &stat);
    numAllocated += createList(nStudy, &direction);
    numAllocated += createList(nStudy, &p);
    numAllocated += createList(nStudy, &cov);
    numAllocated += createList(nStudy, &pos);
    numAllocated += createList(nStudy, &anno);

    int npos = iter->second.size();
    for (int j = 0; j < nStudy; ++j) {
      SEXP t;
      numAllocated += createStringArray(npos, &t);
      initStringArray(t);
      SET_VECTOR_ELT(ref, j, t);

      numAllocated += createStringArray(npos, &t);
      initStringArray(t);
      SET_VECTOR_ELT(alt, j, t);

      numAllocated += createIntArray(npos, &t);
      initIntArray(t);
      SET_VECTOR_ELT(n, j, t);

      numAllocated += createDoubleArray(npos, &t);
      initDoubleArray(t);
      SET_VECTOR_ELT(maf, j, t);

      numAllocated += createDoubleArray(npos, &t);
      initDoubleArray(t);
      SET_VECTOR_ELT(stat, j, t);

      numAllocated += createIntArray(npos, &t);
      initIntArray(t);
      SET_VECTOR_ELT(direction, j, t);

      numAllocated += createDoubleArray(npos, &t);
      initDoubleArray(t);
      SET_VECTOR_ELT(p, j, t);

      // Rprintf("Create double array %d for study %d\n", npos * npos, j);
      numAllocated += createDoubleArray(npos*npos, &t);
      numAllocated += setDim(npos, npos, &t);
      initDoubleArray(t);
      SET_VECTOR_ELT(cov, j, t);
    }
    numAllocated += createStringArray(npos, &pos);
    initStringArray(pos);

    numAllocated += createStringArray(npos, &anno);
    initStringArray(anno);

    SET_VECTOR_ELT(s, RET_REF_INDEX, ref);
    SET_VECTOR_ELT(s, RET_ALT_INDEX, alt);
    SET_VECTOR_ELT(s, RET_NSAMPLE_INDEX, n);
    SET_VECTOR_ELT(s, RET_MAF_INDEX, maf);
    SET_VECTOR_ELT(s, RET_STAT_INDEX, stat);
    SET_VECTOR_ELT(s, RET_DIRECTION_INDEX, direction);
    SET_VECTOR_ELT(s, RET_PVAL_INDEX, p);
    SET_VECTOR_ELT(s, RET_COV_INDEX, cov);
    SET_VECTOR_ELT(s, RET_POS_INDEX, pos);
    SET_VECTOR_ELT(s, RET_ANNO_INDEX, anno);

    SET_VECTOR_ELT(ret, i, s);
  };

  std::map< std::string, std::set<std::string> > posAnnotationMap;

  Rprintf("Read score tests...\n");
  // read pval file and fill in values
  std::string p;
  double tempDouble;
  int tempInt;
  for (int study = 0; study < nStudy; ++study) {
    Rprintf("In study %d\n", study);
    LineReader lr(FLAG_pvalFile[study].c_str());
    std::vector<std::string> fd;
    std::set< std::string> processedSite;
    while (lr.readLineBySep(&fd, " \t")){
      if (fd.size() < PVAL_FILE_FIELD_LEN) {
        Rprintf("File format is incorrect [ %s ], skipped\n", FLAG_pvalFile[study].c_str());
        break;
      }
      p = fd[PVAL_FILE_CHROM_COL];
      p += ':';
      p += fd[PVAL_FILE_POS_COL];
      if (processedSite.count(p) == 0) {
        processedSite.insert(p);
      } else{
        Rprintf("Position %s appeared more than once, skipping...\n", p.c_str());
        continue;
      }

      SEXP u, v, s;
      std::string& gene = locationGeneMap[p];
      if (FLAG_gene.count(gene) == 0) continue;
      if (geneLocationMap.count(gene) == 0 ||
          geneLocationMap[gene].count(p) == 0) continue; // skip non existing position
      int idx = geneLocationMap[gene][p];

      // Rprintf("working on index %d\n", idx);
      u = VECTOR_ELT(ret, geneIndex[gene]);
      v = VECTOR_ELT(u, RET_REF_INDEX);
      s = VECTOR_ELT(v, study); // ref
      SET_STRING_ELT(s, idx, mkChar(fd[PVAL_FILE_REF_COL].c_str()));

      u = VECTOR_ELT(ret, geneIndex[gene]);
      v = VECTOR_ELT(u, RET_ALT_INDEX);
      s = VECTOR_ELT(v, study); // alt
      SET_STRING_ELT(s, idx, mkChar(fd[PVAL_FILE_ALT_COL].c_str()));
      
      if ( str2int(fd[PVAL_FILE_NSAMPLE_COL], &tempInt) ) {
        u = VECTOR_ELT(ret, geneIndex[gene]);
        v = VECTOR_ELT(u, RET_NSAMPLE_INDEX);
        s = VECTOR_ELT(v, study); // n
        INTEGER(s)[idx] = tempInt;
      }

      if ( str2double(fd[PVAL_FILE_MAF_COL], &tempDouble) ) {
        v = VECTOR_ELT(u, RET_MAF_INDEX);
        s = VECTOR_ELT(v, study); // maf
        REAL(s)[idx] = tempDouble;
      }

      if ( str2double( fd[PVAL_FILE_STAT_COL], & tempDouble)) {
        v = VECTOR_ELT(u, RET_STAT_INDEX);
        s = VECTOR_ELT(v, study); // stat
        REAL(s)[idx] = tempDouble;
      }

      v = VECTOR_ELT(u, RET_DIRECTION_INDEX);
      s = VECTOR_ELT(v, study); // direction
      if (fd[PVAL_FILE_DIRECTION_COL] == "+") {
        INTEGER(s)[idx] = 1;
      } else if (fd[PVAL_FILE_DIRECTION_COL] == "-") {
        INTEGER(s)[idx] = -1;
      }

      if ( str2double(fd[PVAL_FILE_PVAL_COL], & tempDouble)) {
        // Rprintf("Set pval index");
        
        v = VECTOR_ELT(u, RET_PVAL_INDEX);
        s = VECTOR_ELT(v, study); // pval
        REAL(s)[idx] = tempDouble;
      };

      // // Rprintf("Set pos index");
      // v = VECTOR_ELT(u, RET_POS_INDEX);
      // s = VECTOR_ELT(v, study); // position
      // SET_STRING_ELT(s, idx, mkChar(p.c_str()));

      // // Rprintf("Set anno index");
      // if (fd.size() >= PVAL_FILE_ANNO_COL + 1) {
      //   v = VECTOR_ELT(u, RET_ANNO_INDEX);
      //   s = VECTOR_ELT(v, study); // anno
      //   SET_STRING_ELT(s, idx, mkChar(fd[PVAL_FILE_ANNO_COL].c_str()));
      // }
      if (fd.size() >= PVAL_FILE_ANNO_COL + 1) {
        const std::string& s = fd[PVAL_FILE_ANNO_COL];
        if (!posAnnotationMap[p].count(s))  {
          posAnnotationMap[p].insert(s);
        }
      }
    }; // end while
    Rprintf("Done read file: %s\n", FLAG_pvalFile[study].c_str());
  }

  // fill in position and annotation
  iter = geneLocationMap.begin();
  for ( int i = 0;
        iter != geneLocationMap.end() ; ++iter, ++i){
    std::map<std::string, int>& loc2idx = iter->second;
    SEXP u = VECTOR_ELT(ret, geneIndex[iter->first]);
    SEXP pos = VECTOR_ELT(u, RET_POS_INDEX);
    SEXP anno = VECTOR_ELT(u, RET_ANNO_INDEX);

    for (std::map<std::string, int>::const_iterator it = loc2idx.begin();
         it != loc2idx.end();
         ++it) {
      int idx = it->second;
      if (locations.count(it->first)) {
        SET_STRING_ELT(pos, idx, mkChar(it->first.c_str()));
      }
      if (posAnnotationMap.count(it->first)) {
        std::string ret;
        set2string(posAnnotationMap[it->first], &ret, ',');
        SET_STRING_ELT(anno, idx, mkChar(ret.c_str()));
      };
    }
  }
  // read cov file and record pos2idx
  Rprintf("Read cov files ... \n");
  for (int study = 0; study < nStudy; ++study) {
    LineReader lr(FLAG_covFile[study].c_str());
    std::vector<std::string> fd;
    std::vector<std::string> pos;
    std::vector<double> cov;
    std::string p;

    while (lr.readLineBySep(&fd, " \t")){
      std::string& g = fd[COV_FILE_GENE_COL];
      if (geneLocationMap.count(g) == 0) continue;
      // Rprintf("begin parse pos, cov\n");
      std::string& chr = fd[COV_FILE_CHROM_COL];
      // Rprintf("pos: %s\n", fd[COV_FILE_COV_COL].c_str());
      stringNaturalTokenize(fd[COV_FILE_COV_COL], ',', &pos); /// temporary use variable pos to store cov
      cov.resize(pos.size());
      for (size_t i = 0; i < pos.size(); ++i){
        cov[i] = atof(pos[i]);
      }
      stringNaturalTokenize(fd[COV_FILE_POS_COL], ',', &pos);
      for (size_t i = 0; i < pos.size(); ++i) {
        pos[i] = chr + ":" + pos[i];
      };

      int k = 0;
      int covLen = geneLocationMap[g].size();
      SEXP u, v, s;
      u = VECTOR_ELT(ret, geneIndex[g]);
      v = VECTOR_ELT(u, RET_COV_INDEX);  // cov is the 6th element in the list
      // Rprintf("Gene [%s] has %d locations \n", g.c_str(), covLen);
      for (size_t i = 0; i < pos.size(); i++) {
        for (size_t j = i; j < pos.size(); j++) {
          std::string& pi = pos[i];
          std::string& pj = pos[j];
          int posi = geneLocationMap[g][pi];
          int posj = geneLocationMap[g][pj];

          s = VECTOR_ELT(v, study);

          REAL(s) [posi * covLen + posj] = cov[k];
          REAL(s) [posj * covLen + posi] = cov[k];
          ++k;
        }
      }
    }
  }
  Rprintf("Finished calculation.\n");
  UNPROTECT(numAllocated);
  return ret;
} // SEXP rvMetaRead2List(SEXP arg_pvalFile, SEXP arg_covFile) {

/**
 * @return a covariance matrix
 */
SEXP impl_readCovByRange(SEXP arg_covFile, SEXP arg_range) {
  int numAllocated = 0;
  SEXP ret = R_NilValue;

  // load by position
  std::string FLAG_covFile;
  extractString(arg_covFile, &FLAG_covFile);
  std::string FLAG_range;
  extractString(arg_range, &FLAG_range);

  // Rprintf("open %s\n", FLAG_covFile.c_str());
  tabix_t* t = ti_open(FLAG_covFile.c_str(), 0);
  if (t == 0) {
    REprintf("Cannot open %s file!\n", FLAG_covFile.c_str());
    return ret;
  }
  // Rprintf("open OK\n");

  if (ti_lazy_index_load(t) < 0) {
    REprintf("[tabix] failed to load the index file.\n");
    return ret;
  }

  std::string line;
  std::vector<std::string> fd;
  std::string chrom;
  std::map< std::string, int> pos2idx;
  std::vector<std::string> pos;
  std::vector<double> cov;

  int lineNo = 0;
  ti_iter_t iter;
  int tid, beg, end;
  const char* s;
  int len;
  // Rprintf("begin parse %s ..\n", FLAG_range.c_str());
  if (ti_parse_region(t->idx, FLAG_range.c_str(), &tid, &beg, &end) == 0) {
    // Rprintf("parse OK\n");
    iter = ti_queryi(t, tid, beg, end);
    while ((s = ti_read(t, iter, &len)) != 0) {
      // fputs(s, stdout); fputc('\n', stdout);
      // Rprintf("%s\n", s);
      line = s;
      stringTokenize(line, "\t", &fd);
      if (fd[3] == "1") continue; // only self covariance
      lineNo ++;

      if (chrom.empty()){
        chrom = fd[0];
      } else{
        if (chrom != fd[0]){
          REprintf("chromosome does not match %s and %s.\n", chrom.c_str(), fd[0].c_str());
          return ret;
        }
      }
      // begPos = fd[1];
      // if (pos2idx.count(begPos) == 0) {
      //   pos2idx[begPos] = pos2idx.size() - 1; // we want 0-based index, and pos2idx[begPos] is already created when call size()
      // }
      stringTokenize(fd[4], ',', &pos);
      for (size_t i = 0; i < pos.size(); ++i){
        // Rprintf("pos[%zu] = %s\n", i, pos[i].c_str());
        if (pos2idx.count(pos[i]) == 0)
          pos2idx[pos[i]] = pos2idx.size() - 1;
      }
    }
    ti_iter_destroy(iter);
    // Rprintf("parse end\n");
  }  else {
    REprintf("invalid region: unknown target name or minus interval.\n");
    return ret;
  }
  int retDim = pos2idx.size();
  Rprintf("Total %d line loaded, now put them to matrix [ %d x %d ] in R ...\n", lineNo, retDim, retDim);
  // // Rprintf("pos2idx.size() = %zu \n", pos2idx.size());
  // for (std::map<std::string, int>::const_iterator iter = pos2idx.begin();
  //      iter != pos2idx.end();
  //      ++iter) {
  //   // Rprintf("%s -> %d\n", iter->first.c_str(), iter->second);
  // }

  // init R matrix
  numAllocated += createDoubleArray(retDim * retDim, &ret);
  initDoubleArray(ret);

  // store results
  int row, col;
  if (ti_parse_region(t->idx, FLAG_range.c_str(), &tid, &beg, &end) == 0) {
    // Rprintf("parse OK\n");
    iter = ti_queryi(t, tid, beg, end);
    while ((s = ti_read(t, iter, &len)) != 0) {
      line = s;
      stringTokenize(line, "\t", &fd);
      if (fd[3] == "1" && pos2idx.count(fd[1]) == 0) continue; // only itslef covariance
      // begPos = fd[1];
      // if (pos2idx.find(begPos) == pos2idx.end()){
      //   REprintf("Cannot find row position: %s\n", begPos.c_str());
      // }
      // row = pos2idx[begPos];

      stringTokenize(fd[5], ",", &pos);
      cov.resize(pos.size());
      for (size_t i = 0; i < pos.size(); ++i){
        cov[i] = atof(pos[i]);
      }
      stringTokenize(fd[4], ",", &pos);
      if (pos.size() != cov.size()) {
        REprintf("mismatch lenght\n");
      };

      row = pos2idx[pos[0]];
      for (size_t i = 0; i < pos.size(); ++i){
        if (pos2idx.find(pos[i]) == pos2idx.end()){
          REprintf("Cannot find col position: %s\n", pos[i].c_str());
        }
        col = pos2idx[pos[i]];

        // Rprintf("%d: %d %d = %g\n", i, row, col, cov[i]);
        REAL(ret)[row * retDim + col] = cov[i];
        REAL(ret)[col * retDim + row] = cov[i];
      };
    }
    ti_iter_destroy(iter);
  }
  ti_close(t);

  // set dim info
  numAllocated += setDim(retDim, retDim, &ret);
  // set matrix label
  SEXP rowName;
  PROTECT(rowName=allocVector(STRSXP, retDim));
  numAllocated += 1;
  std::string label;
  for (std::map<std::string, int>::const_iterator iter = pos2idx.begin();
       iter != pos2idx.end();
       ++iter) {
    label = chrom;
    label += ':';
    label += iter->first;
    SET_STRING_ELT(rowName, iter->second, mkChar(label.c_str()));
  }

  SEXP dimnames;
  PROTECT(dimnames = allocVector(VECSXP, 2));
  numAllocated += 1;
  SET_VECTOR_ELT(dimnames, 0, rowName);
  SET_VECTOR_ELT(dimnames, 1, rowName);
  setAttrib(ret, R_DimNamesSymbol, dimnames);

  UNPROTECT(numAllocated);
  return ret;
}

/**
 * @return 0 if success; -1 if not oK
 */
int parsePosition(const std::string& range,
                  std::string* chrom,
                  int* beg,
                  int* end){

  std::string r;
  r = chopChr(range);
  size_t i = r.find(':');
  if (i == std::string::npos)
    return -1;

  *chrom = r.substr(0, i);

  size_t j = r.find('-', i+1);
  if (j == std::string::npos) { // 1:100
    *beg = atoi(r.substr(i, r.size() - i));
    *end = INT_MAX;
    return 0;
  }
  // 1:100-200
  //  ^   ^   ^
  //  i   j
  //  1   5   8
  *beg = atoi(r.substr(i+1, j - i - 1));
  *end = atoi(r.substr(j+1, r.size() - j));
  return 0;
}
/**
 * @param arg_scoreFile: single variant score test file
 * @return a data frame of
 */
SEXP impl_readScoreByRange(SEXP arg_scoreFile, SEXP arg_range) {
  int numAllocated = 0;
  SEXP ret = R_NilValue;

  // load by position
  std::string FLAG_scoreFile;
  extractString(arg_scoreFile, &FLAG_scoreFile);
  std::string FLAG_range;
  extractString(arg_range, &FLAG_range);

  // parse region
  std::string chrom;
  int beg;
  int end;
  parsePosition(FLAG_range, &chrom, &beg, &end);
  // Rprintf("chrom = %s, beg = %d, end = %d", chrom.c_str(), beg, end);
  
  // set up return values

  std::vector <std::string> fd;
  LineReader lr(FLAG_scoreFile);
  int lineNo = 0;
  int fieldLen = -1;
  std::vector <int> position;
  std::vector <std::string> ref;
  std::vector <std::string> alt;
  std::vector <int> nsample;
  std::vector <double> maf;
  std::vector <double> stat;
  std::vector <std::string> direction;
  std::vector <double> pval;
  std::vector <std::string> anno;
  std::vector <std::string> annoFull;
  while (lr.readLineBySep(&fd, "\t ")) {
    lineNo ++;
    if (fd[0] != chrom) continue;
    int pos = atoi(fd[1]);
    if ( pos < beg || pos > end) continue;
    if (fieldLen < 0) {
      fieldLen = fd.size();
    } else{
      if (fieldLen <= PVAL_FILE_PVAL_COL ) {
        REprintf("skip insufficient fields at line %d.\n", lineNo);
        continue;
      }
      if (fieldLen != (int)fd.size()) {
        REprintf("Inconsistent field length at line %d.\n", lineNo);
        return ret;
      }
    }
    position.push_back(pos);
    ref.push_back(fd[PVAL_FILE_REF_COL]);
    alt.push_back(fd[PVAL_FILE_ALT_COL]);
    nsample.push_back(atoi(fd[PVAL_FILE_NSAMPLE_COL]));
    maf.push_back(atof(fd[PVAL_FILE_MAF_COL]));
    stat.push_back(atof(fd[PVAL_FILE_STAT_COL]));
    direction.push_back(fd[PVAL_FILE_DIRECTION_COL]);
    pval.push_back(atof(fd[PVAL_FILE_PVAL_COL]));

    if (fd.size() > PVAL_FILE_ANNOFULL_COL) {
      anno.push_back(fd[PVAL_FILE_ANNO_COL]);
      annoFull.push_back(fd[PVAL_FILE_ANNOFULL_COL]);
    }
  };

  if (fieldLen < 0) {
    if (lineNo == 0) {
      REprintf("No valid input line read, please check your file [ %s ]\n", FLAG_scoreFile.c_str());
    }
    return ret;
  };

  int retListLen;
  if (anno.size() ) {
    retListLen = 10; // hard coded number
  } else {
    retListLen = 8;
  }
  PROTECT(ret = allocVector(VECSXP, retListLen));
  numAllocated ++;

  std::vector<std::string> listNames;
  int retListIdx = 0;
  numAllocated += storeResult(position, ret, retListIdx++);  
  numAllocated += storeResult(ref, ret, retListIdx++);
  numAllocated += storeResult(alt, ret, retListIdx++);
  numAllocated += storeResult(nsample, ret, retListIdx++);
  numAllocated += storeResult(maf, ret, retListIdx++);
  numAllocated += storeResult(stat, ret, retListIdx++);
  numAllocated += storeResult(direction, ret, retListIdx++);
  numAllocated += storeResult(pval, ret, retListIdx++);
  listNames.push_back("pos");
  listNames.push_back("ref");  
  listNames.push_back("alt");
  listNames.push_back("nSample");  
  listNames.push_back("maf");
  listNames.push_back("stat");
  listNames.push_back("direction");
  listNames.push_back("pVal");
  
  if (anno.size() ) {
    numAllocated += storeResult(anno, ret, retListIdx++);
    numAllocated += storeResult(annoFull, ret, retListIdx++);
    listNames.push_back("anno");
    listNames.push_back("annoFull");
  }
  SEXP sListNames;
  PROTECT(sListNames = allocVector(STRSXP, listNames.size()));
  numAllocated ++;
  for (unsigned int i = 0; i != listNames.size(); ++i){
    SET_STRING_ELT(sListNames, i, mkChar(listNames[i].c_str()));
  }
  setAttrib(ret, R_NamesSymbol, sListNames);

  UNPROTECT(numAllocated);
  return ret;
}
