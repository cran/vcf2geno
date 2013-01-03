#ifndef _R_CPP_INTERFACE_H_
#define _R_CPP_INTERFACE_H_

#include <R.h>
#include <Rdefines.h> // define SEXP

#include <iostream>
#include <string>
#include <vector>
#include <set>



void extractString(SEXP in, std::string* out);

/**
 * NOTE
 */
void extractStringArray(SEXP in, std::vector<std::string>* out);

void extractStringSet(SEXP in, std::set<std::string>* out);

/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char *str);

void dump(std::vector<std::string> & s);

int storeResult(const std::vector<std::string>& in , SEXP ret, int idx) ;

int storeResult(const std::vector<int>& in , SEXP& ret, int idx);

int storeIntResult(const std::vector<std::string>& in , SEXP& ret, int idx);

int storeResult(const std::vector<double>& in , SEXP& ret, int idx);

int storeDoubleResult(const std::vector<std::string>& in , SEXP& ret, int idx); 
int storeResult(const std::string& key, const std::vector<std::string>& val , SEXP ret, int idx); 

int storeResult(const std::string& key, const std::vector<int>& val , SEXP& ret, int idx);

int setDim(int nrow, int ncol, SEXP* s);

/**
 * Set dim attributes for ret[idx]
 */
int setDim(int nrow, int ncol, SEXP ret, int idx);

int createList(int n, SEXP* s); 

int createStringArray(int n, SEXP* s);

int createDoubleArray(int n, SEXP* s);

int createIntArray(int n, SEXP* s);

int setListNames(std::vector<std::string>& names, SEXP* s);

void initDoubleArray(SEXP s); 

void initIntArray(SEXP s); 

void initStringArray(SEXP s); 


#endif /* _R_CPP_INTERFACE_H_ */
