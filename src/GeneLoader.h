#include <string>
#include <set>
#include <vector>
#include <map>
void loadGeneFile(const std::string& geneFile, const std::set<std::string>& geneName, std::map< std::string, std::string>* geneRangePtr);
void loadGeneFile(const std::string& geneFile, const std::string& geneName, std::map< std::string, std::string>* geneRangePtr);
void loadGeneFile(const std::string& geneFile, const std::vector<std::string>& geneName, std::map< std::string, std::string>* geneRangePtr);
