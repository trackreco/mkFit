#ifndef _buildtest_
#define _buildtest_

#include <map>
#include <string>

class TTree;
class TH1F;

unsigned int getPhiPartition(float phi);
void runBuildingTest(bool saveTree, unsigned int nevts);
void runBuildingTest(bool saveTree, TTree *tree,unsigned int& tk_nhits, float& tk_chi2, std::map<std::string,TH1F*>& validation_hists);

#endif
