#ifndef _buildtest_
#define _buildtest_

#include "Track.h"
#include "Matrix.h"
#include "TFile.h"
#include "TTree.h"

bool sortByPhi(Hit hit1,Hit hit2);
unsigned int getPhiPartition(float phi);
void runBuildingTest(bool saveTree, unsigned int nevts);
void runBuildingTest(bool saveTree, TTree *tree,unsigned int& tk_nhits, float& tk_chi2);

#endif
