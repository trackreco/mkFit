void plotFit(){

  gROOT->ProcessLine(".L plotTree.cpp++");
  gSystem->Load("plotTree_cpp.so"); 
  plotTree();

}
