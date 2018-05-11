#include "StackValidation.hh"
#include "StackValidation.cpp"

void makeValidation(const TString & label = "", const TString & extra = "", const Bool_t cmsswComp = false)
{
  StackValidation Stacks(label,extra,cmsswComp);
  Stacks.MakeValidationStacks();
}
