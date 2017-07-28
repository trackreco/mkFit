#include "TTreeValidation.h"

Validation* Validation::make_validation(const std::string& fileName)
{
#ifndef NO_ROOT
  if (Config::root_val || Config::fit_val || Config::cmssw_val) {
    return new TTreeValidation(fileName);
  }
#endif
  return new Validation();
}

Validation::Validation() {}
