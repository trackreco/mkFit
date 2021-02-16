#ifndef ConfigWrapper_h
#define ConfigWrapper_h

namespace mkfit {
  /**
   * The purpose of this namespace is to hide the header of Config.h
   * from CMSSW. This header contain uses of the build-time
   * configuration macros, that should remain as internal details of
   * MkFit package.
   */
  namespace ConfigWrapper {
    void initializeForCMSSW(bool silent);

    void setNTotalLayers(int nTotalLayers);
  }
}

#endif
