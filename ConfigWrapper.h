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
    enum class SeedCleaningOpts {
      noCleaning,
      cleanSeedsN2
    };
    enum class BackwardFit {
      noFit,
      toFirstLayer,
      toPCA
    };

    void initializeForCMSSW(SeedCleaningOpts seedClean, BackwardFit backfit, bool silent);
    void setRemoveDuplicates(bool removeDuplicates);

    void setNTotalLayers(int nTotalLayers);
  }
}

#endif
