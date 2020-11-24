#include "MkBuilderWrapper.h"
#include "MkBuilder.h"

namespace mkfit {
  MkBuilderWrapper::MkBuilderWrapper():
    builder_(MkBuilder::make_builder())
  {}

  MkBuilderWrapper::~MkBuilderWrapper() {}

  void MkBuilderWrapper::populate() {
    MkBuilder::populate();
  }
}
