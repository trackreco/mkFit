#include "Rtypes.h"
#include "vector"

#ifdef __CINT__
#pragma link C++ class vector<vector<int> >+ ;
#ifdef G__VECTOR_HAS_CLASS_ITERATOR
#pragma link C++ operators vector<vector<int> >::iterator;
#pragma link C++ operators vector<vector<int> >::const_iterator;
#pragma link C++ operators vector<vector<int> >::reverse_iterator;
#endif
#endif
