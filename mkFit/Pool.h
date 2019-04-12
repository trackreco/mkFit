#ifndef Pool_h
#define Pool_h
#include "tbb/concurrent_queue.h"

namespace mkfit {

template <typename TT>
struct Pool
{
  typedef std::function<TT*()>     CFoo_t;
  typedef std::function<void(TT*)> DFoo_t;

  CFoo_t m_create_foo  = []()     { return new (_mm_malloc(sizeof(TT), 64)) TT; };
  DFoo_t m_destroy_foo = [](TT* x){ x->~TT(); _mm_free(x); };

  tbb::concurrent_queue<TT*> m_stack;

  size_t size() { return m_stack.unsafe_size(); }

  void populate(int threads = Config::numThreadsFinder)
  {
    for (int i = 0; i < threads; ++i)
    {
      m_stack.push(m_create_foo());
    }
  }

  Pool() {}
  Pool(CFoo_t cf, DFoo_t df) : m_create_foo(cf), m_destroy_foo(df) {}

  ~Pool()
  {
    TT *x;
    while (m_stack.try_pop(x))
    {
      m_destroy_foo(x);
    }
  }

  void SetCFoo(CFoo_t cf) { m_create_foo  = cf; }
  void SetDFoo(DFoo_t df) { m_destroy_foo = df; }

  TT* GetFromPool()
  {
    TT *x;
    if (m_stack.try_pop(x)) {
      return x;
    } else {
      return m_create_foo();
    }
  }

  void ReturnToPool(TT *x)
  {
    m_stack.push(x);
  }
};


} // end namespace mkfit
#endif
