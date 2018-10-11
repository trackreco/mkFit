#ifndef MatriplexPackers_h
#define MatriplexPackers_h

#include "Matrix.h"

#include "Hit.h"
#include "Track.h"

namespace mkfit {


//==============================================================================
// MatriplexErrParPackerSlurpIn
//==============================================================================

template<typename D>
class MatriplexPackerSlurpIn
{
protected:
   alignas(64) int m_idx[NN];

   const D *m_base;
   int      m_pos;

public:
   MatriplexPackerSlurpIn() :
      m_base (0),
      m_pos  (0)
   {}

   void Reset()        { m_pos = 0; }

   void ResetBase()    { m_pos = 0; m_base = 0; }

   void AddNullInput() { m_idx[m_pos++] = 0; }

   void AddInput(const D& item)
   {
      if (m_base == 0)
      {
         m_base = & item;
      }

      // Could issue prefetch requests here.

      m_idx[m_pos] = & item - m_base;

      ++m_pos;
   }

   void AddInputAt(int pos, const D& item)
   {
      while (m_pos < pos)
      {
         // We might not care about initialization / reset to 0.
         // Or we could be building an additional mask (on top of N_proc).
         m_idx[m_pos++] = 0;
      }

      AddInput(item);
   }

   template<typename TM>
   void Pack(TM &mplex, int base_offset)
   {
      assert (m_pos <= NN);

      if (m_pos == 0)
      {
         // throw? ... this is probably not expected to happen.
         return;
      }

#if defined(GATHER_INTRINSICS)
      GATHER_IDX_LOAD(vi, m_idx);
      mplex.SlurpIn(m_base + base_offset, vi, sizeof(D), m_pos);
#else
      mplex.SlurpIn(m_base + base_offset, m_idx, m_pos);
#endif
   }
};


//==============================================================================
// MatriplexErrParPackerSlurpIn
//==============================================================================

// T - input class (Track or Hit), D - data type (float)

// Would this actually be prettier without the inheriatance (and all this-> mess)?

template<typename T, typename D>
class MatriplexErrParPackerSlurpIn : public MatriplexPackerSlurpIn<D>
{
   // alignas(64) int m_idx[NN];

   // const D *m_base;
   // int      m_pos;
   // int      m_off_error;
   int      m_off_param;

public:
   MatriplexErrParPackerSlurpIn() : MatriplexPackerSlurpIn<D>()
      // m_base (0),
      // m_pos  (0)
   {}

   // void Reset()        { m_pos = 0; }

   // void ResetBase()    { m_pos = 0; m_base = 0; }

   // void AddNullInput() { m_idx[m_pos++] = 0; }

   void AddInput(const T& item)
   {
      if (this->m_base == 0)
      {
         // One could use item.errArray() as base and then we'd only need
         // m_off_param. I'm avoiding this as I'm not sure if this
         // requires fetching of item (also when base != 0).
         // Nope ... this compiles into a constant.

         // this->m_base = (D*) &item;
         // m_off_error  = item.errArray() - this->m_base;

         this->m_base = item.errArray();

         m_off_param  = item.posArray() - this->m_base;
      }

      // Could issue L1 prefetch requests here.

      this->m_idx[this->m_pos] = item.errArray() - this->m_base;

      ++this->m_pos;
   }

   void AddInputAt(int pos, const T& item)
   {
      while (this->m_pos < pos)
      {
         // We might not care about initialization / reset to 0.
         // Or we could be building an additional mask (on top of N_proc).
         this->m_idx[this->m_pos++] = 0;
      }

      AddInput(item);
   }

   template<typename TMerr, typename TMpar>
   void Pack(TMerr &err, TMpar &par)
   {
      assert (this->m_pos <= NN);

      if (this->m_pos == 0)
      {
         // throw? ... this is probably not expected to happen.
         return;
      }

#if defined(GATHER_INTRINSICS)
      GATHER_IDX_LOAD(vi, this->m_idx);
      err.SlurpIn(this->m_base,               vi, sizeof(D), this->m_pos);
      par.SlurpIn(this->m_base + m_off_param, vi, sizeof(D), this->m_pos);
#else
      err.SlurpIn(this->m_base,               this->m_idx, this->m_pos);
      par.SlurpIn(this->m_base + m_off_param, this->m_idx, this->m_pos);
#endif
   }
};


//==============================================================================
// MatriplexTrackPackerPlexify
//==============================================================================

class MatriplexTrackPackerPlexify // : public MatriplexTrackPackerBase
{
public:
   MatriplexTrackPackerPlexify() // : MatriplexPlexifyPackerBase(N_proc)
   {}

   void AddInput(const Track &track)
   {
   }

   void Pack(MPlexLS &err, MPlexLV &par)
   {
   }
};


//==============================================================================
// Packer Selection
//==============================================================================

// Optionally ifdef with defines from Makefile.config

using MatriplexHitPacker   = MatriplexErrParPackerSlurpIn<Hit,   float>;
using MatriplexTrackPacker = MatriplexErrParPackerSlurpIn<Track, float>;

using MatriplexHoTPacker   = MatriplexPackerSlurpIn<HitOnTrack>;
}

#endif
