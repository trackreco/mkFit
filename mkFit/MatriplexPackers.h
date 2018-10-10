#ifndef MatriplexPackers_h
#define MatriplexPackers_h

#include "Matrix.h"

#include "Hit.h"
#include "Track.h"

namespace mkfit {


//==============================================================================
// MatriplexErrParPackerSlurpIn
//==============================================================================

// Extract base class to also support HitOnTrack slurping.

// T - input class, D - data type

template<typename T, typename D>
class MatriplexErrParPackerSlurpIn // : public MatriplexPackerSlurpInBase
{
   D    *m_base;     // Template to matriplex type, use as scale.
                     // Need typedef in Matrix.h (and pass it in Packer Selection
                     // section below).
   int   m_pos;
   int   m_off_error;
   int   m_off_param;

   int   m_idx[NN]  __attribute__((aligned(64)));

public:
   MatriplexErrParPackerSlurpIn() :
      m_base (0),
      m_pos  (0)
   {}

   void Reset()      { m_pos = 0; }

   void ResetBase()  { m_pos = 0; m_base = 0; }

   void AddInput(const T& item)
   {
      if (m_base == 0)
      {
         m_base      = (D*) &item;
         m_off_error = item.errArray() - m_base;
         m_off_param = item.posArray() - m_base;
      }

      // Could issue prefetch requests here.

      m_idx[m_pos] = (D*) &item - m_base;

      ++m_pos;
   }

   void AddInputAt(int pos, const T& item)
   {
      while (m_pos < pos)
      {
         // We might not care about initialization / reset to 0.
         // Or we could be building an additional mask (on top of N_proc).
         m_idx[m_pos++] = 0;
      }

      AddInput(item);
   }

   template<typename TMerr, typename TMpar>
   void Pack(TMerr &err, TMpar &par)
   {
      assert (m_pos <= NN);

      if (m_pos == 0)
      {
         // throw ... this is probably not expected to happen.

         return;
      }

#ifdef MIC_INTRINSICS
      __m512i vi = _mm512_load_epi32(m_idx);
      err.SlurpIn(m_base + m_off_error, vi, sizeof(D), m_pos);
      par.SlurpIn(m_base + m_off_param, vi, sizeof(D), m_pos);
#else
      err.SlurpIn(m_base + m_off_error, m_idx, m_pos);
      par.SlurpIn(m_base + m_off_param, m_idx, m_pos);
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

using MatriplexTrackPacker = MatriplexErrParPackerSlurpIn<Track, float>;

}

#endif
