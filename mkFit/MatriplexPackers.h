#ifndef MatriplexPackers_h
#define MatriplexPackers_h

#include "Matrix.h"

#include "Hit.h"
#include "Track.h"

namespace mkfit {


//==============================================================================
// MatriplexTrackPacker - do we even need a baseclass?
//==============================================================================

class MatriplexTrackPackerBase
{
protected:
   const int m_n_proc;

public:
   MatriplexTrackPackerBase(int N_proc) :
      m_n_proc(N_proc)
   {
      assert(m_n_proc <= NN);
   }
};


//==============================================================================
// MatriplexTrackPackerSlurpIn
//==============================================================================

class MatriplexTrackPackerSlurpIn : public MatriplexTrackPackerBase
{
   char *m_base;     // Template to matriplex type, use as scale.
                     // Need typedef in Matrix.h (and pass it in Packer Selection
                     // section below).
   int   m_n_pos;
   int   m_off_error;
   int   m_off_param;

   int   m_idx[NN]  __attribute__((aligned(64)));

public:
   MatriplexTrackPackerSlurpIn(int N_proc) :
      MatriplexTrackPackerBase(N_proc),
      m_n_pos(0)
   {}

   void Reset()
   {
      m_n_pos = 0;
   }

   void AddInput(const Track &track)
   {
      if (m_n_pos == 0)
      {
         m_base      = (char*) &track;
         m_off_error = (char*) track.errors().Array()     - m_base;
         m_off_param = (char*) track.parameters().Array() - m_base;
      }

      m_idx[m_n_pos] = (char*) &track - m_base;

      ++m_n_pos;

      assert(m_n_pos <= m_n_proc);
   }

   void Pack(MPlexLS &err, MPlexLV &par)
   {
#ifdef MIC_INTRINSICS
      __m512i vi = _mm512_load_epi32(m_idx);
      err.SlurpIn(m_base + m_off_error, vi, m_n_proc);
      par.SlurpIn(m_base + m_off_param, vi, m_n_proc);
#else
      err.SlurpIn(m_base + m_off_error, m_idx, m_n_proc);
      par.SlurpIn(m_base + m_off_param, m_idx, m_n_proc);
#endif
   }
};


//==============================================================================
// MatriplexTrackPackerPlexify
//==============================================================================

class MatriplexTrackPackerPlexify : public MatriplexTrackPackerBase
{
public:
   MatriplexTrackPackerPlexify(int N_proc) :
      MatriplexTrackPackerBase(N_proc)
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

using MatriplexTrackPacker = MatriplexTrackPackerSlurpIn;

}

#endif
