#include <TROOT.h>
#include <TFile.h>

const int Mdet = 7, Mlay = 10;
const int Xlay[] = { -1, 5, 4, 5, 4, 7, 10 };

const bool isbrl[] = { 0, 1, 0, 1, 0, 1, 0 };

struct RZBox
{
  int   m_cnt = 0;
  float m_minr, m_maxr;
  float m_minz, m_maxz;

  void  fill(float r, float z)
  {
    if (m_cnt == 0)
    {
      m_minr = m_maxr = r;
      m_minz = m_maxz = z;
    }
    else
    {
      if      (r < m_minr) m_minr = r;
      else if (r > m_maxr) m_maxr = r;
      if      (z < m_minz) m_minz = z;
      else if (z > m_maxz) m_maxz = z;
    }
    ++m_cnt;
  }

  void print(bool nl=true)
  {
    printf("r~[%.4f, %.4f] z~[%.4f, %.4f]%s",
           m_minr, m_maxr, m_minz, m_maxz, nl ? "\n" : "");
  }

  RZBox Extend(float eps=0.05) const
  {
    float dr = eps * (m_maxr - m_minr);
    float dz = eps * (m_maxz - m_minz);

    RZBox e = *this;
    e.m_minr -= dr;  e.m_maxr += dr;
    e.m_minz -= dz;  e.m_maxz += dz;

    return e;
  }

  RZBox Round(double base=100) const
  {
    RZBox e;
    e.m_cnt = m_cnt;

    e.m_minr = std::floor(m_minr * base) / base;
    e.m_maxr = std::ceil (m_maxr * base) / base;
    e.m_minz = std::floor(m_minz * base) / base;
    e.m_maxz = std::ceil (m_maxz * base) / base;

    return e;
  }

  void MergeWith(const RZBox &o)
  {
    if (m_cnt == 0) *this = o;
    if (o.m_cnt == 0) return;

    m_cnt += o.m_cnt;

    m_minr = std::min(m_minr, o.m_minr);
    m_maxr = std::max(m_maxr, o.m_maxr);
    m_minz = std::min(m_minz, o.m_minz);
    m_maxz = std::max(m_maxz, o.m_maxz);
  }
};

struct BBS
{
  int   cnt[Mdet][Mlay] = { 0 };
  RZBox b[Mdet][Mlay];
  RZBox p[Mdet][Mlay];
  RZBox n[Mdet][Mlay];

  RZBox& select_rzbox(int det, int lay, float z)
  {
    if (isbrl[det])  return b[det][lay];
    return (z > 0) ? p[det][lay] : n[det][lay];
  }

  void reset()
  {
    for (int d = 1; d < Mdet; ++d)
    {
      for (int l = 1; l < Mlay; ++l)
      {
        cnt[d][l] = 0;
        b[d][l].m_cnt = 0;
        p[d][l].m_cnt = 0;
        n[d][l].m_cnt = 0;
      }
    }
  }

  void print()
  {
    for (int d = 1; d < Mdet; ++d)
    {
      printf("det %d, is_brl=%d:\n", d, isbrl[d]);
      for (int l = 1; l < Xlay[d]; ++l)
      {
        if (isbrl[d])
        {
          printf("  lay %d: cnt=%5d, ", l, cnt[d][l]); b[d][l].print();
        }
        else
        {
          printf("  pos %d: cnt=%5d, ", l, p[d][l].m_cnt); p[d][l].print();
          printf("  neg %d: cnt=%5d, ", l, n[d][l].m_cnt); n[d][l].print();
        }
      }
      printf("\n");
    }
  }

  void save(const char* fname="bbs.root")
  {
    auto f = TFile::Open(fname, "RECREATE");
    f->WriteObject(this, "bbs");
    f->Close();
    delete f;
  }

  void load(const char* fname="bbs.root")
  {
    auto f = TFile::Open(fname);
    if ( ! f) throw std::runtime_error("tnc::load could not open file");

    BBS *b = 0;
    gFile->GetObject("bbs", b);
    if ( ! b) throw std::runtime_error("tnc::load could not read bbs");

    b->print();
    memcpy(this, b, sizeof(BBS));

    delete b;
    f->Close();
    delete f;
  }
};

void tncdefs()
{
  printf("tncdefs() loaded\n");
}

#pragma link C++ class RZBox;
#pragma link C++ class BBS;
