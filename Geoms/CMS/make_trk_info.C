/*
  Sibling layers are needed for seeding from MC tracks.
  As we do not do seeding we also do not fill them here.
*/

/*
  Some "leftovers" of CMSSW from Config.h

  constexpr float cmsAvgRads[13] = {4.42,7.31,10.17,25.58,33.98,41.79,49.78,60.78,69.2,77.96,86.80,96.53,108.00}; // cms average radii, noSplit version
  constexpr float cmsAvgRads[17] = {4.42, 7.31, 10.17,
                                    25.58,25.58, 33.98,33.98, 41.79, 49.78,
                                    60.57,61.00, 69.41,68.98, 77.96, 86.80, 96.53, 108.00}; // cms average radii, split version
  constexpr float cmsDeltaRad = 2.5; //fixme! using constant 2.5 cm, to be taken from layer properties

  constexpr float cmsAvgZs[26]         = {35.5,48.5,  79.8,79.8,92.6,92.6,105.6,105.6,  131.3,131.3,145.3,145.3,159.3,159.3,173.9,173.9,187.8,187.8,205.4,205.4,224.0,224.0,244.4,244.4,266.3,266.3}; // cms average z
  constexpr float cmsDiskMinRs[26]     = { 5.7, 5.7,  23.1,22.8,23.1,22.8, 23.1, 22.8,  23.3, 23.0, 23.3, 23.0, 23.3, 23.0, 31.6, 34.4, 31.6, 34.4, 31.6, 34.4, 59.9, 38.8, 59.9, 38.8, 59.9, 49.9};
  constexpr float cmsDiskMaxRs[26]     = {14.7,14.7,  50.8,42.0,50.8,42.0, 50.8, 42.0,  76.1,110.0, 76.1,110.0, 76.1,110.0, 75.9,109.7, 75.9,109.7, 75.9,109.7, 75.9,109.4, 75.9,109.4, 75.9,109.4};
  constexpr float cmsDiskMinRsHole[26] = { 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,   42.0,  0.0, 42.0,  0.0, 42.0,  0.0, 42.1,  0.0, 42.1,  0.0, 42.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};
  constexpr float cmsDiskMaxRsHole[26] = {999.,999.,  999.,999.,999.,999., 999., 999.,  59.9, 999., 59.9, 999., 59.9, 999., 59.7, 999., 59.7, 999., 59.7, 999., 999., 999., 999., 999., 999., 999.};

  const float g_disk_dr[]      = {   1,   1,  20,  20,  20,  20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20,   20};
  const float g_layer_zwidth[] = { 30, 30, 30, 70, 70, 70, 70, 70, 70, 110, 110, 110, 110, 110, 110, 110, 110 };
  const float g_layer_dz[]     = { 1, 1, 1, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20 };
*/

#include "tncdefs.h"

#include <cstdio>

BBS   BBB;
FILE *OOO = 0;
int   DEP = 0;

#define NLN      fprintf(OOO,"\n")
#define PRN(...) fprintf(OOO,"%*s",DEP,""); fprintf(OOO, __VA_ARGS__); fprintf(OOO,"\n");

#define SCOPE_BEG PRN("{"); DEP += 2
#define SCOPE_END DEP -= 2; PRN("}")

#define ASSF(var, val) PRN("li.%s = %.3f;", #var, val)
#define ASSI(var, val) PRN("li.%s = %d;", #var, val)
#define ASSB(var, val) PRN("li.%s = %s;", #var, val ? "true" : "false")

void add_barrel(int &lid, int det, int lay, bool is_pix,
                int necp, int necn)
{
  RZBox b = BBB.b[det][lay].Round(100);

  SCOPE_BEG;

  PRN("LayerInfo & li  = ti.m_layers[%d];", lid);
  PRN("li.m_layer_type = LayerInfo::Barrel;");
  PRN("li.set_limits(%.3f, %.3f, %.3f, %.3f);", b.m_minr, b.m_maxr, b.m_minz, b.m_maxz);
  PRN("li.m_propagate_to = li.r_mean();");
  PRN("li.set_next_layers(%d, %d, %d);", lid < 17 ? lid + 1 : -1, necp, necn);
  ASSB(m_is_outer, lid == 17);
  if (is_pix)
  {
    ASSF(m_q_bin, 2.0);
    PRN("li.set_selection_limits(0.01, 0.05, 1.0, 2.0);");
  }
  else
  {
    ASSF(m_q_bin, 20.0);
    PRN("li.set_selection_limits(0.01, 0.2, 10.0, 20.0);");
  }
      
  SCOPE_END;

  ++lid;
}

void add_ecap(int &lid, int det, int lay, bool is_pix,
              int nbrl)
{
  int lid_store = lid;
  
  SCOPE_BEG;
  {
    RZBox b = BBB.p[det][lay].Round(100);

    PRN("LayerInfo & li  = ti.m_layers[%d];", lid);
    PRN("li.m_layer_type = LayerInfo::EndCapPos;");
    PRN("li.set_limits(%.3f, %.3f, %.3f, %.3f);", b.m_minr, b.m_maxr, b.m_minz, b.m_maxz);
    PRN("li.m_propagate_to = li.z_mean();");
    PRN("li.set_next_layers(%d, %d, %d);", nbrl, lid < 44 ? lid + 1 : -1, -1);
    ASSB(m_is_outer, lid == 44);
    if (is_pix) {
      ASSF(m_q_bin, 1.0);
      PRN("li.set_selection_limits(0.01, 0.05, 0.8, 1.6);");
    } else {
      ASSF(m_q_bin, 20.0);
      PRN("li.set_selection_limits(0.01, 0.2, 10.0, 20.0);");
    }
  }
  SCOPE_END;
  SCOPE_BEG;
  {
    lid += 27;
    RZBox b = BBB.n[det][lay].Round(100);

    PRN("LayerInfo & li  = ti.m_layers[%d];", lid);
    PRN("li.m_layer_type = LayerInfo::EndCapNeg;");
    PRN("li.set_limits(%.3f, %.3f, %.3f, %.3f);", b.m_minr, b.m_maxr, b.m_minz, b.m_maxz);
    PRN("li.m_propagate_to = li.z_mean();");
    PRN("li.set_next_layers(%d, %d, %d);", nbrl, lid < 71 ? lid + 1 : -1, -1);
    ASSB(m_is_outer, lid == 71);
    if (is_pix) {
      ASSF(m_q_bin, 1.0);
      PRN("li.set_selection_limits(0.01, 0.05, 0.8, 1.6);");
    } else {
      ASSF(m_q_bin, 20.0);
      PRN("li.set_selection_limits(0.01, 0.2, 10.0, 20.0);");
    }
  }
  SCOPE_END;

  lid = lid_store + 1;
}

void print_trk_info()
{
  OOO = fopen("../CMS-2017.acc", "w");
  // OOO = stdout;

  PRN("void Create_CMS_2017_AutoGen(TrackerInfo& ti)");
  SCOPE_BEG;

  int lid = 0;

  PRN("// PIXB\n");
  add_barrel(lid, 1, 1, true, 18, 45); NLN;
  add_barrel(lid, 1, 2, true, 18, 45); NLN;
  add_barrel(lid, 1, 3, true, 18, 45); NLN;
  add_barrel(lid, 1, 4, true, 21, 48); NLN;
  NLN;
  
  PRN("// TIB\n");
  add_barrel(lid, 3, 1, false, 21, 48); NLN;
  add_barrel(lid, 3, 1, false, 21, 48); NLN;
  add_barrel(lid, 3, 2, false, 21, 48); NLN;
  add_barrel(lid, 3, 2, false, 21, 48); NLN;
  add_barrel(lid, 3, 3, false, 21, 48); NLN;
  add_barrel(lid, 3, 4, false, 27, 54); NLN;
  NLN;

  PRN("// TOB\n");
  add_barrel(lid, 5, 1, false, 27, 54); NLN;
  add_barrel(lid, 5, 1, false, 27, 54); NLN;
  add_barrel(lid, 5, 2, false, 27, 54); NLN;
  add_barrel(lid, 5, 2, false, 27, 54); NLN;
  add_barrel(lid, 5, 3, false, 27, 54); NLN;
  add_barrel(lid, 5, 4, false, 27, 54); NLN;
  add_barrel(lid, 5, 5, false, 27, 54); NLN;
  add_barrel(lid, 5, 6, false, -1, -1); NLN;
  NLN;

  PRN("// PIXE +/-\n");
  add_ecap(lid, 2, 1, true, 4);
  add_ecap(lid, 2, 2, true, 4);
  add_ecap(lid, 2, 3, true, 4);
  NLN;

  PRN("// TID +/-\n");
  add_ecap(lid, 4, 1, false, 10);
  add_ecap(lid, 4, 1, false, 10);
  add_ecap(lid, 4, 2, false, 10);
  add_ecap(lid, 4, 2, false, 10);
  add_ecap(lid, 4, 3, false, 10);
  add_ecap(lid, 4, 3, false, 10);
  NLN;

  PRN("// TID +/-\n");
  add_ecap(lid, 6, 1, false, -1);
  add_ecap(lid, 6, 1, false, -1);
  add_ecap(lid, 6, 2, false, -1);
  add_ecap(lid, 6, 2, false, -1);
  add_ecap(lid, 6, 3, false, -1);
  add_ecap(lid, 6, 3, false, -1);
  add_ecap(lid, 6, 4, false, -1);
  add_ecap(lid, 6, 4, false, -1);
  add_ecap(lid, 6, 5, false, -1);
  add_ecap(lid, 6, 5, false, -1);
  add_ecap(lid, 6, 6, false, -1);
  add_ecap(lid, 6, 6, false, -1);
  add_ecap(lid, 6, 7, false, -1);
  add_ecap(lid, 6, 7, false, -1);
  add_ecap(lid, 6, 8, false, -1);
  add_ecap(lid, 6, 8, false, -1);
  add_ecap(lid, 6, 9, false, -1);
  add_ecap(lid, 6, 9, false, -1);
  NLN;

  SCOPE_END;

  if (OOO != stdout)
  {
    fclose(OOO);
    OOO = 0;
  }
}

BBS* make_trk_info()
{
  BBB.load();

  printf("make_trk_info() and BBS loaded\n");

  print_trk_info();

  return &BBB;
}
