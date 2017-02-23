#include "CylCowWLids.h"

#include "TrackerInfo.h"

#include <cmath>

namespace
{
  float getTheta(float r, float z)
  {
    return std::atan2(r,z);
  }

  float getEta(float r, float z)
  {
    return -1.0f * std::log( std::tan(getTheta(r,z)/2.0f) );
  }

  float getEta(float theta)
  {
    return -1.0f * std::log( std::tan(theta/2.0f) );
  }

  float getTgTheta(float eta)
  {
    return std::tan(2.0*std::atan(std::exp(-eta)));
  }

  class CylCowWLidsCreator
  {
    TrackerInfo &m_trkinfo;

    float        m_det_half_thickness = 0.005; // for 100 micron total

    int          m_n_barrel   = 10;
    int          m_n_endcap   =  9;
    int          m_first_ecap = 10;

    //------------------------------------------------------------------------------

    void add_barrel(int lid, float r, float z, float eta)
    {
      // printf("Adding barrel layer r=%.3f z=%.3f eta_t=%.3f\n", r, z, eta);

      LayerInfo & li  = m_trkinfo.m_layers[lid];

      li.m_rin  = r - m_det_half_thickness;
      li.m_rout = r + m_det_half_thickness;
      li.m_zmin = -z;
      li.m_zmax =  z;

      li.m_next_barrel   = lid < 9 ? lid + 1 : -1;
      li.m_next_ecap_pos = lid < 9 ? lid + 1 +  9 : -1;
      li.m_next_ecap_neg = lid < 9 ? lid + 1 + 18 : -1;

      li.m_is_barrel = true;
      li.m_is_outer  = (lid == 9);
    }

    void add_barrel_r_eta(int lid, float r, float eta)
    {
      float z = r / getTgTheta(eta);

      add_barrel(lid, r, z, eta);
    }

    void add_barrel_r_z(int lid, float r, float z)
    {
      float eta = getEta(r, z);

      add_barrel(lid, r, z, eta);
    }

    void add_endcap(int lid, float r, float z, float eta)
    {
      float r_end = z * getTgTheta(eta);
  
      // printf("Adding endcap layer r=%.3f z=%.3f r_l=%.3f eta_l=%.3f\n", r, z, r_end, eta);

      {
        LayerInfo & li  = m_trkinfo.m_layers[lid];

        li.m_rin  = r_end;
        li.m_rout = r;
        li.m_zmin = z - m_det_half_thickness;
        li.m_zmax = z + m_det_half_thickness;

        li.m_next_barrel   = lid < 18 ? lid - 10 + 2 : -1;
        li.m_next_ecap_pos = lid < 18 ? lid + 1 : -1;
        li.m_next_ecap_neg = -1;

        li.m_is_barrel = false;
        li.m_is_outer  = (lid == 18);
      }
      {
        lid += 9;
        LayerInfo & li  = m_trkinfo.m_layers[lid];

        li.m_rin  = r_end;
        li.m_rout = r;
        li.m_zmin = -z - m_det_half_thickness;
        li.m_zmax = -z + m_det_half_thickness;

        li.m_next_barrel   = lid < 27 ? lid - 19 + 2 : -1;
        li.m_next_ecap_pos = -1;
        li.m_next_ecap_neg = lid < 27 ? lid + 1 : -1;

        li.m_is_barrel = false;
        li.m_is_outer  = (lid == 27);
      }
    }

    //------------------------------------------------------------------------------

  public:

    CylCowWLidsCreator(TrackerInfo& ti) :
      m_trkinfo(ti)
    {}

    void FillTrackerInfo()
    {
      // Actual coverage for tracks with z = 3cm is 2.4
      float full_eta = 2.5;
      float full_eta_pix_0 = 2.55; // To account for BS z-spread
      float full_eta_ec_in[] = { 0, 2.525, 2.515 };

      float pix_0  =  4, pix_sep   = 6;
      float pix_z0 = 24, pix_zgrow = 6;

      float sct_sep = 10;
      float sct_0   = pix_0 + 2 * pix_sep   + sct_sep;
      float sct_zgrow = 10;
      float sct_z0 = pix_z0 + 2 * pix_zgrow + sct_zgrow;

      float pix_ec_zgap   = 2;
      float pix_ec_rextra = 2;

      float sct_ec_zgap   = 4;
      float sct_ec_rextra = 4;

      for (int i = 0; i < 10; ++i) m_trkinfo.new_barrel_layer();
      for (int i = 0; i <  9; ++i) m_trkinfo.new_ecap_pos_layer();
      for (int i = 0; i <  9; ++i) m_trkinfo.new_ecap_neg_layer();
      
      add_barrel_r_eta(0, pix_0, full_eta_pix_0);

      add_barrel_r_z  (1, pix_0 + 1 * pix_sep, pix_z0 + 1 * pix_zgrow);
      add_barrel_r_z  (2, pix_0 + 2 * pix_sep, pix_z0 + 2 * pix_zgrow);

      for (int i = 0; i < 7; ++i)
      {
        add_barrel_r_z(3 + i, sct_0  + i * sct_sep, sct_z0 + i * sct_zgrow);
      }

      for (int i = 1; i < 3; ++i)
      {
        add_endcap(9 + i,
                   pix_0  + i * pix_sep   + pix_ec_rextra,
                   pix_z0 + i * pix_zgrow + pix_ec_zgap,
                   full_eta_ec_in[i]);
      }
      for (int i = 0; i < 7; ++i)
      {
        add_endcap(12 + i,
                   sct_0  + i * sct_sep   + sct_ec_rextra,
                   sct_z0 + i * sct_zgrow + sct_ec_zgap,
                   full_eta);
      }
      // + endcap disks at -z
    }
  };
}

//==============================================================================

void Create_TrackerInfo(TrackerInfo& ti, bool verbose)
{
  CylCowWLidsCreator creator(ti);

  creator.FillTrackerInfo();

  if (verbose) {
    printf("==========================================================================================\n");
  }
  printf("CylCowWLids -- Create_TrackerInfo finished\n");

  if (verbose) {
    printf("==========================================================================================\n");
    for (auto &i : ti.m_layers)  i.print_layer();
    printf("==========================================================================================\n");
  }
}
