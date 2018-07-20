/*
  Sibling layers are needed for seeding from MC tracks.
  As we do not do seeding we also do not fill them here.
*/

/*
  Some "leftovers" of CMSSW from Config.h. For 2016 layout!

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


// For TID and TEC only:
constexpr float cmsDiskMinRs[24]     = {23.1,22.8,23.1,22.8,23.1,22.8,    23.3, 23.0, 23.3, 23.0, 23.3, 23.0, 31.6, 34.4, 31.6, 34.4, 31.6, 34.4,  59.9, 38.8, 59.9, 38.8, 59.9, 49.9};
constexpr float cmsDiskMaxRs[24]     = {50.8,42.0,50.8,42.0,50.8,42.0,    76.1,110.0, 76.1,110.0, 76.1,110.0, 75.9,109.7, 75.9,109.7, 75.9,109.7,  75.9,109.4, 75.9,109.4, 75.9,109.4};
constexpr float cmsDiskMinRsHole[24] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,    42.0,  0.0, 42.0,  0.0, 42.0,  0.0, 42.1,  0.0, 42.1,  0.0, 42.1,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0};
constexpr float cmsDiskMaxRsHole[24] = {999.,999.,999.,999.,999.,999.,    59.9, 999., 59.9, 999., 59.9, 999., 59.7, 999., 59.7, 999., 59.7, 999.,  999., 999., 999., 999., 999., 999.};
// 1. Note -- TEC stereo has index pefore the correcponding mono / full sub-layer!
// 2. For TID we only reduce the stereo layer R limits.
// 3. For TEC we also store the hole extent.
// 4. Other disks limits are taken from the tracking ntuple bounding boxes.
//    It should be possible to make this automatic.

constexpr int   N_barrel     = 18;
constexpr int   N_endcap     = 27;
constexpr int   N_pix_endcap =  3;

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

void setHitSelDynamicFactors(int id){

  // config on hit selection windows
  // geometry and track boundaries
  int brl_tibId[2] = {4,9};
  int brl_tobId[2] = {10,17};
  int ecp_stripId[2] = {21,44};
  int ecn_stripId[2] = {48,71};
  // phi dynamic factor
  float phif_ptlow_brl_mono       = 3.0;
  float phif_ptlow_brl_stereo     = 2.0;
  float phif_ptlow_treg_ec_mono   = 3.0;
  float phif_ptlow_treg_ec_stereo = 2.0;
  float phif_ptlow_ec_mono        = 2.0;
  float phif_treg_ec_mono         = 1.5;
  // q dynamic factors
  float qf_treg_tib = 1.5;
  float qf_treg_tob = 1.25;
  
  bool is_stereo_lyr = ( (id<21 && (id==5 || id==7 || id==11 || id==13)) || (id>=21 && id<45 && id%2==0) || (id>=48 && id<72 && id%2>0) ); // In ECP(N), even (odd) layers are stereo
  
  if( id>=brl_tibId[0] && id<=brl_tobId[1] )
    {
    
      if( id<brl_tobId[0] ){
	ASSF(m_qf_treg, qf_treg_tib);
      }
      else{
	ASSF(m_qf_treg, qf_treg_tob);
      }

      if( is_stereo_lyr ){
	ASSF(m_phif_lpt_brl, phif_ptlow_brl_stereo);
      }
      else{
	ASSF(m_phif_lpt_brl, phif_ptlow_brl_mono);
      }

    }
  else if( (id>=ecp_stripId[0] && id<=ecp_stripId[1]) || (id>=ecn_stripId[0] && id<=ecn_stripId[1]) )
    {
      
      if( is_stereo_lyr ){
	ASSF(m_phif_lpt_treg, phif_ptlow_treg_ec_stereo);
      }
      else{
	ASSF(m_phif_lpt_treg, phif_ptlow_treg_ec_mono);
	ASSF(m_phif_lpt_ec, phif_ptlow_ec_mono);
	ASSF(m_phif_treg, phif_treg_ec_mono);
      }

    }
   
} 


void assignSubDetector(int id){

  int pixb_lyrs[2]={0,3};
  int tib_lyrs[2]={4,9};
  int tob_lyrs[2]={10,17};
  int pixep_lyrs[2]={18,20};
  int pixen_lyrs[2]={45,47};
  int tidp_lyrs[2]={21,26};
  int tidn_lyrs[2]={48,53};
  int tecp_lyrs[2]={27,44};
  int tecn_lyrs[2]={54,71};

  if(id>=pixb_lyrs[0] && id<=pixb_lyrs[1]){ 
    ASSB(m_is_pixb_lyr, 1);
  }
  else if(id>=tib_lyrs[0] && id<=tib_lyrs[1]){ 
    ASSB(m_is_tib_lyr, 1);
  }
  else if(id>=tob_lyrs[0] && id<=tob_lyrs[1]){
    ASSB(m_is_tib_lyr, 1);
  }  
  else if((id>=pixep_lyrs[0] && id<=pixep_lyrs[1])||(id>=pixen_lyrs[0] && id<=pixen_lyrs[1])){
    ASSB(m_is_pixe_lyr, 1);
  }
  else if((id>=tidp_lyrs[0] && id<=tidp_lyrs[1])||(id>=tidn_lyrs[0] && id<=tidn_lyrs[1])){
    ASSB(m_is_tid_lyr, 1);
  }
  else if((id>=tecp_lyrs[0] && id<=tecp_lyrs[1])||(id>=tecn_lyrs[0] && id<=tecn_lyrs[1])){
    ASSB(m_is_tec_lyr, 1);
  }

}

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
  ASSB(m_is_stereo_lyr, ( (lid<21 && (lid==5 || lid==7 || lid==11 || lid==13)) || (lid>=21 && lid<45 && lid%2==0) || (lid>=48 && lid<72 && lid%2>0) )); // In ECP(N), even (odd) layers are stereo
  assignSubDetector(lid);
  setHitSelDynamicFactors(lid);
  if (is_pix)
  {
    ASSB(m_is_seed_lyr, 1);
    ASSF(m_q_bin, 2.0);
    PRN("li.set_selection_limits(0.01, 0.05, 1.0, 2.0);");
  }
  else if (lid==4)
  {
    ASSF(m_q_bin, 6.0);
    PRN("li.set_selection_limits(0.01, 0.015, 6.0, 12.0);");
  }
  else if (lid==5)
  {
    ASSF(m_q_bin, 6.0);
    PRN("li.set_selection_limits(0.023, 0.03, 6.0, 12.0);");
  }
  else if (lid==6)
  {
    ASSF(m_q_bin, 6.0);
    PRN("li.set_selection_limits(0.01, 0.015, 6.0, 12.0);");
  }
  else if (lid==7)
  {
    ASSF(m_q_bin, 6.0);
    PRN("li.set_selection_limits(0.016, 0.03, 6.0, 12.0);");
  }
  else if (lid==8)
  {
    ASSF(m_q_bin, 6.0);
    PRN("li.set_selection_limits(0.01, 0.015, 6.0, 12.0);");
  }
  else if (lid==9)
  {
    ASSF(m_q_bin, 6.0);
    PRN("li.set_selection_limits(0.01, 0.015, 6.0, 12.0);");
  }
  else if (lid==10)
  {
    ASSF(m_q_bin, 9.5);
    PRN("li.set_selection_limits(0.01, 0.015, 9.5, 19.0);");
  }
  else if (lid==11)
  {
    ASSF(m_q_bin, 9.5);
    PRN("li.set_selection_limits(0.016, 0.03, 9.5, 19.0);");
  }
  else if (lid==12)
  {
    ASSF(m_q_bin, 9.5);
    PRN("li.set_selection_limits(0.01, 0.015, 9.5, 19.0);");
  }
  else if (lid==13)
  {
    ASSF(m_q_bin, 9.5);
    PRN("li.set_selection_limits(0.013, 0.03, 9.5, 19.0);");
  }
  else if (lid==14)
  {
    ASSF(m_q_bin, 9.5);
    PRN("li.set_selection_limits(0.01, 0.015, 9.5, 19.0);");
  }
  else if (lid==15)
  {
    ASSF(m_q_bin, 9.5);
    PRN("li.set_selection_limits(0.01, 0.015, 9.5, 19.0);");
  }
  else if (lid==16)
  {
    ASSF(m_q_bin, 9.5);
    PRN("li.set_selection_limits(0.01, 0.015, 9.5, 19.0);");
  }
  else if (lid==17)
  {
    ASSF(m_q_bin, 9.5);
    PRN("li.set_selection_limits(0.01, 0.015, 9.5, 19.0);");
  }
  else
  {
    ASSF(m_q_bin, 20.0);
    PRN("li.set_selection_limits(0.01, 0.2, 10.0, 20.0);");
  }
      
  SCOPE_END;

  ++lid;
}

void add_ecap(int &lid, int det, int lay, bool is_pix, bool stereo_hack,
              int next_brl)
{
  int lid_store = lid;

  bool hole_hack = false;
  int shi = -1;
  if (stereo_hack)
  {
    shi = lid - N_barrel - N_pix_endcap;
    if (shi > 5)
    {
      shi = shi - 1; // TEC order reversed in Giuseppe's map!
      if (shi < 18)  // No holes on last three layers
      {
        hole_hack = true;
      }
    }
  }

  SCOPE_BEG;
  {
    RZBox b = BBB.p[det][lay].Round(100);

    float min_r = stereo_hack ? cmsDiskMinRs[shi] : b.m_minr;
    float max_r = stereo_hack ? cmsDiskMaxRs[shi] : b.m_maxr;

    PRN("LayerInfo & li  = ti.m_layers[%d];", lid);
    PRN("li.m_layer_type = LayerInfo::EndCapPos;");
    PRN("li.set_limits(%.3f, %.3f, %.3f, %.3f);", min_r, max_r, b.m_minz, b.m_maxz);
    PRN("li.m_propagate_to = li.z_mean();");
    PRN("li.set_next_layers(%d, %d, %d);", next_brl, lid < 44 ? lid + 1 : -1, -1);
    ASSB(m_is_outer, lid == 44);
    ASSB(m_is_stereo_lyr, ( (lid<21 && (lid==5 || lid==7 || lid==11 || lid==13)) || (lid>=21 && lid<45 && lid%2==0) || (lid>=48 && lid<72 && lid%2>0) )); // In ECP(N), even (odd) layers are stereo
    assignSubDetector(lid);
    setHitSelDynamicFactors(lid);
    if (is_pix) {
      ASSB(m_is_seed_lyr, 1);
      ASSF(m_q_bin, 1.0);
      PRN("li.set_selection_limits(0.01, 0.015, 0.8, 1.6);");
    } 
    else {
      if(lid==21){
	ASSF(m_q_bin, 5.5);
	PRN("li.set_selection_limits(0.01, 0.015, 5.5, 11.0);");
      }
      else if(lid==22){
	ASSF(m_q_bin, 5.5);
	PRN("li.set_selection_limits(0.01, 0.03, 5.5, 11.0);");
      }
      else if(lid==23){
	ASSF(m_q_bin, 5.5);
	PRN("li.set_selection_limits(0.01, 0.015, 5.5, 11.0);");
      }
      else if(lid==24){
	ASSF(m_q_bin, 5.5);
	PRN("li.set_selection_limits(0.01, 0.03, 5.5, 11.0);");
      }
      else if(lid==25){
	ASSF(m_q_bin, 5.5);
	PRN("li.set_selection_limits(0.01, 0.015, 5.5, 11.0);");
      }
      else if(lid==26){
	ASSF(m_q_bin, 5.5);
	PRN("li.set_selection_limits(0.01, 0.03, 5.5, 11.0);");
      }
      else if(lid==27){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==28){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==29){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==30){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==31){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==32){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==33){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==34){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==35){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==36){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==37){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==38){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==39){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==40){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==41){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==42){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==43){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==44){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else{
	ASSF(m_q_bin, 20.0);
	PRN("li.set_selection_limits(0.01, 0.2, 10.0, 20.0);");
      }
    }
    if (hole_hack) {
      PRN("li.set_r_hole_range(%.3f, %.3f);", cmsDiskMinRsHole[shi], cmsDiskMaxRsHole[shi]);
    }
  }
  SCOPE_END;
  SCOPE_BEG;
  {
    lid += N_endcap;
    RZBox b = BBB.n[det][lay].Round(100);

    float min_r = stereo_hack ? cmsDiskMinRs[shi] : b.m_minr;
    float max_r = stereo_hack ? cmsDiskMaxRs[shi] : b.m_maxr;

    PRN("LayerInfo & li  = ti.m_layers[%d];", lid);
    PRN("li.m_layer_type = LayerInfo::EndCapNeg;");
    PRN("li.set_limits(%.3f, %.3f, %.3f, %.3f);", min_r, max_r, b.m_minz, b.m_maxz);
    PRN("li.m_propagate_to = li.z_mean();");
    PRN("li.set_next_layers(%d, %d, %d);", next_brl, -1, lid < 71 ? lid + 1 : -1);
    ASSB(m_is_outer, lid == 71);
    ASSB(m_is_stereo_lyr, ( (lid<21 && (lid==5 || lid==7 || lid==11 || lid==13)) || (lid>=21 && lid<45 && lid%2==0) || (lid>=48 && lid<72 && lid%2>0) )); // In ECP(N), even (odd) layers are stereo
    assignSubDetector(lid);
    setHitSelDynamicFactors(lid);
    if (is_pix) {
      ASSB(m_is_seed_lyr, 1);
      ASSF(m_q_bin, 1.0);
      PRN("li.set_selection_limits(0.01, 0.015, 0.8, 1.6);");
    } 
    else {
      if(lid==48){
	ASSF(m_q_bin, 5.5);
	PRN("li.set_selection_limits(0.01, 0.015, 5.5, 11.0);");
      }
      else if(lid==49){
	ASSF(m_q_bin, 5.5);
	PRN("li.set_selection_limits(0.01, 0.03, 5.5, 11.0);");
      }
      else if(lid==50){
	ASSF(m_q_bin, 5.5);
	PRN("li.set_selection_limits(0.01, 0.015, 5.5, 11.0);");
      }
      else if(lid==51){
	ASSF(m_q_bin, 5.5);
	PRN("li.set_selection_limits(0.01, 0.03, 5.5, 11.0);");
      }
      else if(lid==52){
	ASSF(m_q_bin, 5.5);
	PRN("li.set_selection_limits(0.01, 0.015, 5.5, 11.0);");
      }
      else if(lid==53){
	ASSF(m_q_bin, 5.5);
	PRN("li.set_selection_limits(0.01, 0.03, 5.5, 11.0);");
      }
      else if(lid==54){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==55){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==56){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==57){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==58){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==59){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==60){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==61){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==62){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==63){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==64){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==65){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==66){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==67){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==68){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==69){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else if(lid==70){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.015, 10.0, 20.0);");
      }
      else if(lid==71){
	ASSF(m_q_bin, 10.0);
	PRN("li.set_selection_limits(0.01, 0.03, 10.0, 20.0);");
      }
      else{
	ASSF(m_q_bin, 20.0);
	PRN("li.set_selection_limits(0.01, 0.2, 10.0, 20.0);");
      }
    }
    if (hole_hack) {
      PRN("li.set_r_hole_range(%.3f, %.3f);", cmsDiskMinRsHole[shi], cmsDiskMaxRsHole[shi]);
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
  add_ecap(lid, 2, 1, true, false, 4);
  add_ecap(lid, 2, 2, true, false, 4);
  add_ecap(lid, 2, 3, true, false, 4);
  NLN;

  PRN("// TID +/-\n");
  add_ecap(lid, 4, 1, false, false, 10);
  add_ecap(lid, 4, 1, false, true,  10);
  add_ecap(lid, 4, 2, false, false, 10);
  add_ecap(lid, 4, 2, false, true,  10);
  add_ecap(lid, 4, 3, false, false, 10);
  add_ecap(lid, 4, 3, false, true,  10);
  NLN;

  PRN("// TEC +/-\n");
  add_ecap(lid, 6, 1, false, false, -1);
  add_ecap(lid, 6, 1, false, true,  -1);
  add_ecap(lid, 6, 2, false, false, -1);
  add_ecap(lid, 6, 2, false, true,  -1);
  add_ecap(lid, 6, 3, false, false, -1);
  add_ecap(lid, 6, 3, false, true,  -1);
  add_ecap(lid, 6, 4, false, false, -1);
  add_ecap(lid, 6, 4, false, true,  -1);
  add_ecap(lid, 6, 5, false, false, -1);
  add_ecap(lid, 6, 5, false, true,  -1);
  add_ecap(lid, 6, 6, false, false, -1);
  add_ecap(lid, 6, 6, false, true,  -1);
  add_ecap(lid, 6, 7, false, false, -1);
  add_ecap(lid, 6, 7, false, true,  -1);
  add_ecap(lid, 6, 8, false, false, -1);
  add_ecap(lid, 6, 8, false, true,  -1);
  add_ecap(lid, 6, 9, false, false, -1);
  add_ecap(lid, 6, 9, false, true,  -1);
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
