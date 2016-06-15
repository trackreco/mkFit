#include "seedtestMPlex.h"
#include "seedtest.h"
#include "Hit.h"

// #define DEBUG
#include "Debug.h"

void findSeedsByRoadSearch(TrackVec& evt_seed_tracks, std::vector<LayerOfHits>& evt_lay_hits, Event *& ev)
{
  bool debug(false);

  // use this to initialize tracks
  TrackState dummystate;
  for (int i = 0; i < 6; i++) {
    dummystate.parameters[i] = 0.f;
  }
  dummystate.errors = ROOT::Math::SMatrixIdentity();

  // 0 = first layer, 1 = second layer, 2 = third layer
  const LayerOfHits& lay1_hits = evt_lay_hits[1];
  LayerOfHits& lay0_hits = evt_lay_hits[0];
  LayerOfHits& lay2_hits = evt_lay_hits[2];

  int seedID = 0;
  for (int ihit1 = 0; ihit1 < lay1_hits.m_capacity; ++ihit1) { 
    dprint("seedID: " << seedID);

    const Hit & hit1   = lay1_hits.m_hits[ihit1];
    const float hit1_z = hit1.z();

    dprint(" ihit1: " << ihit1 << " mcTrackID: " << hit1.mcTrackID(ev->simHitsInfo_) << " phi: " << hit1.phi() << " z: " << hit1.z());
    dprint("  predphi: " << hit1.phi() << "+/-" << Config::lay01angdiff << " predz: " << hit1.z()/2.0f << "+/-" << Config::seed_z0cut/2.0f << std::endl);

    std::vector<int> cand_hit0_indices; // pass by reference
    lay0_hits.SelectHitIndices(hit1_z/2.0f,hit1.phi(),Config::seed_z0cut/2.0f,Config::lay01angdiff,cand_hit0_indices,true,false);
    for (auto&& ihit0 : cand_hit0_indices){
      const Hit & hit0   = lay0_hits.m_hits[ihit0]; 
      const float hit0_z = hit0.z();
      const float hit0_x = hit0.x(); const float hit0_y = hit0.y();
      const float hit1_x = hit1.x(); const float hit1_y = hit1.y();
      const float hit01_r2 = getRad2(hit0_x-hit1_x,hit0_y-hit1_y);

      const float quad = std::sqrt((4.0f*Config::maxCurvR*Config::maxCurvR - hit01_r2) / hit01_r2);
    
      // center of negative curved track
      const float aneg = 0.5f*((hit0_x+hit1_x)-(hit0_y-hit1_y)*quad);
      const float bneg = 0.5f*((hit0_y+hit1_y)+(hit0_x-hit1_x)*quad);

      // negative points of intersection with third layer
      float lay2_negx = 0.0f, lay2_negy = 0.0f;
      intersectThirdLayer(aneg,bneg,hit1_x,hit1_y,lay2_negx,lay2_negy);
      const float lay2_negphi = getPhi(lay2_negx,lay2_negy);

      // center of positive curved track
      const float apos = 0.5f*((hit0_x+hit1_x)+(hit0_y-hit1_y)*quad);
      const float bpos = 0.5f*((hit0_y+hit1_y)-(hit0_x-hit1_x)*quad);
      
      // positive points of intersection with third layer
      float lay2_posx = 0.0f, lay2_posy = 0.0f;
      intersectThirdLayer(apos,bpos,hit1_x,hit1_y,lay2_posx,lay2_posy);
      const float lay2_posphi = getPhi(lay2_posx,lay2_posy);

      std::vector<int> cand_hit2_indices;
      lay2_hits.SelectHitIndices((2.0f*hit1_z-hit0_z),(lay2_posphi+lay2_negphi)/2.0f,
				 Config::seed_z2cut,(lay2_posphi-lay2_negphi)/2.0f,
				 cand_hit2_indices,true,false);

      dprint(" ihit0: " << ihit0 << " mcTrackID: " << hit0.mcTrackID(ev->simHitsInfo_) << " phi: " << hit0.phi() << " z: " << hit0.z());
      dprint("  predphi: " << (lay2_posphi+lay2_negphi)/2.0f << "+/-" << (lay2_posphi-lay2_negphi)/2.0f << " predz: " << 2.0f*hit1_z-hit0_z << "+/-" << Config::seed_z2cut << std::endl);

      for (auto&& ihit2 : cand_hit2_indices){ // loop over candidate second layer hits
	const Hit & hit2 = lay2_hits.m_hits[ihit2];

	const float lay1_predz = (hit0_z + hit2.z()) / 2.0f;
	// filter by residual of second layer hit
	if (std::abs(lay1_predz-hit1_z) > Config::seed_z1cut) continue;

	const float hit2_x = hit2.x(); const float hit2_y = hit2.y();

	// now fit a circle, extract pT and d0 from center and radius
	const float mr = (hit1_y-hit0_y)/(hit1_x-hit0_x);
	const float mt = (hit2_y-hit1_y)/(hit2_x-hit1_x);
	const float a  = (mr*mt*(hit2_y-hit0_y) + mr*(hit1_x+hit2_x) - mt*(hit0_x+hit1_x))/(2.0f*(mr-mt));
	const float b  = -1.0f*(a-(hit0_x+hit1_x)/2.0f)/mr + (hit0_y+hit1_y)/2.0f;
	const float r  = getHypot(hit0_x-a,hit0_y-b);

	// filter by d0 cut 5mm, pT cut 0.5 GeV (radius of 0.5 GeV track)
	if ((r < Config::maxCurvR) || (std::abs(getHypot(a,b)-r) > Config::seed_d0cut)) continue; 
	
	// create a track object
	int hitIndices[3] = {ihit0,ihit1,ihit2};
	Track seed(dummystate,0.0f,seedID,Config::nlayers_per_seed,hitIndices); // argh! super not type-safe with dummystate

	dprint(" ihit2: " << ihit2 << " mcTrackID: " << hit2.mcTrackID(ev->simHitsInfo_) << " phi: " << hit2.phi() << " z: " << hit2.z()); 

	// estimate and set the charge
	seed.setCharge(((hit2_y-hit0_y)*(hit2_x-hit1_x)>(hit2_y-hit1_y)*(hit2_x-hit0_x)?1:-1));

	// save the track and track extra
	evt_seed_tracks.push_back(seed);

	// increment dummy counter for seedID
	seedID++; 
      } // end loop over third layer matches
    } // end loop over first layer matches
  } // end loop over second layer hits
  // still need to do conformal fit + KF fit --> done in fitSeeds()
}
