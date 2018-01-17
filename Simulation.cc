#include <cmath>

#include "TrackerInfo.h"
#include "Simulation.h"
#include "Event.h"
//#define DEBUG
#include "Debug.h"

//#define SOLID_SMEAR
#define SCATTER_XYZ

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk,
                       HitVec& hits, MCHitInfoVec& hitinfos, Event& ev,
                       int itrack, int& charge,
                       const Geometry& geom, TSVec& initTSs)
{
#ifdef DEBUG
  bool debug = true;
#endif
  // This input, 0.5<pt<10 GeV (below ~0.5 GeV does not make 10 layers)
  float pt = Config::minSimPt + g_unif(g_gen)*(Config::maxSimPt - Config::minSimPt);
  pos = SVector3(Config::beamspotX*g_gaus(g_gen), Config::beamspotY*g_gaus(g_gen), Config::beamspotZ*g_gaus(g_gen));

  dprint("production point x=" << pos[0] << " y=" << pos[1] << " z=" << pos[2] << std::endl);

  if (charge==0) {
    if (g_unif(g_gen) > 0.5) charge = -1;
    else charge = 1;
  }

  float phi = 0.5*Config::PI*g_unif(g_gen); // make an angle between 0 and pi/2

  float px = pt * cos(phi);
  float py = pt * sin(phi);

  if (g_unif(g_gen) > 0.5) px*=-1.;
  if (g_unif(g_gen) > 0.5) py*=-1.;

  dprint("phi= " << phi << std::endl);

roll_eta_dice:

  // this generates flat in eta
  float eta = Config::minSimEta + (Config::maxSimEta - Config::minSimEta) * g_unif(g_gen);

  // XXXXMT Hardhack ... exclude transition region eta
  if (Config::TrkInfo.is_transition(eta, 0.1f))
    goto roll_eta_dice;

  float pz  = pt*(1./(std::tan(2*std::atan(std::exp(-eta)))));
  // XXXXMT Commented this out to get ecap_pos only !!!!
  //if (g_unif(g_gen) > 0.5) pz *= -1.;
  dprint("pz="<<pz<<", eta="<<eta);
  mom=SVector3(px,py,pz);
  covtrk=ROOT::Math::SMatrixIdentity();
  //initial covariance can be tricky
  for (int r=0; r<6; ++r) {
    for (int c=0; c<6; ++c) {
      if (r==c) {
        if (r<3) covtrk(r,c)=pow(1.0*pos[r],2);//100% uncertainty on position
        else covtrk(r,c)=pow(1.0*mom[r-3],2);  //100% uncertainty on momentum
      } else {
        covtrk(r,c)=0.;                   //no covariance
      }
    }
  }

  dprint("track with px: " << px << " py: " << py << " pz: " << pz << " pt: " << sqrt(px*px+py*py) << " p: " << sqrt(px*px+py*py+pz*pz) << std::endl);

  TrackState initState;
  initState.parameters=SVector6(pos[0],pos[1],pos[2],mom[0],mom[1],mom[2]);
  initState.errors=covtrk;
  initState.charge=charge;

  TrackState tmpState = initState;

  // useful info for loopers/overlaps

  std::vector<int> layer_counts(Config::nTotalLayers);
  for (int ilayer = 0; ilayer < Config::nTotalLayers; ++ilayer)
  {
    layer_counts[ilayer] = 0;
  }
  int simLayer = -1;

  hits.reserve(Config::nTotHit);
  initTSs.reserve(Config::nTotHit);

  // XXMT4M - This should become while not in outer layer (or at least propState.state == invalid)
  // Try the invalid thingy first ... but would be good to also know what layers are final.

  const PropagationFlags pflags(PF_use_param_b_field);

  for (int ihit = 0; ihit < Config::nTotHit; ++ihit)
  {
    dprintf("\n================================================================================\n");
    dprintf("=== Going for next hit %d from layer %d, xyrzphi=(%f,%f,%f,%f,%f)\n",
            ihit, simLayer, tmpState.x(), tmpState.y(), tmpState.posR(), tmpState.z(), tmpState.posPhi());
    dprintf("--------------------------------------------------------------------------------\n");

    auto propState = propagateHelixToNextSolid(tmpState, geom, pflags);

    float initX   = propState.parameters.At(0);
    float initY   = propState.parameters.At(1);
    float initZ   = propState.parameters.At(2);
    float initPhi = getPhi(initX,initY);
    float initRad = std::sqrt(getRad2(initX,initY));

    UVector3 init_point(initX,initY,initZ);
    simLayer = geom.LayerIndex(init_point);

    // check to see if propagation actually landed inside a layer!
    const auto theInitSolid = geom.InsideWhat(init_point);
    if ( ! theInitSolid ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid." <<std::endl;
      std::cerr << "itrack = " << itrack << ", ihit = " << ihit << ", r = " << initRad << ", r*4cm = " << 4*ihit << ", phi = " << initPhi << ", eta = " << getEta(initRad,initZ) << std::endl;
      std::cerr << "initX = " << initX << ", initY = " << initY << ", initZ = " << initZ << std::endl;
      float outpt  = std::sqrt(getRad2(propState.parameters[3],propState.parameters[4]));
      float outphi = getPhi(propState.parameters[3],propState.parameters[4]);
      std::cerr << "pt = " << outpt << ", pz = " << propState.parameters[5] << ", track phi = " << outphi << ", track eta = " << getEta(pt,pz) << std::endl;

      continue;
    }
    bool is_barrel = theInitSolid->is_barrel_;

scatter_and_smear:

#ifdef SCATTERING
    // PW START
    // multiple scattering. Follow discussion in PDG Ch 32.3
    // this generates a scattering distance yplane and an angle theta_plane in a coordinate
    // system normal to the initial direction of the incident particle ...
    // question: do I need to generate two scattering angles in two planes (theta_plane) and add those or 
    // can I do what I did, generate one (theta_space) and rotate it by another random numbers
    const float z1 = g_gaus(g_gen);
    const float z2 = g_gaus(g_gen);
    const float phismear = g_unif(g_gen)*Config::TwoPI; // random rotation of scattering plane
    // will update for tilt down a few lines
    const float p = sqrt(propState.parameters.At(3)*propState.parameters.At(3)+
                         propState.parameters.At(4)*propState.parameters.At(4)+
                         propState.parameters.At(5)*propState.parameters.At(5));
    UVector3 pvec(propState.parameters.At(3)/p, propState.parameters.At(4)/p, propState.parameters.At(5)/p);

    // now we need to transform to the right coordinate system
    // in this coordinate system z -> z'=z i.e., unchanged (this is a freedom I have since the 
    // definition of the plane defined by nhat is invariant under rotations about nhat
    // x, y change
    // so the transformation is a rotation about z and a translation
    // if the normal vector is given by n == x' = (x1,y1,0) [NB no component along z here]
    // then the 
    // y axis in the new coordinate system is given by y' = z' x x' = (-y1,x1,0)
    UVector3 init_xprime; // normal on surface
    bool init_good = theInitSolid->Normal(init_point, init_xprime);

    if ( ! init_good ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find normal vector at " << init_point <<std::endl;
      break;
    }
    assert(std::abs(init_xprime[2])<1e-10); // in this geometry normal vector is in xy plane

    const float x = Config::xr / std::abs(init_xprime.Dot(pvec)); // take care of dip angle
    const float betac = sqrt(p*p+(.135*.135))/(p*p); // m=130 MeV, pion
    const float theta_0 = 0.0136/(betac*p)*sqrt(x/Config::X0)*(1+0.038*log(x/Config::X0));// eq 32.14
    const float y_plane = z1*x*theta_0/sqrt(12.)+ z2*x*theta_0/2.;
    const float theta_plane = z2*theta_0;
    const float theta_space = sqrt(2)*theta_plane;
    dprint("yplane, theta_space = " << y_plane << ", " << theta_space);

    UVector3 yprime(-init_xprime[1],init_xprime[0],0); // result of cross product with zhat
    //const double phi0 = atan2(xpime[1], init_xprime[0]);

#ifdef SCATTER_XYZ // --> particle will have 3D scatter constrained to detector (including momentum vector)
    const float scatteredX = initX + y_plane *(-init_xprime[1]*cos(phismear)); 
    const float scatteredY = initY + y_plane *(+init_xprime[0]*cos(phismear));
    const float scatteredZ = initZ + y_plane *(           sin(phismear));
    const float scatteredPhi = getPhi(scatteredX,scatteredY);
    const float scatteredRad = std::sqrt(getRad2(scatteredX,scatteredY));
#else // --> particle only has momentum vector scattered 
    const float scatteredX = initX;
    const float scatteredY = initY;
    const float scatteredZ = initZ;
    const float scatteredPhi = getPhi(initX,initY);
    const float scatteredRad = std::sqrt(getRad2(initX,initY));
#endif  // SCATTER_XYZ
    UVector3 pvecprime;
    
    const float v0 = sqrt(2+pow((pvec[1]+pvec[2])/pvec[0],2));
    const float v1 = sqrt(2-pow((pvec[1]-pvec[2]),2));
    const float a = pvec[0]; 
    const float b = pvec[1]; 
    const float c = pvec[2];

    auto sign = [] (const float a) { 
      if ( a > 0 ) 
        return 1;
      else if ( a < 0 ) 
        return -1;
      else
        return 0;
    };
                    
    pvecprime[0] = a*cos(theta_space) + ((b + c)*cos(phismear)*sin(theta_space))/(a*v0) + 
      (a*(b - c)*sin(theta_space)*sin(phismear))/(v1*sign(1 + (b - c)*c));
    pvecprime[1] = b*cos(theta_space) - (cos(phismear)*sin(theta_space))/v0 + 
      ((-1 + pow(b,2) - b*c)*sin(theta_space)*sin(phismear))/(v1*sign(1 + (b - c)*c));
    pvecprime[2] = c*cos(theta_space) - (cos(phismear)*sin(theta_space))/v0 + 
      (std::abs(1 + (b - c)*c)*sin(theta_space)*sin(phismear))/v1; 
    assert(pvecprime.Mag()<=1.0001);

    dprint("theta_space, phismear = " << theta_space << ", " << phismear << std::endl
        << "init_xprime = " << init_xprime << std::endl
        << "yprime = " << yprime << std::endl
        << "phat      = " << pvec << "\t" << pvec.Mag() << std::endl
        << "pvecprime = " << pvecprime << "\t" << pvecprime.Mag() << std::endl
        << "angle     = " << pvecprime.Dot(pvec) << "(" << cos(theta_space) << ")" << std::endl
        << "pt, before and after: " << pvec.Perp()*p << ", "<< pvecprime.Perp()*p);

    pvecprime.Normalize();

    // now update propstate with the new scattered results

    propState.parameters[0] = scatteredX;
    propState.parameters[1] = scatteredY;
    propState.parameters[2] = scatteredZ;

    propState.parameters[3] = pvecprime[0]*p;
    propState.parameters[4] = pvecprime[1]*p;
    propState.parameters[5] = pvecprime[2]*p;
    
    // PW END
    float hitZ    = Config::hitposerrZ*g_gaus(g_gen)+scatteredZ;  // smear in Z after scatter
    float hitPhi  = ((Config::hitposerrXY/scatteredRad)*g_gaus(g_gen))+scatteredPhi; // smear in phi after scatter
    float hitRad  = (Config::hitposerrR)*g_gaus(g_gen)+scatteredRad; // smear in rad after scatter
#else // no Scattering --> use this position for smearing
    // XXMT4K - reuse Config z/r for endcaps in longitudinal/perpendicular sense;
    // XY is still good? Argh ... too late to think.
    // Similar hack below for variance ...
    float poserr_Z = is_barrel ? Config::hitposerrZ : Config::hitposerrR;
    float poserr_R = is_barrel ? Config::hitposerrR : Config::hitposerrZ;

    float hitZ    = poserr_Z * g_gaus(g_gen) + initZ;
    float hitRad  = poserr_R * g_gaus(g_gen) + initRad;
    float hitPhi  = ((Config::hitposerrXY/initRad) * g_gaus(g_gen)) + initPhi;
    // OLD CODE:
    // float hitZ    = Config::hitposerrZ*g_gaus(g_gen)+initZ;
    // float hitPhi  = ((Config::hitposerrXY/initRad)*g_gaus(g_gen))+initPhi;
    // float hitRad  = (Config::hitposerrR)*g_gaus(g_gen)+initRad;
#endif // SCATTERING

#ifdef SOLID_SMEAR
    UVector3 scattered_point(scatteredX,scatteredY,scatteredZ);
    const auto theScatteredSolid = geom.InsideWhat(scattered_point);
    if ( ! theScatteredSolid ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid AFTER scatter." << std::endl;
      std::cerr << "itrack = " << itrack << ", ihit = " << ihit << ", r = " << sqrt(scatteredX*scatteredX + scatteredY*scatteredY) << ", r*4cm = " << 4*ihit << ", phi = " << getPhi(scatteredX,scatteredY) << std::endl;
      std::cerr << "initX = " << initX << ", initY = " << initY << ", initZ = " << initZ << std::endl;
      std::cerr << "scatteredX = " << scatteredX << ", scatteredY = " << scatteredY << ", scatteredZ = " << scatteredZ << std::endl << std::endl;       
      //    continue;
    }

    UVector3 scattered_xprime; // normal on surface
    bool scattered_good = theScatteredSolid->Normal(scattered_point, scattered_xprime);
      
    if ( ! scattered_good ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find normal vector at " << scattered_point <<std::endl;
    }
    assert(std::abs(scattered_xprime[2])<1e-10); // in this geometry normal vector is in xy plane
    
    // smear along scattered yprime (where yprime = z_prime x x_prime, zprime = (0,0,1)

    float xyres_smear = g_gaus(g_gen)*Config::hitposerrXY;

    float hitX = scatteredX - (scattered_xprime[1] * xyres_smear); 
    float hitY = scatteredY + (scattered_xprime[0] * xyres_smear); 
#else // No solid smearing --> use old r-phi smear, whether scattered or not
    float hitRad2 = hitRad*hitRad;

    float hitX    = hitRad*cos(hitPhi);
    float hitY    = hitRad*sin(hitPhi);
    
    float varPhi = Config::varXY/hitRad2;
#endif

#ifdef DEBUG
    if (debug) {
      UVector3 hit_point(hitX,hitY,hitZ);
      const auto theHitSolid = geom.InsideWhat(hit_point);
      if ( ! theHitSolid ) {
        dmutex_guard;
        std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid AFTER scatter+smear." << std::endl;
        std::cerr << "itrack = " << itrack << ", ihit = " << ihit << ", r = " << sqrt(hitX*hitX + hitY*hitY) << ", r*4cm = " << 4*ihit << ", phi = " << getPhi(hitX,hitY) << std::endl;
        std::cerr << "initX = " << initX << ", initY = " << initY << ", initZ = " << initZ << std::endl;
#ifdef SCATTERING
        std::cerr << "scatteredX = " << scatteredX << ", scatteredY = " << scatteredY << ", scatteredZ = " << scatteredZ << std::endl;
#endif
        std::cerr << "hitX = " << hitX << ", hitY = " << hitY << ", hitZ = " << hitZ << std::endl << std::endl; 
      }
    }
#endif

    SMatrixSym33 covXYZ = ROOT::Math::SMatrixIdentity();

#ifdef SOLID_SMEAR
    covXYZ(0,0) = Config::varXY; // yn^2 / (xn^2 + yn^2) * delx^2 + xn^2 / (xn^2 + yn^2) * dely^2
    covXYZ(1,1) = Config::varXY; // xn^2 / (xn^2 + yn^2) * delx^2 + yn^2 / (xn^2 + yn^2) * dely^2
    // covXYZ(0,1) -> -(xn * yn) / (xn^2 + yn^2) * delx^2 + (xn * yn) / (xn^2 + yn^2) * dely^2 
    // covXYZ(1,0)  = covXYZ(0,1)    
    covXYZ(2,2) = Config::varZ;
#else //  SOLID_SMEAR --> covariance for pure cylindrical geometry with smearing
    // XXMT4K - continuation of z/r hack.
    // Still unsure if there is more interplay between varXY and varR ...
    float var_R = is_barrel ? Config::varR : Config::varZ;
    float var_Z = is_barrel ? Config::varZ : Config::varR;

    covXYZ(0,0) = hitX*hitX * var_R / hitRad2 + hitY*hitY*varPhi;
    covXYZ(1,1) = hitX*hitX * varPhi + hitY*hitY * var_R / hitRad2;
    covXYZ(2,2) = var_Z;
    covXYZ(0,1) = hitX*hitY*(var_R / hitRad2 - varPhi);
    covXYZ(1,0) = covXYZ(0,1); // MT: this should be redundant, it's a Sym matrix type.

    // OLD CODE:
    // covXYZ(0,0) = hitX*hitX*Config::varR/hitRad2 + hitY*hitY*varPhi;
    // covXYZ(1,1) = hitX*hitX*varPhi + hitY*hitY*Config::varR/hitRad2;
    // covXYZ(2,2) = Config::varZ;
    // covXYZ(0,1) = hitX*hitY*(Config::varR/hitRad2 - varPhi);
    // covXYZ(1,0) = covXYZ(0,1); // MT: this should be redundant, it's a Sym matrix type.

    dprint("initPhi: " << initPhi << " hitPhi: " << hitPhi << " initRad: " << initRad  << " hitRad: " << hitRad << std::endl
        << "initX: " << initX << " hitX: " << hitX << " initY: " << initY << " hitY: " << hitY << " initZ: " << initZ << " hitZ: " << hitZ << std::endl 
        << "cov(0,0): " << covXYZ(0,0) << " cov(1,1): " << covXYZ(1,1) << " Config::varZ: " << Config::varZ << " cov(2,2): " << covXYZ(2,2) << std::endl 
        << "cov(0,1): " << covXYZ(0,1) << " cov(1,0): " << covXYZ(1,0) << std::endl);
#endif

    UVector3 x1(hitX,hitY,hitZ);

    // XXXXMT4K Are we introducing bias here?
    // 1. I understand smearing needs to change position ... but I do not understand
    // why multiple scattering needs to. Why doesn't it just change momentum?
    // 2. Almost all hits are produced withing epsilon of the entry face. should
    // we try to move it towards the center of the detector?

    if (theInitSolid->Inside(x1) == VUSolid::eOutside)
    {
      dprintf("re-scattering/smearing after hit got placed outside of layer %d\n", simLayer);
      goto scatter_and_smear;
    }

    MCHitInfo hitinfo(itrack, simLayer, layer_counts[simLayer], ev.nextMCHitID());
    hits.emplace_back(x1, covXYZ, hitinfo.mcHitID_);
    hitinfos.emplace_back(hitinfo);

    initTSs.emplace_back(propState); // if no scattering, will just parameters from prop to next layer

    tmpState = propState;

    dprint("hit1Id: " << hitinfo.mcHitID_ <<std::endl
	   << "ihit: " << ihit << " layer: " << simLayer << " counts: " << layer_counts[simLayer]);

    ++layer_counts[simLayer]; // count the number of times passed into layer

    // End track if we hit one of the outer layers.
    if (theInitSolid->is_outer_)
    {
      dprintf("Outer layer reached, track finished.\n");
      break;
    }
  } // end loop over nHitsPerTrack
}
