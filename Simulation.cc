#include <cmath>

#include "Simulation.h"

//#define SOLID_SMEAR
#define SCATTER_XYZ

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, HitVec& hits, unsigned int itrack,
                       int& charge, float pt, Geometry *theGeom, HitVec* initHits)
{
  bool dump = false;
  //assume beam spot width 25um in xy and 5cm in z 
  //  pos=SVector3(0.0,0.0,0.0); // no beamspot smearing
  pos=SVector3(0.0025*g_gaus(g_gen), 0.0025*g_gaus(g_gen), 1.0*g_gaus(g_gen));

  /*
  std::cout << "Simulation Track: " << itrack << std::endl;
  std::cout << "MC vrx: " << pos[0] << " vry: " << pos[1] << std::endl;
  std::cout << "MC IP:  " << sqrt(pos[0]*pos[0] + pos[1]*pos[1]) << std::endl;
  std::cout << "MC Pos Phi (x/y) "  << atan2(pos[0],pos[1]) << " (y/x) " << atan2(pos[1],pos[0]) << std::endl;
  */
  if (dump) std::cout << "pos x=" << pos[0] << " y=" << pos[1] << " z=" << pos[2] << std::endl;

  if (charge==0) {
    if (g_unif(g_gen) > 0.5) charge = -1;
    else charge = 1;
  }

  //float phi = 0.5*TMath::Pi()*(1-g_unif(g_gen)); // make an angle between 0 and pi/2 //fixme
  float phi = 0.5*TMath::Pi()*g_unif(g_gen); // make an angle between 0 and pi/2

  //  std::cout << "MC Gen Phi "  << phi << std::endl;

  float px = pt * cos(phi);
  float py = pt * sin(phi);

  if (g_unif(g_gen)>0.5) px*=-1.;
  if (g_unif(g_gen)>0.5) py*=-1.;

  //  std::cout << "MC Mom Phi (x/y) "  << atan2(px,py) << " (y/x) " << atan2(py,px) << std::endl;
  float pz = pt*(2.3*(g_unif(g_gen)-0.5));//so that we have -1<eta<1

  mom=SVector3(px,py,pz);
  covtrk=ROOT::Math::SMatrixIdentity();
  //initial covariance can be tricky
  for (unsigned int r=0; r<6; ++r) {
    for (unsigned int c=0; c<6; ++c) {
      if (r==c) {
      if (r<3) covtrk(r,c)=pow(1.0*pos[r],2);//100% uncertainty on position
      else covtrk(r,c)=pow(1.0*mom[r-3],2);  //100% uncertainty on momentum
      } else {
        covtrk(r,c)=0.;                   //no covariance
      }
    }
  }

  if (dump) std::cout << "track with px: " << px << " py: " << py << " pz: " << pz << " pt: " << sqrt(px*px+py*py) << " p: " << sqrt(px*px+py*py+pz*pz) << std::endl << std::endl;

  const float hitposerrXY = 0.01;//assume 100mum uncertainty in xy coordinate
  const float hitposerrZ = 0.1;//assume 1mm uncertainty in z coordinate
  const float hitposerrR = hitposerrXY/10.;

  const float varXY  = hitposerrXY*hitposerrXY;
  const float varZ   = hitposerrZ*hitposerrZ;

  TrackState initState;
  initState.parameters=SVector6(pos[0],pos[1],pos[2],mom[0],mom[1],mom[2]);
  initState.errors=covtrk;
  initState.charge=charge;

  TrackState tmpState = initState;

  // useful info for loopers/overlaps

  unsigned int nLayers = theGeom->CountLayers();
  unsigned int layer_counts[nLayers];
  for (unsigned int ilayer=0;ilayer<nLayers;++ilayer){
    layer_counts[ilayer]=0;
  }

  unsigned int nTotHit = nLayers; // can tune this number!
  // to include loopers, and would rather add a break on the code if layer ten exceeded
  // if block BREAK if hit.Layer == theGeom->CountLayers() 
  // else --> if (NMAX TO LOOPER (maybe same as 10?) {break;} else {continue;}
  
  unsigned int simLayer = 0;

  for (unsigned int ihit=0;ihit<nTotHit;++ihit) {  // go to first layer in radius using propagation.h
    //TrackState propState = propagateHelixToR(tmpState,4.*float(ihit+1));//radius of 4*ihit
    auto propState = propagateHelixToNextSolid(tmpState,theGeom);

    float initX   = propState.parameters.At(0);
    float initY   = propState.parameters.At(1);
    float initZ   = propState.parameters.At(2);
    float initPhi = atan2(initY,initX);
    float initRad = sqrt(initX*initX+initY*initY);

    UVector3 init_point(initX,initY,initZ);
    simLayer = theGeom->LayerIndex(init_point);

#ifdef SCATTERING
    // PW START
    // multiple scattering. Follow discussion in PDG Ch 32.3
    // this generates a scattering distance yplane and an angle theta_plane in a coordinate
    // system normal to the initial direction of the incident particle ...
    // question: do I need to generate two scattering angles in two planes (theta_plane) and add those or 
    // can I do what I did, generate one (theta_space) and rotate it by another random numbers
    const float z1 = g_gaus(g_gen);
    const float z2 = g_gaus(g_gen);
    const float phismear = g_unif(g_gen)*TMath::TwoPi(); // random rotation of scattering plane
    const float X0 = 9.370; // cm, from http://pdg.lbl.gov/2014/AtomicNuclearProperties/HTML/silicon_Si.html
    //const float X0 = 0.5612; // cm, for Pb
    float x = 0.0025; //.1 cm  -assumes radial impact. This is bigger than what we have in main
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
    const auto theInitSolid = theGeom->InsideWhat(init_point);
    if ( ! theInitSolid ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid." <<std::endl;
      std::cerr << "itrack = " << itrack << ", ihit = " << ihit << ", r = " << initRad << ", r*4cm = " << 4*ihit << ", phi = " << initPhi << std::endl;
      std::cerr << "initX = " << initX << ", initY = " << initY << ", initZ = " << initZ << std::endl;
      float pt = sqrt(propState.parameters[3]*propState.parameters[3]+
                      propState.parameters[4]*propState.parameters[4]);
      std::cerr << "pt = " << pt << ", pz = " << propState.parameters[5] << std::endl;

      continue;
    }
    UVector3 init_xprime; // normal on surface
    bool init_good = theInitSolid->Normal(init_point, init_xprime);
      
    if ( ! init_good ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find normal vector at " << init_point <<std::endl;
      break;
    }
    assert(std::abs(init_xprime[2])<1e-10); // in this geometry normal vector is in xy plane

    x /= std::abs(init_xprime.Dot(pvec)); // take care of dip angle
    const float betac = sqrt(p*p+(.135*.135))/(p*p); // m=130 MeV, pion
    const float theta_0 = 0.0136/(betac*p)*sqrt(x/X0)*(1+0.038*log(x/X0));// eq 32.14
    const float y_plane = z1*x*theta_0/sqrt(12.)+ z2*x*theta_0/2.;
    const float theta_plane = z2*theta_0;
    const float theta_space = sqrt(2)*theta_plane;
    if ( dump ) 
      std::cout << "yplane, theta_space = " << y_plane << ", " << theta_space << std::endl;

    UVector3 yprime(-init_xprime[1],init_xprime[0],0); // result of cross product with zhat
    //const double phi0 = atan2(xpime[1], init_xprime[0]);

#ifdef SCATTER_XYZ // --> particle will have 3D scatter constrained to detector (including momentum vector)
    const float scatteredX = initX + y_plane *(-init_xprime[1]*cos(phismear)); 
    const float scatteredY = initY + y_plane *(+init_xprime[0]*cos(phismear));
    const float scatteredZ = initZ + y_plane *(           sin(phismear));
    const float scatteredPhi = atan2(scatteredY,scatteredX);
    const float scatteredRad = sqrt(scatteredX*scatteredX+scatteredY*scatteredY);
#else // --> particle only has momentum vector scattered 
    const float scatteredX = initX;
    const float scatteredY = initY;
    const float scatteredZ = initZ;
    const float scatteredPhi = atan2(initY,initX);
    const float scatteredRad = sqrt(initX*initX+initY*initY);
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
    if ( dump ) {
      std::cout << "theta_space, phismear = " << theta_space << ", " << phismear << std::endl;
      std::cout << "init_xprime = " << init_xprime << std::endl;
      std::cout << "yprime = " << yprime << std::endl;
      std::cout << "phat      = " << pvec << "\t" << pvec.Mag() << std::endl;
      std::cout << "pvecprime = " << pvecprime << "\t" << pvecprime.Mag() << std::endl;
      std::cout << "angle     = " << pvecprime.Dot(pvec) << "(" << cos(theta_space) << ")" << std::endl;
      std::cout << "pt, before and after: " << pvec.Perp()*p << ", "<< pvecprime.Perp()*p << std::endl;
    }
    pvecprime.Normalize();

    // now update propstate with the new scattered results

    propState.parameters[0] = scatteredX;
    propState.parameters[1] = scatteredY;
    propState.parameters[2] = scatteredZ;

    propState.parameters[3] = pvecprime[0]*p;
    propState.parameters[4] = pvecprime[1]*p;
    propState.parameters[5] = pvecprime[2]*p;
    
    // PW END
    float hitZ    = hitposerrZ*g_gaus(g_gen)+scatteredZ;  // smear in Z after scatter
    float hitPhi  = ((hitposerrXY/scatteredRad)*g_gaus(g_gen))+scatteredPhi; // smear in phi after scatter
    float hitRad  = (hitposerrR)*g_gaus(g_gen)+scatteredRad; // smear in rad after scatter
#else // no Scattering --> use this position for smearing
    float hitZ    = hitposerrZ*g_gaus(g_gen)+initZ;
    float hitPhi  = ((hitposerrXY/initRad)*g_gaus(g_gen))+initPhi;
    float hitRad  = (hitposerrR)*g_gaus(g_gen)+initRad;
#endif // SCATTERING

#ifdef SOLID_SMEAR
    UVector3 scattered_point(scatteredX,scatteredY,scatteredZ);
    const auto theScatteredSolid = theGeom->InsideWhat(scattered_point);
    if ( ! theScatteredSolid ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid AFTER scatter." << std::endl;
      std::cerr << "itrack = " << itrack << ", ihit = " << ihit << ", r = " << sqrt(scatteredX*scatteredX + scatteredY*scatteredY) << ", r*4cm = " << 4*ihit << ", phi = " << atan2(scatteredY,scatteredX) << std::endl;
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

    float xyres_smear = g_gaus(g_gen)*hitposerrXY;

    float hitX = scatteredX - (scattered_xprime[1] * xyres_smear); 
    float hitY = scatteredY + (scattered_xprime[0] * xyres_smear); 
#else // No solid smearing --> use old r-phi smear, whether scattered or not
    float hitRad2 = hitRad*hitRad;
    float hitX    = hitRad*cos(hitPhi);
    float hitY    = hitRad*sin(hitPhi);

    float varPhi = varXY/hitRad2;
    float varR   = hitposerrR*hitposerrR;
#endif

    if (dump) {
      UVector3 hit_point(hitX,hitY,hitZ);
      const auto theHitSolid = theGeom->InsideWhat(hit_point);
      if ( ! theHitSolid ) {
        std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid AFTER scatter+smear." << std::endl;
        std::cerr << "itrack = " << itrack << ", ihit = " << ihit << ", r = " << sqrt(hitX*hitX + hitY*hitY) << ", r*4cm = " << 4*ihit << ", phi = " << atan2(hitY,hitX) << std::endl;
        std::cerr << "initX = " << initX << ", initY = " << initY << ", initZ = " << initZ << std::endl;
#ifdef SCATTERING
        std::cerr << "scatteredX = " << scatteredX << ", scatteredY = " << scatteredY << ", scatteredZ = " << scatteredZ << std::endl;
#endif
        std::cerr << "hitX = " << hitX << ", hitY = " << hitY << ", hitZ = " << hitZ << std::endl << std::endl; 
      }
    }

    SVector3 x1(hitX,hitY,hitZ);
    SMatrixSym33 covXYZ = ROOT::Math::SMatrixIdentity();

#ifdef SOLID_SMEAR
    covXYZ(0,0) = varXY; // yn^2 / (xn^2 + yn^2) * delx^2 + xn^2 / (xn^2 + yn^2) * dely^2
    covXYZ(1,1) = varXY; // xn^2 / (xn^2 + yn^2) * delx^2 + yn^2 / (xn^2 + yn^2) * dely^2
    // covXYZ(0,1) -> -(xn * yn) / (xn^2 + yn^2) * delx^2 + (xn * yn) / (xn^2 + yn^2) * dely^2 
    // covXYZ(1,0)  = covXYZ(0,1)    
    covXYZ(2,2) = varZ;
#else //  SOLID_SMEAR --> covariance for pure cylindrical geometry with smearing
    covXYZ(0,0) = hitX*hitX*varR/hitRad2 + hitY*hitY*varPhi;
    covXYZ(1,1) = hitX*hitX*varPhi + hitY*hitY*varR/hitRad2;
    covXYZ(2,2) = varZ;
    covXYZ(0,1) = hitX*hitY*(varR/hitRad2 - varPhi);
    covXYZ(1,0) = covXYZ(0,1);

    if (dump) {
      std::cout << "initPhi: " << initPhi << " hitPhi: " << hitPhi << " initRad: " << initRad  << " hitRad: " << hitRad << std::endl
                << "initX: " << initX << " hitX: " << hitX << " initY: " << initY << " hitY: " << hitY << " initZ: " << initZ << " hitZ: " << hitZ << std::endl 
                << "cov(0,0): " << covXYZ(0,0) << " cov(1,1): " << covXYZ(1,1) << " varZ: " << varZ << " cov(2,2): " << covXYZ(2,2) << std::endl 
                << "cov(0,1): " << covXYZ(0,1) << " cov(1,0): " << covXYZ(1,0) << std::endl << std::endl;
    }
#endif

    SVector3 initVecXYZ(initX,initY,initZ);
    Hit initHitXYZ(initVecXYZ,covXYZ,itrack,simLayer,layer_counts[simLayer]); 
    if (nullptr != initHits) {
      initHits->push_back(initHitXYZ);
    }

    Hit hit1(x1,covXYZ,initHitXYZ.mcHitInfo());    // will want to make ihit == ilayer, and last number is number of times in layer -- > set to zero++ for now
    hits.push_back(hit1);
    tmpState = propState;

    if (dump){
      std::cout << "initHitId: " << initHitXYZ.hitID() << " hit1Id: " << hit1.hitID() <<std::endl;
      std::cout << "ihit: " << ihit << " layer: " << simLayer << " counts: " << layer_counts[simLayer] << std::endl;
    }

    ++layer_counts[simLayer]; // count the number of times passed into layer

  } // end loop over nHitsPerTrack
}
