#include <cmath>

#include "Simulation.h"

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, HitVec& hits,
                       int& charge, float pt, Geometry *theGeom, HitVec* initHits)
{
  bool dump = false;

  unsigned int nTotHit = theGeom->CountLayers();

  //assume beam spot width 1mm in xy and 1cm in z
  pos=SVector3(0.1*g_gaus(g_gen), 0.1*g_gaus(g_gen), 1.0*g_gaus(g_gen));

  if (dump) std::cout << "pos x=" << pos[0] << " y=" << pos[1] << " z=" << pos[2] << std::endl;

  if (charge==0) {
    if (g_unif(g_gen) > 0.5) charge = -1;
    else charge = 1;
  }

  //float phi = 0.5*TMath::Pi()*(1-g_unif(g_gen)); // make an angle between 0 and pi/2 //fixme
  float phi = 0.5*TMath::Pi()*g_unif(g_gen); // make an angle between 0 and pi/2
  float px = pt * cos(phi);
  float py = pt * sin(phi);

  if (g_unif(g_gen)>0.5) px*=-1.;
  if (g_unif(g_gen)>0.5) py*=-1.;
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

  TrackState initState;
  initState.parameters=SVector6(pos[0],pos[1],pos[2],mom[0],mom[1],mom[2]);
  initState.errors=covtrk;
  initState.charge=charge;

  TrackState tmpState = initState;

  // go to first layer in radius using propagation.h
  for (unsigned int nhit=1;nhit<=nTotHit;++nhit) {
    //TrackState propState = propagateHelixToR(tmpState,4.*float(nhit));//radius of 4*nhit
    auto propState = propagateHelixToNextSolid(tmpState,theGeom);

    float initX   = propState.parameters.At(0);
    float initY   = propState.parameters.At(1);
    float initZ   = propState.parameters.At(2);
    float initPhi = atan2(initY,initX);
    float initRad = sqrt(initX*initX+initY*initY);
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
    float x = 0.1; // cm  -assumes radial impact. This is bigger than what we have in main
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
    UVector3 point(initX,initY,initZ);
    const auto theSolid = theGeom->InsideWhat(point);
    if ( ! theSolid ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find solid." <<std::endl;
      std::cerr << "nhit = " << nhit << ", r = " << initRad << ", r*4cm = " << 4*nhit 
                << ", phi = " << initPhi ;
      float pt = sqrt(propState.parameters[3]*propState.parameters[3]+
                      propState.parameters[4]*propState.parameters[4]);
      std::cerr << ", pt = " << pt << ", pz = " << propState.parameters[5] << std::endl;

      continue;
    }
    UVector3 xprime; // normal on surface
    bool good = theSolid->Normal(point, xprime);
      
    if ( ! good ) {
      std::cerr << __FILE__ << ":" << __LINE__ << ": failed to find normal vector." <<std::endl;
      continue;
    }
    assert(std::abs(xprime[2])<1e-10); // in this geometry normal vector is in xy plane

    x /= std::abs(xprime.Dot(pvec)); // take care of dip angle
    const float betac = sqrt(p*p+(.135*.135))/(p*p); // m=130 MeV, pion
    const float theta_0 = 0.0136/(betac*p)*sqrt(x/X0)*(1+0.038*log(x/X0));// eq 32.14
    const float y_plane = z1*x*theta_0/sqrt(12.)+ z2*x*theta_0/2.;
    const float theta_plane = z2*theta_0;
    const float theta_space = sqrt(2)*theta_plane;
    if ( dump ) 
      std::cout << "yplane, theta_space = " << y_plane << ", " << theta_space << std::endl;

    UVector3 yprime(-xprime[1],xprime[0],0); // result of dot product with zhat
    //const double phi0 = atan2(xpime[1], xprime[0]);
    
    const float scatteredX = initX + y_plane *(-xprime[1]*cos(phismear)); 
    const float scatteredY = initY + y_plane *(+xprime[0]*cos(phismear));
    const float scatteredZ = initZ + y_plane *(           sin(phismear));
    const float scatteredPhi = atan2(scatteredY,scatteredX);
    const float scatteredRad = sqrt(scatteredX*scatteredX+scatteredY*scatteredY);

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
      std::cout << "xprime = " << xprime << std::endl;
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
    float hitZ    = hitposerrZ*g_gaus(g_gen)+scatteredZ;
    float hitPhi  = ((hitposerrXY/scatteredRad)*g_gaus(g_gen))+scatteredPhi;
    float hitRad  = (hitposerrR)*g_gaus(g_gen)+scatteredRad;
#else // SCATTERING
    float hitZ    = hitposerrZ*g_gaus(g_gen)+initZ;
    float hitPhi  = ((hitposerrXY/initRad)*g_gaus(g_gen))+initPhi;
    float hitRad  = (hitposerrR)*g_gaus(g_gen)+initRad;
#endif // SCATTERING
    float hitRad2 = hitRad*hitRad;
    float hitX    = hitRad*cos(hitPhi);
    float hitY    = hitRad*sin(hitPhi);

    float varXY  = hitposerrXY*hitposerrXY;
    float varPhi = varXY/hitRad2;
    float varR   = hitposerrR*hitposerrR;
    float varZ   = hitposerrZ*hitposerrZ;

    SVector3 x1(hitX,hitY,hitZ);
    SMatrixSym33 covXYZ = ROOT::Math::SMatrixIdentity();

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

    Hit hit1(x1,covXYZ);    
    hits.push_back(hit1);  
    tmpState = propState;

    SVector3 initVecXYZ(initX,initY,initZ);
    Hit initHitXYZ(initVecXYZ,covXYZ);
    if (nullptr != initHits) {
      initHits->push_back(initHitXYZ);
    }
  }
}
