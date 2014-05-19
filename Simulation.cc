#include <cmath>

#include "Simulation.h"

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, std::vector<Hit>& hits, int& charge, float pt) {

  unsigned int nTotHit = 10;

  //assume beam spot width 1mm in xy and 1cm in z
  pos=SVector3(0.1*g_gaus(g_gen), 0.1*g_gaus(g_gen), 1.0*g_gaus(g_gen));

  if (charge==0) {
    if (g_unif(g_gen) > 0.5) charge = -1;
    else charge = 1;
  }

  float phi = 0.5*TMath::Pi()*g_unif(g_gen); // make an angle between 0 and pi/2
  float px = pt * cos(phi);
  float py = pt * sin(phi);

  if (g_unif(g_gen)>0.5) px*=-1.;
  if (g_unif(g_gen)>0.5) py*=-1.;
  float pz = pt*(2.3*(g_unif(g_gen)-0.5));//so that we have -1<eta<1

  mom=SVector3(px,py,pz);
  covtrk=ROOT::Math::SMatrixIdentity();
  //initial covariance can be tricky
  for (unsigned int r=0;r<6;++r) {
    for (unsigned int c=0;c<6;++c) {
      if (r==c) covtrk(r,c)=1;
      else covtrk(r,c)=0.5;
    }
  }

  //std::cout << "track with p=" << px << " " << py << " " << pz << " pt=" << sqrt(px*px+py*py) << " p=" << sqrt(px*px+py*py+pz*pz) << std::endl;

  float hitposerrXY = 0.01;//assume 100mum uncertainty in xy coordinate
  float hitposerrZ = 0.1;//assume 1mm uncertainty in z coordinate
  float k=charge*100./(-0.299792458*3.8);
  float curvature = pt*k;
  float ctgTheta=mom.At(2)/pt;

  TrackState initState;
  initState.parameters=SVector6(pos[0],pos[1],pos[2],mom[0],mom[1],mom[2]);
  initState.errors=covtrk;
  initState.charge=charge;

  TrackState tmpState = initState;

  //do 4 cm in radius using propagation.h
  for (unsigned int nhit=1;nhit<=nTotHit;++nhit) {
    TrackState propState = propagateHelixToR(tmpState,4.*float(nhit));//radius of 4*nhit

    float hitx = hitposerrXY*g_gaus(g_gen)+propState.parameters.At(0);
    float hity = hitposerrXY*g_gaus(g_gen)+propState.parameters.At(1);
    //float hity = sqrt((pos.At(0) + k*(px*sinAP-py*(1-cosAP)))*(pos.At(0) + k*(px*sinAP-py*(1-cosAP)))+
    //          	(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))*(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))-
    //	   	        hitx*hitx);//try to get the fixed radius
    float hitz = hitposerrZ*g_gaus(g_gen)+propState.parameters.At(2);

    //std::cout << "hit#" << nhit << " " << hitx << " " << hity << " " << hitz << std::endl;
    SVector3 x1(hitx,hity,hitz);
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=hitposerrXY*hitposerrXY; 
    covx1(1,1)=hitposerrXY*hitposerrXY;
    covx1(2,2)=hitposerrZ*hitposerrZ;
    Hit hit1(x1,covx1);    
    hits.push_back(hit1);  
    tmpState = propState;
  }
  
  /*
  //do 4 cm along path
  for (unsigned int nhit=1;nhit<=nTotHit;++nhit) {
    float distance = 4.*float(nhit);//~4 cm distance along curvature between each hit
    float angPath = distance/curvature;
    float cosAP=cos(angPath);
    float sinAP=sin(angPath);
    float hitx = gRandom->Gaus(0,hitposerr)+(pos.At(0) + k*(px*sinAP-py*(1-cosAP)));
    float hity = gRandom->Gaus(0,hitposerr)+(pos.At(1) + k*(py*sinAP+px*(1-cosAP)));
    //float hity = sqrt((pos.At(0) + k*(px*sinAP-py*(1-cosAP)))*(pos.At(0) + k*(px*sinAP-py*(1-cosAP)))+
    //          	(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))*(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))-
    //	   	        hitx*hitx);//try to get the fixed radius
    float hitz = gRandom->Gaus(0,hitposerr)+(pos.At(2) + distance*ctgTheta);    
    //std::cout << "hit#" << nhit << " " << hitx << " " << hity << " " << hitz << std::endl;
    SVector3 x1(hitx,hity,hitz);
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=hitposerr*hitposerr; 
    covx1(1,1)=hitposerr*hitposerr;
    covx1(2,2)=hitposerr*hitposerr;
    Hit hit1(x1,covx1);    
    hits.push_back(hit1);
  }
  */

}
