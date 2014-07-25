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
      if (r==c) {
	if (r<3) covtrk(r,c)=pow(1.0*pos[r],2);//100% uncertainty on position
	else covtrk(r,c)=pow(1.0*mom[r-3],2);  //100% uncertainty on momentum
      } else covtrk(r,c)=0.;                   //no covariance
    }
  }

  //std::cout << "track with p=" << px << " " << py << " " << pz << " pt=" << sqrt(px*px+py*py) << " p=" << sqrt(px*px+py*py+pz*pz) << std::endl;

  // For phi smearing comment out hit poserrXY -- unclear, needed for covariance matrix

  float hitposerrXY = 0.01;//assume 100mum uncertainty in xy coordinate
  float hitposerrZ = 0.1;//assume 1mm uncertainty in z coordinate

  TrackState initState;
  initState.parameters=SVector6(pos[0],pos[1],pos[2],mom[0],mom[1],mom[2]);
  initState.errors=covtrk;
  initState.charge=charge;

  TrackState tmpState = initState;

  //do 4 cm in radius using propagation.h
  for (unsigned int nhit=1;nhit<=nTotHit;++nhit) {
    TrackState propState = propagateHelixToR(tmpState,4.*float(nhit));//radius of 4*nhit

    // Uncomment these lines for phi smearing
    /*
    float initX = propState.parameters.At(0);
    float initY = propState.parameters.At(1);
    float initZ  = propState.parameters.At(2);
    
    float initPhi = convertXYtoPhi(initX,initY);
    float hitPhi = convertXYtoPhiErr(initX,initY,hitposerrXY)*g_gaus(g_gen)+initPhi;
    float hitx = convertPhitoX(hitPhi,initX,initY);
    float hity = convertPhitoY(hitPhi,initX,initY);
    float hitz = hitposerrZ*g_gaus(g_gen)+initZ;
    */
    

    // Comment these lines for phi smearing if using xy smear

    float hitx = hitposerrXY*g_gaus(g_gen)+propState.parameters.At(0);
    float hity = hitposerrXY*g_gaus(g_gen)+propState.parameters.At(1);
    float hitz = hitposerrZ*g_gaus(g_gen)+propState.parameters.At(2);
    //   //float hity = sqrt((pos.At(0) + k*(px*sinAP-py*(1-cosAP)))*(pos.At(0) + k*(px*sinAP-py*(1-cosAP)))+
    //          	(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))*(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))-
    //	   	        hitx*hitx);//try to get the fixed radius


    //std::cout << "hit#" << nhit << " " << hitx << " " << hity << " " << hitz << std::endl;
    //
    SVector3 x1(hitx,hity,hitz);
    
    // Uncomment this line for phi smear cov
    //    SMatrixSym33 covx1 = calcSimCov(hitx,hity,hitz,hitposerrXY,hitposerrZ);
    //    SMatrixSym33 covx2 = calcSimCov(initX,initY,initZ,hitposerrXY,hitposerrZ);

    // Comment out these lines if not using xy smear
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=hitposerrXY*hitposerrXY; 
    covx1(1,1)=hitposerrXY*hitposerrXY;
    covx1(2,2)=hitposerrZ*hitposerrZ;
    
    Hit hit1(x1,covx1);    
    hits.push_back(hit1);  
    tmpState = propState;
    /*
    std::cout << "Cov for New, SmearedXY: " << std::endl;
    std::cout << "0,0: " << covx1(0,0) << " 1,1: " << covx1(1,1) << " 2,2: " << covx1(2,2) << std::endl;
    std::cout << "     " << "        " << " 0,1: " << covx1(0,1) << " 1,0: " << covx1(1,0) << std::endl;
    std::cout << "does this work? error 0,0: " << hit1.error()[0][0] << std::endl << std::endl << std::endl << std::endl;
    
    std::cout << "Cov for New, UnsmearedXY: " << std::endl;
    std::cout << "0,0: " << covx2(0,0) << " 1,1: " << covx2(1,1) << " 2,2: " << covx2(2,2) << std::endl;
    std::cout << "     " << "        " << " 0,1: " << covx2(0,1) << " 1,0: " << covx2(1,0) << std::endl;

    std::cout << "Cov for old method: " << std::endl;
    std::cout << "0,0: " << covx3(0,0) << " 1,1: " << covx3(1,1) << " 2,2: " << covx3(2,2) << std::endl;
    */
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

// Uncomment these lines for phi smearing 

float convertXYtoPhi(float iniX, float iniY){
  //  return TMath::ATan2(iniY,iniX);
  return atan2(iniY,iniX);
}

float convertXYtoPhiErr(float iniX, float iniY, float hitXYerr){
  //  float phiErr = TMath::Sqrt(1.0 / (iniX*iniX + iniY*iniY));
  float phiErr = sqrt(1.0 / (iniX*iniX + iniY*iniY));
  // Assume error on x,y is 100 um, or 0.01 cm
  return phiErr * hitXYerr;
}

float convertPhitoX(float iniPhi, float iniX, float iniY){
  //  return TMath::Sqrt(iniX*iniX + iniY*iniY) * TMath::Cos(iniPhi);
  return sqrt(iniX*iniX + iniY*iniY) * cos(iniPhi);
}

float convertPhitoY(float iniPhi, float iniX, float iniY){
  //  return TMath::Sqrt(iniX*iniX + iniY*iniY) * TMath::Sin(iniPhi);
  return sqrt(iniX*iniX + iniY*iniY) * sin(iniPhi);
}

SMatrixSym33 calcSimCov(float hitX, float hitY, float hitZ, float hitPosErrXY, float hitPosErrZ){
  SMatrixSym33 covXYZ = ROOT::Math::SMatrixIdentity();

  float hitX2 = hitX*hitX;
  float hitY2 = hitY*hitY;
  float rad2  = hitX2 + hitY2;
  float varXY = hitPosErrXY*hitPosErrXY;
  float varZ  = hitPosErrZ*hitPosErrZ;
  float rad3  = rad2*rad2*rad2;

  

 
  covXYZ(0,0) = varXY * (hitX2*hitX2 + hitY2 + hitX2*hitY2) / (rad2*rad2);
  covXYZ(1,1) = varXY * (hitY2*hitY2 + hitX2 + hitX2*hitY2) / (rad2*rad2);
  covXYZ(0,1) = varXY * hitX*hitY * (rad2 - 1) / (rad2*rad2);
  covXYZ(1,0) = covXYZ(0,1);
  covXYZ(2,2) = varZ;
  
  
  return covXYZ;
}
