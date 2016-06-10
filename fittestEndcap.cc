#include "fittestEndcap.h"
#include "Propagation.h"
#include "KalmanUtils.h"

void fittestEndcap() {

  std::cout << "fittestEndcap" << std::endl;

  for (int nt = 0; nt<10; ++nt) {

    float pt = Config::minSimPt+g_unif(g_gen)*(Config::maxSimPt-Config::minSimPt);
    SVector3 pos(Config::beamspotX*g_gaus(g_gen), Config::beamspotY*g_gaus(g_gen), Config::beamspotZ*g_gaus(g_gen));

    int charge = 1;
    if (g_unif(g_gen) > 0.5) charge = -1;

    float phi = 0.5*Config::PI*g_unif(g_gen); // make an angle between 0 and pi/2

    float px = pt * cos(phi);
    float py = pt * sin(phi);

    if (g_unif(g_gen)>0.5) px*=-1.;
    if (g_unif(g_gen)>0.5) py*=-1.;

    //now flat eta between 1.5 and 2.5
    float eta = 1.5+g_unif(g_gen);
    float pz  = pt*(1./(std::tan(2*std::atan(std::exp(-eta)))));
    if (g_unif(g_gen)>0.5) pz*=-1.;

    SVector3 mom=SVector3(1./pt, getPhi(px, py), getTheta(pt, pz));
    SMatrixSym66 covtrk=ROOT::Math::SMatrixIdentity();
    //initial covariance can be tricky
    for (int r=0; r<6; ++r) {
      for (int c=0; c<6; ++c) {
	if (r==c) {
	  if (r<3) covtrk(r,c)=pow(1.0*pos[r],2);//100% uncertainty on position
	  else covtrk(r,c)=pow(1.0*mom[r-3],2);  //100% uncertainty on momentum
	} else {
	  covtrk(r,c)=0.;                      //no covariance
	}
      }
    }

    std::cout << "\nsimulate track with px: " << px << " py: " << py << " pz: " << pz << " pt: " << pt << " eta=" << eta << " p: " << sqrt(px*px+py*py+pz*pz) 
	      << " ipt=" << mom[0] << " phi=" << mom[1] << " theta=" << mom[2] << std::endl;

    TrackState initState;
    initState.parameters=SVector6(pos[0],pos[1],pos[2],mom[0],mom[1],mom[2]);
    initState.errors=covtrk;
    initState.charge=charge;

    TrackState tmpState = initState;
    std::vector<Hit> hits;

    const float hitposerrZ    = 0.001;
    const float hitposerrRPhi = 0.01;
    const float hitposerrR    = 0.1;

    int ndisks = 10;
    float disk_z[ndisks] = {5.,10.,15.,20.,25.,30.,35.,40.,45.,50};
    for (int id = 0; id<ndisks; ++id) {

      tmpState = propagateHelixToZ(tmpState, disk_z[id]);

      SVector3 intersection = tmpState.position();
      // std::cout << std::endl << "intersection x=" << intersection[0] << " y=" << intersection[1] << " z=" << intersection[2] 	      
      // 		<< " r=" << getHypot(intersection.At(0), intersection.At(1)) << " phi=" << getPhi(intersection.At(0), intersection.At(1)) << std::endl;

      float hitZ    = hitposerrZ*g_gaus(g_gen)+intersection.At(2);
      float hitRad  = hitposerrR*g_gaus(g_gen)+getHypot(intersection.At(0), intersection.At(1));
      float hitPhi  = ((hitposerrRPhi/getHypot(intersection.At(0), intersection.At(1)))*g_gaus(g_gen))+getPhi(intersection.At(0), intersection.At(1));

      float hitX    = hitRad*cos(hitPhi);
      float hitY    = hitRad*sin(hitPhi);

      SVector3 hitpos(hitX,hitY,hitZ);

      float varPhi = hitposerrRPhi*hitposerrRPhi/(hitRad*hitRad);
      SMatrixSym33 covXYZ = ROOT::Math::SMatrixIdentity();
      covXYZ(0,0) = hitX*hitX*hitposerrR*hitposerrR/(hitRad*hitRad) + hitY*hitY*varPhi;
      covXYZ(1,1) = hitX*hitX*varPhi + hitY*hitY*hitposerrR*hitposerrR/(hitRad*hitRad);
      covXYZ(2,2) = hitposerrZ*hitposerrZ;
      covXYZ(0,1) = hitX*hitY*(hitposerrR*hitposerrR/(hitRad*hitRad) - varPhi);
      covXYZ(1,0) = covXYZ(0,1);

      Hit hit(hitpos,covXYZ,id);
      hits.push_back(hit);
      // std::cout << "simulate hit with x=" << hitpos[0] << " y=" << hitpos[1] << " z=" << hitpos[2] 
      // 		<< " r=" << hitRad << " phi=" << hitPhi << std::endl;
    }

    tmpState = initState;
    for (int i=0;i<hits.size();++i) {  
      TrackState propState = propagateHelixToZ(tmpState, hits[i].z());
      // std::cout << std::endl << "propagate to hit#" << i << std::endl;
      // std::cout << propState.parameters << std::endl;
      // std::cout << propState.errors << std::endl;

      tmpState = updateParametersEndcap(propState, hits[i].measurementState());
      // std::cout << "update" << std::endl;
      // std::cout << tmpState.parameters << std::endl;
      // std::cout << tmpState.errors << std::endl;

      float chi2 = computeChi2Endcap(propState, hits[i].measurementState());
      // std::cout << "chi2=" << chi2 << std::endl;
    }
  
    std::cout << "found track with parameters : " << tmpState.parameters << std::endl;

  }

  return;
}
