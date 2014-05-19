#ifndef _track_
#define _track_

#include "Hit.h"
#include "Matrix.h"
#include <vector>

struct TrackState {
public:
  SVector6 parameters;
  SMatrixSym66 errors;
  int charge;
};

class Track {

public:
  Track(TrackState state,std::vector<Hit> hits, float chi2) {
    state_=state;
    hits_=hits;
    chi2_=chi2;
  }
  Track(int charge,SVector3 position, SVector3 momentum, SMatrixSym66 errors,std::vector<Hit> hits, float chi2) {
    state_.charge=charge;
    state_.errors=errors;
    state_.parameters = SVector6(position.At(0),position.At(1),position.At(2),momentum.At(0),momentum.At(1),momentum.At(2));
    hits_=hits;
    chi2_=chi2;
  }
  Track(int& charge,SVector6& parameters, SMatrixSym66& errors,std::vector<Hit> hits, float chi2) {
    state_.charge=charge;
    state_.errors=errors;
    state_.parameters = parameters;
    hits_=hits;
    chi2_=chi2;
  }

  ~Track(){}

  int& charge() {return state_.charge;}
  SVector3 position() {return SVector3(state_.parameters[0],state_.parameters[1],state_.parameters[2]);}
  SVector3 momentum() {return SVector3(state_.parameters[3],state_.parameters[4],state_.parameters[5]);}
  SVector6& parameters() {return state_.parameters;}
  SMatrixSym66& errors() {return state_.errors;}
  TrackState state() {return state_;}
  std::vector<Hit>& hitsVector() {return hits_;}
  float& chi2() {return chi2_;}
  void resetHits() {hits_.clear();}
  void addHit(Hit hit,float chi2) {hits_.push_back(hit);chi2_+=chi2;}
  unsigned int nHits() {hits_.size();}
  Track clone() {return Track(state_,hits_,chi2_);}

private:
  TrackState state_;
  std::vector<Hit> hits_;
  float chi2_;

};

#endif
