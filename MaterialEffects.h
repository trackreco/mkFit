#ifndef _material_effects_
#define _material_effects_

#include "Config.h"

namespace mkfit {

inline int getZbinME(const float z){return (std::abs(z) * Config::nBinsZME)/(Config::rangeZME);}
inline int getRbinME(const float r){return (         r  * Config::nBinsRME)/(Config::rangeRME);}
inline float getRlVal(const int zb, const int rb){return Config::RlgridME[zb][rb];}
inline float getXiVal(const int zb, const int rb){return Config::XigridME[zb][rb];}

inline int getDetId(const float zin, const float r)
{
  const float z = std::abs(zin);
  //pixel barrel
  if (r<4 && z<28) {
    if (z<20) return 0;  //0.018
    else     return 1;   //0.031
  }
  if (r<8 && z<28) {
    if (z<20) return 2; //0.017
    else     return 3; //0.023
  }
  if (r<12 && z<28) {
    if (z<20) return 4; //0.018
    else     return 5; //0.028
  }
  if (r<17 && z<28) {
    if (z<20) return 6; //0.021
    else     return 7; //0.040
  }

  //pixel endcap
  if (r<17 && z<36) {
    if (r>9.5 && z <32) return 8;  //0.067
    else     return 9;   //0.039
  }
  if (r<17 && z<45) {
    if (r>9.5 && z < 40) return 10; //0.070
    else     return 11; //0.040
  }
  if (r<17 && z >=45) {
    if (r>9.5 && z < 49) return 12; //0.103
    else     return 13; //0.097
  }

  //TIB
  if (r<29 && z<70) {
    if (z<22) return 14; //0.028
    else      return 15; //0.038
  }

  if (r<38 && z<70) {
    if (z<25) return 16; //0.025; used to be z < 27
    else      return 17; //0.034
  }
  if (r<46 && z<70) {
    if (z<44) return 18; //0.037; used to be z < 45
    else      return 19; //0.078
  }
  if (r<55 && z<70) {
    if (z<50) return 20; //0.048
    else      return 21; //0.064
  }

  //TID
  if (r<55 && z>70 && z<120) {
    if (r>35 && z<80) return 22; //0.168
    else if (z<86)    return 23; //0.086 (average of 0.084 and 0.088)
    else if (z<90)    return 24; //0.144
    else if (z<98)    return 25; //0.033
    else if (r>35 && z<104) return 26; //0.157
    else              return 27; //0.078
  }

  //TOB
  if (r<65 && z<120) {
    if (z<17)      return 28; //0.014
    else if (z<70) return 29; //0.032
    else           return 30; //0.052
  }
  if (r<75 && z<120) {
    if (z<17)      return 31; //0.012
    else if (z<70) return 32; //0.026
    else           return 33; //0.038
   }
  if (r<82 && z<120) {
    if (z<17)      return 34; //0.015
    else if (z<70) return 35; //0.035
    else           return 36; //0.061
  }
  if (r<90 && z<120) {
    if (z<17)      return 37; //0.015
    else if (z<70) return 38; //0.035
    else           return 39; //0.043
  }
  if (r<100 && z<120) {
    if (z<17)      return 40; //0.015
    else if (z<70) return 41; //0.036
    else           return 42; //0.033
  }
  if (r<120 && z<120) { //TYPO: this is supposed to be r<120? 
    if (z<17)      return 43; //0.010
    else if (z<70) return 44; //0.021
    else           return 45; //0.022
  }

  //TEC
  if (z>120 && r<120) { //TYPO: changed to r<120 instead of 102
    if (z<128) {
      if (r<35)      return 46; //0.093
      else if (r<55) return 47; //0.084
      else if (r<80) return 48; //0.100    
      else           return 49; //0.194
    }
    if (z<132) {
      if (r<45)      return 50; //0.093
      else if (r<70) return 51; //0.108
      else           return 52; //0.200
    }
    if (z<136) {
      if (r<35)      return 53; //0.093
      else if (r<55) return 54; //0.084
      else if (r<80) return 55; //0.100
      else           return 56; //0.194
    }
    if (z<138) {
      if (r<45)      return 57; //0.093
      else if (r<70) return 58; //0.108
      else           return 59; //0.200
    }
    if (z<142) { //was 143
      if (r<35)      return 60; //0.038
      else if (r<55) return 61; //0.075
      else if (r<80) return 62; //0.038
      else           return 63; //0.075
    }
    if (z<146) {
      if (r<45)      return 64; //0.038
      else           return 65; //0.075
    }
    if (z<150) {
      if (r<35)      return 66; //0.038
      else if (r<55) return 67; //0.075
      else if (r<80) return 68; //0.038
      else           return 69; //0.075
    }
    if (z<153) {
      if (r<45)      return 70; //0.038
      else           return 71; //0.075
    }
    if (z<156) { //was 157
      if (r<35)      return 72; //0.039
      else if (r<55) return 73; //0.078
      else if (r<80) return 74; //0.039
      else           return 75; //0.078
    }
    if (z<160) {
      if (r<45)      return 76; //0.039
      else           return 77; //0.078
    }
    if (z<164) {
      if (r<35)      return 78;
      else if (r<55) return 79;
      else if (r<80) return 80;
      else           return 81;
    }
    if (z<167) {
      if (r<45)      return 82;
      else           return 83;
    }

    if (z<170) {
      if (r<55)      return 84; //0.046
      else if (r<80) return 85; //0.023
      else           return 86; //0.046
    }
    if (z<174) { //was 175
      if (r<45)      return 87;
      else           return 88;
    }

    if (z<178) {
      if (r<55)      return 89;
      else if (r<80) return 90;
      else           return 91;
    }
    if (z<180) {
      if (r<45)      return 92;
      else           return 93;
    }

    if (z<184) { //was 185
      if (r<55)      return 94;
      else if (r<80) return 95;
      else           return 96;
    }
    if (z<188) {
      if (r<45)      return 97;
      else           return 98;
    }

    if (z<192) {
      if (r<55)      return 99;
      else if (r<80) return 100;
      else           return 101;
    }
    if (z<195) {
      if (r<45)      return 102;
      else           return 103;
    }

    if (z<202) {
      if (r<55)      return 104;
      else if (r<80) return 105;
      else           return 106;
    }
    if (z<206) { //was 205
      if (r<45)      return 107;
      else           return 108;
    }

    if (z<210) { //was 209
      if (r<55)      return 109;
      else if (r<80) return 110;
      else           return 111;
    }
    if (z<212) {
      if (r<45)      return 112;
      else           return 113;
    }

    if (z<222) { //was 221
      if (r<55)      return 114;
      else if (r<80) return 115;
      else           return 116;
    }
    if (z<224)       return 117;

    if (z<228) {
      if (r<55)      return 118;
      else if (r<80) return 119;
      else           return 120;
    }
    if (z<232)       return 121;

    if (z<242) {
      if (r<55)      return 122;
      else if (r<80) return 123;
      else           return 124;
    }
    if (z<244)       return 125;

    if (z<248) { //was 249
      if (r<55)      return 126;
      else if (r<80) return 127;
      else           return 128;
    }
    if (z<252)       return 129;

    if (z<264) {  //was 265
      if (r<80)      return 130;
      else           return 131;
    }
    if (z<267)       return 132; //TYPO - this was 252 again

    if (z<270) {
      if (r<80)      return 133;
      else           return 134;
    }
    if (z<280)       return 135;
  }
  return -1;
}

inline void fillZRgridME()
{
  for (int zb = 0; zb < Config::nBinsZME; zb++)
  {
    const float zbf = (zb * Config::rangeZME)/Config::nBinsZME;
    for (int rb = 0; rb < Config::nBinsRME; rb++)
    {
      const float rbf = (rb * Config::rangeRME)/Config::nBinsRME;
      const int detid = getDetId(zbf,rbf);
      Config::RlgridME[zb][rb] = (detid>=0?Config::Rl[detid]:0.f);
      Config::XigridME[zb][rb] = (detid>=0?Config::Xi[detid]:0.f);
    }
  }
}

} // end namespace mkfit
#endif
