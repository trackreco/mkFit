#ifndef _material_effects_
#define _material_effects_

#include "Config.h"

inline int getZbinME(const float z){return (std::abs(z) * Config::nBinsZME)/(Config::rangeZME);}
inline int getRbinME(const float r){return (         r  * Config::nBinsRME)/(Config::rangeRME);}
inline float getRlVal(const int zb, const int rb){return Config::RlgridME[zb][rb];}
inline float getXiVal(const int zb, const int rb){return Config::XigridME[zb][rb];}

inline int getDetId(const float zin, const float r)
{
  const float z = std::abs(zin);
  //pixel barrel
  if (r<6 && z<30) return 0;
  if (r<9 && z<30) {
    if (z<7) return 1;
    else     return 2;
  }
  if (r<20 && z<30) {
    if (z<12) return 3;
    else      return 4;
  }
  //pixel endcap
  if (r<20 && z>30) {
    if      (z<40)  return 5;
    else if (z<280) return 6;
  }
  //TIB
  if (r<30 && z<70) {
    if (z<22) return 7;
    else      return 8;
  }
  if (r<38 && z<70) {
    if (z<27) return 9;
    else      return 10;
  }
  if (r<46 && z<70) {
    if (z<45) return 11;
    else      return 12;
  }
  if (r<55 && z<70) {
    if (z<50) return 13;
    else      return 14;
  }
  //TID
  if (r<55 && z>70 && z<120) {
    if (z<81)       return 15;
    else if (z<94)  return 16;
    else if (z<99)  return 17;
    else if (z<103) return 18;
    else            return 19;
  }
  //TOB
  if (r<65 && z<120) {
    if (z<17)      return 20;
    else if (z<70) return 21;
    else           return 22;
  }
  if (r<75 && z<120) {
    if (z<17)      return 23;
    else if (z<70) return 24;
    else           return 25;
  }
  if (r<82 && z<120) {
    if (z<17)      return 26;
    else if (z<70) return 27;
    else           return 28;
  }
  if (r<90 && z<120) {
    if (z<17)      return 29;
    else if (z<70) return 30;
    else           return 31;
  }
  if (r<100 && z<120) {
    if (z<17)      return 32;
    else if (z<70) return 33; 
    else           return 34; 
  }
  if (r<102 && z<120) {
    if (z<17)      return 35; 
    else if (z<70) return 36; 
    else           return 37; 
  }
  //TEC
  if (z>120 && r<102) {
    if (z<128) {
      if (r<55)      return 38; 
      else if (r<80) return 39;
      else           return 40;
    }
    if (z<132) {
      if (r<45)      return 41;
      else if (r<70) return 42;
      else           return 43;
    }
    if (z<136) {
      if (r<55)      return 44;
      else if (r<80) return 45;
      else           return 46;
    }
    if (z<138) {
      if (r<45)      return 47;
      else if (r<70) return 48;
      else           return 49;
    }

    if (z<143) {
      if (r<35)      return 50;
      else if (r<55) return 51; 
      else if (r<80) return 52; 
      else           return 53; 
    }
    if (z<146) {
      if (r<45)      return 54;
      else           return 55;
    }

    if (z<150) {
      if (r<35)      return 56;
      else if (r<55) return 57;
      else if (r<80) return 58;
      else           return 59;
    }
    if (z<153) {
      if (r<45)      return 60;
      else           return 61;
    }

    if (z<157) {
      if (r<35)      return 62;
      else if (r<55) return 63;
      else if (r<80) return 64;
      else           return 65;
    }
    if (z<160) {
      if (r<45)      return 66;
      else           return 67;
    }

    if (z<164) {
      if (r<35)      return 68;
      else if (r<55) return 69;
      else if (r<80) return 70;
      else           return 71;
    }
    if (z<167) {
      if (r<45)      return 72;
      else           return 73;
    }

    if (z<170) {
      if (r<55)      return 74;
      else if (r<80) return 75;
      else           return 76;
    }
    if (z<175) {
      if (r<45)      return 77;
      else           return 78;
    }

    if (z<178) {
      if (r<55)      return 79;
      else if (r<80) return 80;
      else           return 81;
    }
    if (z<180) {
      if (r<45)      return 82;
      else           return 83;
    }

    if (z<185) {
      if (r<55)      return 84;
      else if (r<80) return 85;
      else           return 86;
    }
    if (z<188) {
      if (r<45)      return 87;
      else           return 88;
    }

    if (z<192) {
      if (r<55)      return 89;
      else if (r<80) return 90;
      else           return 91;
    }
    if (z<195) {
      if (r<45)      return 92;
      else           return 93;
    }

    if (z<202) {
      if (r<55)      return 94;
      else if (r<80) return 95;
      else           return 96;
    }
    if (z<205) {
      if (r<45)      return 97;
      else           return 98;
    }

    if (z<209) {
      if (r<55)      return 99;
      else if (r<80) return 100;
      else           return 101;
    }
    if (z<212) {
      if (r<45)      return 102;
      else           return 103;
    }

    if (z<221) {
      if (r<55)      return 104;
      else if (r<80) return 105;
      else           return 106;
    }
    if (z<224)       return 107;

    if (z<228) {
      if (r<55)      return 108;
      else if (r<80) return 109;
      else           return 110;
    }
    if (z<232)       return 111;

    if (z<242) {
      if (r<55)      return 112;
      else if (r<80) return 113;
      else           return 114;
    }
    if (z<244)       return 115;

    if (z<249) {
      if (r<55)      return 116;
      else if (r<80) return 117;
      else           return 118;
    }
    if (z<252)       return 119;

    if (z<265) {
      if (r<80)      return 120;
      else           return 121;
    }
    if (z<252)       return 122;

    if (z<270) {
      if (r<80)      return 123;
      else           return 124;
    }
    if (z<280)       return 125;
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

#endif
