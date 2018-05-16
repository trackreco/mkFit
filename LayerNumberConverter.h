#ifndef LayerNumberConverter_h
#define LayerNumberConverter_h

namespace mkfit {

enum struct TkLayout {phase0 = 0, phase1 = 1};

class LayerNumberConverter {
public:
  LayerNumberConverter(TkLayout layout) : lo_(layout) {}
  unsigned int nLayers() const {
    if (lo_ == TkLayout::phase0) return 69;
    if (lo_ == TkLayout::phase1) return 72;
    return 10;
  }
  int convertLayerNumber(int det, int lay, bool useMatched, int isStereo, bool posZ) const {
    if (det == 1 || det == 3 || det == 5){
      return convertBarrelLayerNumber(det, lay, useMatched, isStereo);
    } else {
      int disk = convertDiskNumber(det, lay, useMatched, isStereo);
      if (disk < 0) return -1;

      int lOffset = 0;
      if (lo_ == TkLayout::phase1) lOffset = 1;
      disk += 17+lOffset;
      if (! posZ) disk += 25+2*lOffset;
      return disk;
    }
    return -1;
  }
  
  int convertBarrelLayerNumber(int cmsswdet, int cmsswlay, bool useMatched, int isStereo) const {
    int lOffset = 0;
    if (lo_ == TkLayout::phase1) lOffset = 1;
    if (cmsswdet==2 || cmsswdet==4 || cmsswdet==6) return -1;//FPIX, TID, TEC
    if (cmsswdet==1) return cmsswlay-1;//BPIX
    if (useMatched) {
      //TIB
      if (cmsswdet==3 && cmsswlay==1 && isStereo==-1) return 3+lOffset;
      if (cmsswdet==3 && cmsswlay==2 && isStereo==-1) return 4+lOffset;
      if (cmsswdet==3 && cmsswlay==3 && isStereo==0 ) return 5+lOffset;
      if (cmsswdet==3 && cmsswlay==4 && isStereo==0 ) return 6+lOffset;
      //TOB
      if (cmsswdet==5 && cmsswlay==1 && isStereo==-1) return 7+lOffset;
      if (cmsswdet==5 && cmsswlay==2 && isStereo==-1) return 8+lOffset;
      if (cmsswdet==5 && cmsswlay==3 && isStereo==0 ) return 9+lOffset;
      if (cmsswdet==5 && cmsswlay==4 && isStereo==0 ) return 10+lOffset;
      if (cmsswdet==5 && cmsswlay==5 && isStereo==0 ) return 11+lOffset;
      if (cmsswdet==5 && cmsswlay==6 && isStereo==0 ) return 12+lOffset;
      return -1;
    } else {
      //TIB
      if (cmsswdet==3 && cmsswlay==1 && isStereo==0) return 3+lOffset;
      if (cmsswdet==3 && cmsswlay==1 && isStereo==1) return 4+lOffset;
      if (cmsswdet==3 && cmsswlay==2 && isStereo==0) return 5+lOffset;
      if (cmsswdet==3 && cmsswlay==2 && isStereo==1) return 6+lOffset;
      if (cmsswdet==3 && cmsswlay==3 && isStereo==0) return 7+lOffset;
      if (cmsswdet==3 && cmsswlay==4 && isStereo==0) return 8+lOffset;
      //TOB
      if (cmsswdet==5 && cmsswlay==1 && isStereo==0) return 9+lOffset;
      if (cmsswdet==5 && cmsswlay==1 && isStereo==1) return 10+lOffset;
      if (cmsswdet==5 && cmsswlay==2 && isStereo==0) return 11+lOffset;
      if (cmsswdet==5 && cmsswlay==2 && isStereo==1) return 12+lOffset;
      if (cmsswdet==5 && cmsswlay==3 && isStereo==0) return 13+lOffset;
      if (cmsswdet==5 && cmsswlay==4 && isStereo==0) return 14+lOffset;
      if (cmsswdet==5 && cmsswlay==5 && isStereo==0) return 15+lOffset;
      if (cmsswdet==5 && cmsswlay==6 && isStereo==0) return 16+lOffset;
      return -1;
    }
  }  
  int convertDiskNumber(int cmsswdet, int cmsswdisk, bool useMatched, int isStereo) const {
    if (cmsswdet==1 || cmsswdet==3 || cmsswdet==5) return -1;//BPIX, TIB, TOB
    if (cmsswdet==2) return cmsswdisk-1;//FPIX
    int lOffset = 0;
    if (lo_ == TkLayout::phase1) lOffset = 1;
    if (useMatched) {
      return -1;
    } else {
      //TID
      if (cmsswdet==4 && cmsswdisk==1 && isStereo==0) return 2+lOffset;
      if (cmsswdet==4 && cmsswdisk==1 && isStereo==1) return 3+lOffset;
      if (cmsswdet==4 && cmsswdisk==2 && isStereo==0) return 4+lOffset;
      if (cmsswdet==4 && cmsswdisk==2 && isStereo==1) return 5+lOffset;
      if (cmsswdet==4 && cmsswdisk==3 && isStereo==0) return 6+lOffset;
      if (cmsswdet==4 && cmsswdisk==3 && isStereo==1) return 7+lOffset;
      //TEC
      if (cmsswdet==6 && cmsswdisk==1 && isStereo==0) return 8+lOffset;
      if (cmsswdet==6 && cmsswdisk==1 && isStereo==1) return 9+lOffset;
      if (cmsswdet==6 && cmsswdisk==2 && isStereo==0) return 10+lOffset;
      if (cmsswdet==6 && cmsswdisk==2 && isStereo==1) return 11+lOffset;
      if (cmsswdet==6 && cmsswdisk==3 && isStereo==0) return 12+lOffset;
      if (cmsswdet==6 && cmsswdisk==3 && isStereo==1) return 13+lOffset;
      if (cmsswdet==6 && cmsswdisk==4 && isStereo==0) return 14+lOffset;
      if (cmsswdet==6 && cmsswdisk==4 && isStereo==1) return 15+lOffset;
      if (cmsswdet==6 && cmsswdisk==5 && isStereo==0) return 16+lOffset;
      if (cmsswdet==6 && cmsswdisk==5 && isStereo==1) return 17+lOffset;
      if (cmsswdet==6 && cmsswdisk==6 && isStereo==0) return 18+lOffset;
      if (cmsswdet==6 && cmsswdisk==6 && isStereo==1) return 19+lOffset;
      if (cmsswdet==6 && cmsswdisk==7 && isStereo==0) return 20+lOffset;
      if (cmsswdet==6 && cmsswdisk==7 && isStereo==1) return 21+lOffset;
      if (cmsswdet==6 && cmsswdisk==8 && isStereo==0) return 22+lOffset;
      if (cmsswdet==6 && cmsswdisk==8 && isStereo==1) return 23+lOffset;
      if (cmsswdet==6 && cmsswdisk==9 && isStereo==0) return 24+lOffset;
      if (cmsswdet==6 && cmsswdisk==9 && isStereo==1) return 25+lOffset;
      return -1;
    }
  }
  TkLayout lo_;
};

} // end namespace mkfit

#endif
