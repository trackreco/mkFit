//
// ********************************************************************
// * This Software is part of the AIDA Unified Solids Library package *
// * See: https://aidasoft.web.cern.ch/USolids                        *
// ********************************************************************
//
// $Id:$
//
// --------------------------------------------------------------------
//
// UTessellatedSolid
//
// Class description:
//
// UTessellatedSolid is a special Geant4 solid defined by a number of
// facets (UVFacet). It is important that the supplied facets shall form a
// fully enclose space which is the solid.
// At the moment only two types of facet can be used for the construction of
// a UTessellatedSolid, i.e. the UTriangularFacet and UQuadrangularFacet.
//
// How to contruct a UTessellatedSolid:
//
// First declare a tessellated solid:
//
//      UTessellatedSolid* solidTarget = new UTessellatedSolid("Solid_name");
//
// Define the facets which form the solid
//
//      double targetSiz = 10*cm ;
//      UTriangularFacet *facet1 = new
//      UTriangularFacet (UVector3(-targetSize,-targetSize,        0.0),
//                         UVector3(+targetSize,-targetSize,        0.0),
//                         UVector3(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      UTriangularFacet *facet2 = new
//      UTriangularFacet (UVector3(+targetSize,-targetSize,        0.0),
//                         UVector3(+targetSize,+targetSize,        0.0),
//                         UVector3(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      UTriangularFacet *facet3 = new
//      UTriangularFacet (UVector3(+targetSize,+targetSize,        0.0),
//                         UVector3(-targetSize,+targetSize,        0.0),
//                         UVector3(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      UTriangularFacet *facet4 = new
//      UTriangularFacet (UVector3(-targetSize,+targetSize,        0.0),
//                         UVector3(-targetSize,-targetSize,        0.0),
//                         UVector3(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      UQuadrangularFacet *facet5 = new
//      UQuadrangularFacet (UVector3(-targetSize,-targetSize,      0.0),
//                           UVector3(-targetSize,+targetSize,      0.0),
//                           UVector3(+targetSize,+targetSize,      0.0),
//                           UVector3(+targetSize,-targetSize,      0.0),
//                           ABSOLUTE);
//
// Then add the facets to the solid:
//
//      solidTarget->AddFacet((UVFacet*) facet1);
//      solidTarget->AddFacet((UVFacet*) facet2);
//      solidTarget->AddFacet((UVFacet*) facet3);
//      solidTarget->AddFacet((UVFacet*) facet4);
//      solidTarget->AddFacet((UVFacet*) facet5);
//
// Finally declare the solid is complete:
//
//      solidTarget->SetSolidClosed(true);
//
// 11.07.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UTessellatedSolid_hh
#define UTessellatedSolid_hh 1

#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "VUSolid.hh"
#include "VUFacet.hh"
#include "UVoxelizer.hh"
#include "UTessellatedSolid.hh"

#ifndef USOLIDSONLY
#include "G4VGraphicsScene.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4VisExtent.hh"

#endif // USOLIDSONLY

struct UVertexInfo
{
  int id;
  double mag2;
};


class UVertexComparator
{
  public:
    bool operator()(const UVertexInfo& l, const UVertexInfo& r)
    {
      return l.mag2 == r.mag2 ? l.id < r.id : l.mag2 < r.mag2;
    }
};

class UTessellatedSolid : public VUSolid
{
  public:  // with description

    UTessellatedSolid();
    virtual ~UTessellatedSolid();

#ifndef USOLIDSONLY
    UTessellatedSolid(const std::string& name);

    UTessellatedSolid(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.
#endif // USOLIDSONLY

    UTessellatedSolid(const UTessellatedSolid& s);
    UTessellatedSolid& operator= (const UTessellatedSolid& s);
    UTessellatedSolid& operator+= (const UTessellatedSolid& right);

    bool AddFacet(VUFacet* aFacet);
    inline VUFacet* GetFacet(int i) const
    {
      return fFacets[i];
    }
    int GetNumberOfFacets() const;

    virtual double GetSurfaceArea();

    virtual VUSolid::EnumInside Inside(const UVector3& p) const;

    virtual bool Normal(const UVector3& p, UVector3& aNormal) const;

    virtual double SafetyFromOutside(const UVector3& p, bool aAccurate = false) const;

    virtual double SafetyFromInside(const UVector3& p, bool aAccurate = false) const;
    virtual UGeometryType GetEntityType() const;

    void SetSolidClosed(const bool t);

    bool GetSolidClosed() const;

    virtual UVector3 GetPointOnSurface() const;

    virtual std::ostream& StreamInfo(std::ostream& os) const;

#ifdef USOLIDSONLY

    UTessellatedSolid(const std::string& name);

    virtual double Capacity()
    {
      return 0;
    }

    virtual double SurfaceArea()
    {
      return GetSurfaceArea();
    }

    inline virtual void GetParametersList(int /*aNumber*/, double* /*aArray*/) const {}
    inline virtual void ComputeBBox(UBBox* /*aBox*/, bool /*aStore = false*/) {}

#endif // USOLIDSONLY

    inline void SetMaxVoxels(int max)
    {
      fVoxels.SetMaxVoxels(max);
    }

    /*
    inline void SetMaxVoxels(const UVector3 &reductionRatio)
    {
      fVoxels.SetMaxVoxels(reductionRatio);
    }
    */

    inline UVoxelizer& GetVoxels()
    {
      return fVoxels;
    }

    virtual VUSolid* Clone() const;

#ifndef USOLIDSONLY

    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const;

    virtual G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v)const;

    virtual G4double DistanceToIn(const G4ThreeVector& p) const;

    virtual G4double DistanceToOut(const G4ThreeVector& p) const;

    virtual G4double DistanceToOut(const G4ThreeVector& p,
                                   const G4ThreeVector& v,
                                   const G4bool calcNorm,
                                   G4bool* validNorm,
                                   G4ThreeVector* norm) const;

    virtual double GetCubicVolume();

//  virtual void ComputeDimensions (UVPVParameterisation* p, const int n, const UVPhysicalVolume* pRep) const;

    virtual bool CalculateExtent(const EAxis pAxis, const UVoxelLimits& pVoxelLimit, const UAffineTransform& pTransform, double& pMin, double& pMax) const;

    // when we would have visualization, these routines would be enabled

    UVector3List* CreateRotatedVertices(const UAffineTransform& pTransform) const;
    // Create the List of transformed vertices in the format required
    // for VUSolid:: ClipCrossSection and ClipBetweenSections.

    // Functions for visualization
    virtual void  DescribeYourselfTo(UVGraphicsScene& scene) const;
//  virtual UVisExtent   GetExtent () const;

#endif // USOLIDSONLY

    double      GetMinXExtent() const;
    double      GetMaxXExtent() const;
    double      GetMinYExtent() const;
    double      GetMaxYExtent() const;
    double      GetMinZExtent() const;
    double      GetMaxZExtent() const;

#ifdef USOLIDSONLY

    virtual inline double DistanceToIn(const UVector3& p, const UVector3& v, double aPstep = UUtils::kInfinity) const
    {
      return DistanceToInCore(p, v, aPstep);
    }

    virtual inline double DistanceToOut(const UVector3& p,
                                        const UVector3& v,
                                        UVector3&       aNormalVector,
                                        bool&           aConvex,
                                        double aPstep = UUtils::kInfinity

                                       ) const
    {
      return DistanceToOutCore(p, v, aNormalVector, aConvex, aPstep);
    }

    void Extent(UVector3& aMin, UVector3& aMax) const;

#endif // USOLIDSONLY

    int AllocatedMemoryWithoutVoxels();

    int AllocatedMemory();

    void DisplayAllocatedMemory();

  private: // without description

    double DistanceToOutNoVoxels(const UVector3& p,
                                 const UVector3& v,
                                 UVector3&       aNormalVector,
                                 bool&           aConvex,
                                 double aPstep = UUtils::kInfinity

                                ) const;

    double DistanceToInCandidates(const std::vector<int>& candidates, const UVector3& aPoint, const UVector3& aDirection /*, double aPstep, const UBits &bits*/) const;

    void DistanceToOutCandidates(const std::vector<int >& candidates, const UVector3& aPoint, const UVector3& direction, double& minDist, UVector3& minNormal, int& minCandidate/*, double aPstep*/ /*,  UBits &bits*/) const;

    double DistanceToInNoVoxels(const UVector3& p, const UVector3& v, double aPstep = UUtils::kInfinity) const;

    void SetExtremeFacets();

    VUSolid::EnumInside InsideNoVoxels(const UVector3& p) const;

    VUSolid::EnumInside InsideVoxels(const UVector3& aPoint) const;

    void Voxelize();

    void CreateVertexList();

    void PrecalculateInsides();

    void SetRandomVectors();

    double DistanceToInCore(const UVector3& p, const UVector3& v, double aPstep = UUtils::kInfinity) const;

    double DistanceToOutCore(const UVector3& p,
                             const UVector3& v,
                             UVector3&       aNormalVector,
                             bool&           aConvex,
                             double aPstep = UUtils::kInfinity

                            ) const;

    int SetAllUsingStack(const std::vector<int>& voxel, const std::vector<int>& max, bool status, UBits& checked);

    void DeleteObjects();
    void CopyObjects(const UTessellatedSolid& s);

    static bool CompareSortedVoxel(const std::pair<int, double>& l, const std::pair<int, double>& r)
    {
      return l.second < r.second;
    }

    double MinDistanceFacet(const UVector3& p, bool simple, VUFacet*& facet) const;


    std::vector<VUFacet*>  fFacets;
    std::set<VUFacet*> fExtremeFacets;  // Does all other facets lie on or behind this surface?
//  std::vector<int> fExtremeFacets2; // Does all other facets lie on or behind this surface?
    UGeometryType           fGeometryType;
    double                 fCubicVolume;
    double                 fSurfaceArea;

    std::vector<UVector3>  fVertexList;

    std::set<UVertexInfo, UVertexComparator> fFacetList;

    UVector3 fMinExtent, fMaxExtent;

    inline bool OutsideOfExtent(const UVector3& p, double tolerance = 0) const
    {
      return (p.x < fMinExtent.x - tolerance || p.x > fMaxExtent.x + tolerance ||
              p.y < fMinExtent.y - tolerance || p.y > fMaxExtent.y + tolerance ||
              p.z < fMinExtent.z - tolerance || p.z > fMaxExtent.z + tolerance);
    }

    bool fSolidClosed;

    static const double dirTolerance;
    std::vector<UVector3> fRandir;

    double fgToleranceHalf;

    int fMaxTries;

    UVoxelizer fVoxels;  // voxelized solid

    UBits fInsides;

    void Initialize();

};

#endif
