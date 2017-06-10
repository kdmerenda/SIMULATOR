#ifndef _TelescopeSimulatorKG_RayTracer_h_
#define _TelescopeSimulatorKG_RayTracer_h_

static const char CVSId_TelescopeSimulatorKG_RayTracer[] = 
  "$Id: RayTracer.h 29739 2016-12-20 00:32:07Z rulrich $";

#include "RTResult.h"

#include <string>

namespace fdet {
  class Telescope;
}

namespace utl {
  class RandomEngine;
  class Photon;
  class Point;
}

class TPolyLine3D;
class TObjArray;

namespace TelescopeSimulatorKG {

  //class Photon;
  class Lens;
  class Mirror;
  class Camera;
  class Filter;
  
  /**
     \class RayTracer
     \brief Simulates all ray tracing inside a telescope
     
     \author Ralf Ulrich
     \date Mon Apr 24 23:27:58 ART 2006
  */
  
  class RayTracer {
    
    RayTracer();
    RayTracer& operator=(const RayTracer&);
    RayTracer(const RayTracer&);
    
  public:
    
    RayTracer(const fdet::Telescope& tel,
	      utl::RandomEngine& random,
	      bool doShadow=true, 
	      bool doShadowSupport=true, 
	      bool hasMercedes=true,
	      bool plotPhotonTracks=false,
	      double drawPhotonsProbabilty=0);
    
    ~RayTracer();
    
    RTResult Trace(const utl::Photon& photonIn, utl::Photon& photonOut,
		   int& nreflections, int& col, int& row);
    
    TPolyLine3D* DrawTrack(const utl::Point& p1, const utl::Point& p2, int color) const;

    void SetCorrectorRing(bool f) { fHasCorrectorRing=f; }

    static int GetDebugLevel() { return fgDebugLevel; }
    static void SetDebugLevel(int l) { fgDebugLevel = l;}
    
  private:
    const fdet::Telescope *fTel;
    const utl::RandomEngine *fRandom;
    bool fDoShadow;
    bool fDoShadowSupport;
    bool fHasMercedes;
    bool fHasCorrectorRing;
    
    // for 3d plotting
    bool fPlotPhotonTracks;
    double fDrawPhotonProbability;
    TObjArray* fObjectsPhotons;

    Lens *fLens;
    Camera *fCamera;
    Mirror *fMirror;
    Filter *fFilter;

    static int fgDebugLevel;
  };
  
} // TelescopeSimulatorKG

#endif // _TelescopeSimulatorKG_RayTracer_h_

// Configure (x)emacs for this file ...
// Local Variables:
// mode:c++
// compile-command: "make -C .. -k"
// End:
