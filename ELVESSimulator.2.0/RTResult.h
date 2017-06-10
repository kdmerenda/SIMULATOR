#ifndef _TelescopeSimulatorKG_RTResult_h_
#define _TelescopeSimulatorKG_RTResult_h_

#include <string>

namespace TelescopeSimulatorKG {

  enum RTResult {
    eOK = 0, // nothing happened to photon
    eMissedDiaphragm,
    eMissedPixels,
    eMissedFocalSurface,
    eShadowed, // light hit absorbing surface within the telescope building (camera backside,...)
    eAbsorbed,
    eReflected,
    eReflectedByFilter,
    eReflectedByLens,
    eAbsorbedByFilter,
    eAbsorbedByLens,
    eAbsorbedByMirror,
    eBackscattered,
    eInvalidMercedes,
    eInvalidHitMercedes,
    eInvalidLens,
  };
  
  // the clear text expression of the above enum. MAKE SURE THEY MATCH!
  const std::string RTResultName[] = {"eOK",
				      "MissedDiaphragm",
				      "eMissedPixels",
				      "eMissedFocalSurface",
				      "eShadowed",
				      "eAbsorbed",
				      "eReflected",
				      "eReflectedByFilter",
				      "eReflectedByLens",
				      "eAbsorbedByFilter",
				      "eAbsorbedByLens",
				      "eAbsorbedByMirror",
				      "eBackscattered",
				      "eInvalidMercedes",
				      "eInvalidHitMercedes",
				      "eInvalidLens"};
			 
  
}

#endif
