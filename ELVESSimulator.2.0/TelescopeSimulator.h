/**
   \file
   Header file of the TelescopeSimulator module

   \author Ralf Ulrich, FZK
   \date Tue Feb 28 12:15:36 CET 2006
*/

#ifndef _TelescopeSimulatorKG_TelescopeSimulator_h_
#define _TelescopeSimulatorKG_TelescopeSimulator_h_

static const char CVSId_TelescopeSimulatorKG_TelescopeSimulator[] =
  "$Id: TelescopeSimulator.h 15752 2010-02-21 05:24:13Z rulrich $";

#include <fwk/VModule.h>

#include <utl/CoordinateSystemPtr.h>
#include <utl/ShadowPtr.h>
#include <utl/Photon.h>

#include <map>
#include <list>
#include <vector>
#include <fstream>


namespace utl {
  class Point;
  class Vector;
  class RandomEngine;
}

namespace fdet {
  class Telescope;
}

namespace evt {
  class Event;
}

namespace TelescopeSimulatorKG {

  /**
     \class TelescopeSimulator

     \brief Simulates the FD telescope

     \author Ralf Ulrich
     \date Fri Apr 21 13:57:20 ART 2006
     \ingroup FDSimModules
  */

  class TelescopeSimulator : public fwk::VModule {

  public:
    typedef std::list<std::pair<utl::Photon, int> > PhotonList;
    typedef PhotonList::iterator PhotonListIterator;
    typedef PhotonList::const_iterator PhotonListConstIterator;

    TelescopeSimulator();

    fwk::VModule::ResultFlag Init();
    fwk::VModule::ResultFlag Run(evt::Event& event);
    fwk::VModule::ResultFlag Finish();

    std::string GetSVNId() const
    { return std::string("$Id: TelescopeSimulator.h 15752 2010-02-21 05:24:13Z rulrich $"); }

    void SetSpotPhotonList(PhotonList& phList);

    /// Returns whether pixel traces are stored for individual light components
    bool StoreLightComponentsAtPixels()
    { return fStoreLightComponentsAtPixels; }

    /// Sets whether pixel traces are stored for individual light components
    void SetStoreLightComponentsAtPixels(const bool store)
    { fStoreLightComponentsAtPixels = store; }

    /// Returns the current verbosity level
    int GetVerbosity() { return fVerbosityLevel; }
    /// Sets the verbosity level
    void SetVerbosity(const int verbosity) { fVerbosityLevel = verbosity; }

  private:
    int fVerbosityLevel;

    // general options
    bool fDoNoShadow;
    bool fDoNoShadowSupport;
    bool fHasNoMercedes;

    // drawing options
    bool fDraw;
    int fDrawNumberOfPhotons;

    // shodow photons output file
    std::string fShadowDataOutName;

    // spot shape mode
    bool fSpotMode;
    PhotonList* fSpotPhotonList;

    bool fStoreLightComponentsAtPixels;

    utl::RandomEngine* fRandomEngine;

    // drum simulations
    bool fDrumMode;
    std::string fGeneralConfigSignature;

    // -- shadow factor calculation --
    std::map<int, unsigned int> fNHit;
    std::map<int, unsigned int> fN_1;
    std::map<int, unsigned int> fN_2;
    std::map<int, unsigned int> fN_3;
    std::map<int, unsigned int> fN_4;
    std::map<int, unsigned int> fN_5;
    std::map<int, unsigned int> fN_6;

    // --------------------------
    // Wed Jul 26 10:51:33 CEST 2006
    // this is for SG spot tables !!!!!!!!!
    static void TransformToLocalCameraCoordinates(const double laz, const double lze, double& caz, double& cze);
    // ------------------------------------------

    REGISTER_MODULE("TelescopeSimulatorKG", TelescopeSimulator);

  };

}


#endif

// Configure (x)emacs for this file ...
// Local Variables:
// mode: c++
// compile-command: "make -C .. -k"
// End:
