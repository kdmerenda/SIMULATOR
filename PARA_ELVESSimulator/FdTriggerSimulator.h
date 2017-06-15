/**
   \file
   Header file of the FdTriggerSimulator module

   \version $Id: FdTriggerSimulator.h 14717 2009-09-17 20:24:36Z lukas $

   \author Ralf Ulrich
   \author Sergio Petrera
   \author Rossana ...
   \date October 2006
*/

#ifndef _FdTriggerSimulatorOG_FdTriggerSimulator_h_
#define _FdTriggerSimulatorOG_FdTriggerSimulator_h_

static const char CVSId_FdTriggerSimulatorOG_FdTriggerSimulator[] =
  "$Id: FdTriggerSimulator.h 14717 2009-09-17 20:24:36Z lukas $";

#include <utl/config.h>
#include <utl/TimeStamp.h>
#include <fwk/VModule.h>
#include <fevt/SLTData.h>

#include <vector>
#include <map>
#include <list>
#include <set>
#include <string>


#include <utl/Trace-fwd.h>

class TEyeEvent;
class TMirrorEvent;


#ifdef FDLIB_V3R3
#include <stdio.h>
class EyeEventClassifier;
class VTltProcessor;
#endif

namespace utl {
  class TimeInterval;
}

namespace fevt {
  class Telescope;
  class Eye;
  class PixelSimData;
}

namespace fdet {
  class Telescope;
  class Eye;
  class Channel;
}

namespace FdTriggerSimulatorOG {

  /**
     \class FdTriggerSimulator

     \brief Simulates the FLT, SLT, TLT and T3 of the telescopes

     \author Sergio Petrera
     \author Rossana ...
     \author Ralf Ulrich
     \date October 2006
     \ingroup FDSimModules
  */

  class FdTriggerSimulator : public fwk::VModule {

  public:
    enum ETriggerLevel {
      eT3,
      eTLT,
      eSLT,
      eFLT,
      eNoTrigger
    };

    FdTriggerSimulator();
    ~FdTriggerSimulator() { }

    /// Init method of the module
    fwk::VModule::ResultFlag Init();

    /// Run method of the module
    fwk::VModule::ResultFlag Run(evt::Event& event);

    /// Finish method of the module
    fwk::VModule::ResultFlag Finish();

    std::string GetSVNId() const
    { return std::string("$Id: FdTriggerSimulator.h 14717 2009-09-17 20:24:36Z lukas $"); }

  private:
    // top level trigger routines
    bool FLTSim(fevt::Telescope& tel);
    int  SLTSim(fevt::Telescope& tel);
    bool TLTSim(TMirrorEvent *tel, fevt::Telescope& evtTel);
    bool T3Sim(TEyeEvent& eye, fevt::Eye& evtEye);

    /// FLT algorithm
    bool T1Trigger(unsigned int chId, unsigned int pixelthreshold,
                   unsigned int nSamp,
                   const utl::TraceI& trace,
                   fevt::PixelSimData& pixelSim,
                   const fdet::Channel& detChannel);

    /// SLT algorithm
    int T2Trigger(const fdet::Telescope& detTel, unsigned int lastFLTbin);
    int ShiftEventToSLT(fevt::Telescope& tel, int sltTimeOffset);

    /// Prepare list of pixel to be read-out
    void DoPixelList(int sltTimeOffset,
                     const fdet::Telescope& detTel);

    /// Method to search patterns in list
    unsigned int TestPattern(const std::set<int>& pixInCycle, int col, int row,
                             int& multipattern,
                             std::vector<int>& usedInSlt);

    TMirrorEvent * MakeMirrorEvent(fevt::Telescope& tel,
                                   const evt::Event& event,
                                   const utl::TimeInterval& offsetFdSd,
                                   int sltTimeOffset);

    bool AddEyeEvent(fevt::Eye& eye, evt::Event& event,
                     TEyeEvent& EyeEvent);

  private:
    // Datacard variables
    ETriggerLevel  fMinRequiredTriggerLevel;
    ETriggerLevel  fMaxSimTriggerLevel;
    unsigned int   fVerbosity;
    int            fColRO;     //< Number of columns (around FLTed pixel) to read out neighbours
    int            fRowRO;     //< Number of rows (around FLTed pixel) to read out neighbours
    double         fT0;        //< Variable for timing offset

    std::list<int>  fChannelFLT;            //< list of triggered pixels
    std::list<int>  fChannelReadOutList;    //< list of pixel to be read out into the FD-event

    // FLT stuff
    typedef std::map<int,int> FltTrace;
    std::map<int,FltTrace> fFLT;          // Flag of pixel ON at each time
    std::map<int,int>      fMultiplicity; // FLTcount;  /// Count of simultaneous pixel on at each time (multiplicity)

    //  SLT stuff
    std::map<int,fevt::SLTData> fSLT;     // SLT word for each column (1-20) and SLT time slice (1000 ns)

    //TLT stuff
    int fTLTPrintLevel;
    utl::TimeStamp fStartMultiplicityTLT;

    int pageCounter;
    
#ifdef FDLIB_V3R3
  private:
    EyeEventClassifier* fT3;
    VTltProcessor* fTltProcessor;
    FILE* fTltLogFile;
#endif

    REGISTER_MODULE("FdTriggerSimulatorOG", FdTriggerSimulator);

  };



} // FdTriggerSimulatorOG


#endif // _FdTriggerSimulatorOG_FdTriggerSimulator_h_

// Configure (x)emacs for this file ...
// Local Variables:
// mode:c++
// compile-command: "make -C .. -k"
// End:
