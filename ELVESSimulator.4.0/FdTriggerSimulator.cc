/**
   \file
   Implementation of the FdTriggerSimulator module

   \author Sergio Petrera
   \author Rossella Caruso
   \author Ralf Ulrich
   \version $Id: FdTriggerSimulator.cc 29669 2016-11-29 13:48:02Z novotny $
   \date 06 May 2004

*/

//#warning Needs some refurbishment to adapt to the new hardware in HEAT (read comment below)
/* =========================================================
  FdTriggerSimulator needs some refurbishment to adapt to
  the new hardware in HEAT (and AN), i.e.:
           old     |   new
        --------------------
   FADC     1000   | 2000
           (100ns) | (50ns)
   FLT-Mult 1000   | 1000
   SLT      1000   | 1000
========================================================= */


static const char CVSId[] =
  "$Id: FdTriggerSimulator.cc 29669 2016-11-29 13:48:02Z novotny $";

#include <iostream>
#include <iomanip>

#include <utl/config.h>
#include <utl/Reader.h>
#include <utl/ErrorLogger.h>
#include <utl/AugerUnits.h>
#include <utl/MathConstants.h>
#include <utl/PhysicalConstants.h>
#include <utl/TabulatedFunction.h>
#include <utl/TabulatedFunctionErrors.h>
#include <utl/Trace.h>
#include <utl/MultiTabulatedFunction.h>
#include <utl/Point.h>
#include <utl/UTMPoint.h>
#include <utl/TimeStamp.h>
#include <utl/TimeInterval.h>
#include <utl/UTCDateTime.h>
#include <utl/AugerException.h>

#include <fwk/CentralConfig.h>
#include <fwk/SVNGlobalRevision.h>
#include <fwk/LocalCoordinateSystem.h>

#include <det/Detector.h>
#include <fdet/FDetector.h>
#include <fdet/Eye.h>
#include <fdet/Camera.h>
#include <fdet/Channel.h>
#include <fdet/Pixel.h>
#include <fdet/Telescope.h>

#include <evt/Event.h>
#include <evt/ShowerSimData.h>
#include <fevt/FEvent.h>
#include <fevt/Eye.h>
#include <fevt/EyeHeader.h>
#include <fevt/EyeTriggerData.h>
#include <fevt/Telescope.h>
#include <fevt/TelescopeSimData.h>
#include <fevt/TelescopeTriggerData.h>
#include <fevt/SLTData.h>
#include <fevt/PixelSimData.h>
#include <fevt/PixelTriggerData.h>
#include <fevt/Pixel.h>
#include <fevt/ChannelSimData.h>
#include <fevt/Channel.h>

#include <io/FDasToOfflineEventConverter.h>
#include <AugerEvent.h>

//-- FD-Event-Lib
#include <FDEventLibVersion.hh>
#include <MiEvent.hh>
#include <MiEventHeader.hh>
#include <MiRun.hh>
#include <EyeEvent.hh>

#include <FdNumbering.hh>

#ifdef FDLIB_V3R3 // FDEventLib above v3r3
// EYE T3
#include <EyeEventClassifier.hh>
// Telescope TLT
#include <PrototypeTltProcessor.hh>

#ifdef FDLIB_V3R5 // FDEventLib above v3r5
#include <MStcCutTlt.hh> //new multiplicity TLT
#endif // IF FDEventLib above v3r5

#endif // IF FDEventLib above v3r3

#include "SltPatternData.h"
#include "FdTriggerSimulator.h"

#include <sstream>

using namespace FdUtil;
using namespace FdTriggerSimulatorOG;
using namespace std;
using namespace utl;
using namespace evt;
using namespace fwk;
using namespace det;
using namespace fevt;
using namespace io;

// T3 settings from Evb
/** settings for standard telescopes */
#define EVENTCLASSIFIER_PARAMETERS_STD { \
    "Standard parameter set", /* fName */ \
  0.15,  /* fNT3ToNT1RatioCut */ \
   250,  /* fMinMaxCut */ \
   200,  /* fSamePixelTimeCut1 */ \
   200,  /* fSamePixelTimeCut2 */ \
  1000,  /* fCloseShowerCut */ \
   0.5,  /* fSamePixelTimeRatioCut */ \
}

#define EVENTCLASSIFIER_PARAMETERS_OLD { \
    "Old parameter set", /* fName */ \
     0,  /* fNT3ToNT1RatioCut => no kLargeEvents */ \
   200,  /* fMinMaxCut */ \
     0,  /* fSamePixelTimeCut1 => no ratio calculation*/ \
   100,  /* fSamePixelTimeCut2 */ \
  1000,  /* fCloseShowerCut */ \
     1,  /* fSamePixelTimeRatioCut => no ratio cut */ \
}

// settings #1 for HEAT never used

/** settings #2 for HEAT telescopes */
#define EVENTCLASSIFIER_PARAMETERS_HEAT_2 { \
    "Parameter set 2 for HEAT", /* fName */ \
  0.15,  /* fNT3ToNT1RatioCut */ \
   125,  /* fMinMaxCut */ \
   100,  /* fSamePixelTimeCut1 */ \
   100,  /* fSamePixelTimeCut2 */ \
  1000,  /* fCloseShowerCut */ \
   0.5,  /* fSamePixelTimeRatioCut */ \
}


FdTriggerSimulator::FdTriggerSimulator()
#ifdef FDLIB_V3R3
  : fT3(0),
    fTltProcessor(0),
    fTltLogFile(0)
#endif
{
}


VModule::ResultFlag
FdTriggerSimulator::Init()
{
  // :( sorry for changing this
  //counter for pages
  pageCounter=0;
  
  CentralConfig* const cc = CentralConfig::GetInstance();
  Branch topB = cc->GetTopBranch("FdTriggerSimulator");

  topB.GetChild("verbosityLevel").GetData(fVerbosity);
  topB.GetChild("colRO").GetData(fColRO);
  topB.GetChild("rowRO").GetData(fRowRO);

  fStartMultiplicityTLT = UTCDateTime(2030,1,1,0,0,0,0).GetTimeStamp();
  if (topB.GetChild("StartMultiplicityTLT"))
    topB.GetChild("StartMultiplicityTLT").GetData(fStartMultiplicityTLT);

  fTLTPrintLevel = 0;
  if (topB.GetChild("TltPrintLevel"))
    topB.GetChild("TltPrintLevel").GetData(fTLTPrintLevel);

  string levelStr;
  topB.GetChild("minimalTriggerLevel").GetData(levelStr);
  if (levelStr=="T3") {
    fMinRequiredTriggerLevel = eT3;
  } else if (levelStr=="TLT") {
    fMinRequiredTriggerLevel = eTLT;
  } else if (levelStr=="SLT") {
    fMinRequiredTriggerLevel = eSLT;
  } else if (levelStr=="FLT") {
    fMinRequiredTriggerLevel = eFLT;
  } else if (levelStr=="all") {
    fMinRequiredTriggerLevel = eNoTrigger;
  } else {
    ostringstream err;
    err << " unkown enumeration value: \"" << levelStr <<
      "\" used. Check your xml/xsd files!";
    ERROR (err.str());
    return eFailure;
  }

  fMaxSimTriggerLevel = eT3;
  if (!topB.GetChild("maximumSimulatedTriggerLevel")==0) {
    string levelStr;
    topB.GetChild("maximumSimulatedTriggerLevel").GetData(levelStr);
    if (levelStr=="T3") {
      fMaxSimTriggerLevel = eT3;
    } else if (levelStr=="TLT") {
      fMaxSimTriggerLevel = eTLT;
    } else if (levelStr=="SLT") {
      fMaxSimTriggerLevel = eSLT;
    } else if (levelStr=="FLT") {
      fMaxSimTriggerLevel = eFLT;
    } else if (levelStr=="all") {
      fMaxSimTriggerLevel = eNoTrigger;
    } else {
      ostringstream err;
      err << " unkown enumeration value: \"" << levelStr <<
        "\" used. Check your xml/xsd files!";
      ERROR (err.str());
      return eFailure;
    }
  }

  // check FdEventLib version
#ifndef FDLIB_V3R3
  if (fMinRequiredTriggerLevel<eSLT) {
    ERROR("Without FDEventLib >=v3r3 there is no support for FD TLT/T3 simulations! You need to upgrade! ");
    return eFailure;
  }
#endif

  // info output
  ostringstream info;
  info << " Version: "
       << GetVersionInfo(VModule::eRevisionNumber) << "\n"
          " Parameters:\n"
          "            verbosity: " << fVerbosity << "\n"
          " add. columns readout: " << fColRO << "\n"
          "    add. rows readout: " << fRowRO << "\n"
          "multiplicityTLT after: " << fStartMultiplicityTLT << "\n"
          "       TLT printlevel: " << fTLTPrintLevel << "\n"
          "    min trigger level: ";
  switch(fMinRequiredTriggerLevel) {
      case eT3: info << "T3 \n"; break;
      case eTLT: info << "TLT \n"; break;
      case eSLT: info << "SLT \n"; break;
      case eFLT: info << "FLT \n"; break;
      case eNoTrigger: info << "no trigger \n"; break;
      default: info << "UNKOWN \n"; break;
  }
  info << "max sim.trigger level: ";
  switch(fMaxSimTriggerLevel) {
      case eT3: info << "T3 \n"; break;
      case eTLT: info << "TLT \n"; break;
      case eSLT: info << "SLT \n"; break;
      case eFLT: info << "FLT \n"; break;
      case eNoTrigger: info << "no trigger \n"; break;
      default: info << "UNKOWN \n"; break;
  }

  INFO(info);

  fTltLogFile = fopen("Tlt.log", "w");

#ifndef FDLIB_V3R3

    WARNING("\n For TLT and T3 simulations you need FDEventLib >= v3r3 \n");

#endif

  return eSuccess;

}// end of Init



VModule::ResultFlag FdTriggerSimulator::Run(evt::Event& event){
  //page counter increments here
  pageCounter++;
  cout << "TRIGGER FOR PAGE " << pageCounter << endl;
  
  if (!event.HasFEvent()) {
    ERROR("Event has no FEvent.");
    return eFailure;
  }
  FEvent&  fEvent = event.GetFEvent();

  // Generated event print-out ////////////////////////////////////////
  if (fVerbosity >= 1) {
    cout <<"****************************************"
         <<"***************************************"<<endl;
    cout << "***** Event no.  "<<fEvent.GetHeader().GetId()<<endl;;
    cout << "Timestamp " << event.GetHeader().GetTime()<< endl;
    cout <<"****************************************"
         <<"***************************************"<<endl;
  }
  // Generated event print-out ////////////////////////////////////////


  // initialize TLT
#ifdef FDLIB_V3R3
    // INIT TLT
  if (fMaxSimTriggerLevel <= eTLT) {

    ostringstream TLTinfo;
    if (event.GetHeader().GetTime() < fStartMultiplicityTLT) { //prototype TLT
      TLTinfo << "Initializing Prototype-TLT processor ("
              << UTCDateTime(event.GetHeader().GetTime())
              << " is before switch-time " << UTCDateTime(fStartMultiplicityTLT) << ")";
      fTltProcessor = new PrototypeTltProcessor();
    }
    else {
#ifdef FDLIB_V3R5
      TLTinfo << "Initializing Multiplicity-TLT processor ("
              << UTCDateTime(event.GetHeader().GetTime())
              << " is after switch-time " << UTCDateTime(fStartMultiplicityTLT) << ")";
      MiReadout::MStcCutTltParametersRec tlt_parrec = MSTC_TLT_DEFAULT_PARAMETERS;
      fTltProcessor = new MStcCutTlt(&tlt_parrec); //new TLT (GAP-2007-118)
#else

      WARNING("\n For multiplicity TLT simulations you need FDEventLib >= v3r5 (switching to prototype TLT!)\n");
      fTltProcessor = new PrototypeTltProcessor();

#endif


    }
    INFO(TLTinfo);

    fTltProcessor->SetLogfile(fTltLogFile);
    // fTltProcessor->SetPrintLevel(2); //default 2
    fTltProcessor->SetPrintLevel(fTLTPrintLevel);
    //cout << " TLT: load default config" << endl;
    //MiReadout::TltParametersRec tltpar = DEFAULT_TLT_PARAMETERS;
    //fTltProcessor->SetParameters( &tltpar );
    //cout << " TLT: load default config. done." << endl;

    /*
      --- modify some of the TLT parameters ...
      //    15,     fCoincidenceWindow
      //    true,   fEnableMultipleRejection
      //    10,     fNumPatternToAnalyse
      //    10 ,    fKeepRejectedFraction
      //    10,     fFirstSltBin
      //    35,     fLastSltBin
      //    120,    fMaxEventSize
      //    0.40,   fCutParameter1
      //    0.90    fCutParameter2
      */

    ofstream gTltOut("Tlt.out");
    fTltProcessor->PrintSettings(gTltOut);
  }

  if (fMaxSimTriggerLevel <= eT3) {
    // INIT T3
    INFO("Initializing T3 processor");
#if FDEVENTLIB_VERSION_CODE >= ModuleVersionCode (4, 0, 6)
    fT3 = &EyeEventClassifier::GetInstance();
    EyeEventClassifier::GetInstance().SetLogLevel(2);
#else
    fT3 = new EyeEventClassifier(); // note: This leaks now and always has
#endif
    fT3->SetLogLevel(2);
  }

#else

  WARNING("\n For TLT and T3 simulations you need FDEventLib >= v3r3 \n");

#endif

  const fdet::FDetector& theFDet = Detector::GetInstance().GetFDetector();

  for (fevt::FEvent::EyeIterator iEye = fEvent.EyesBegin(ComponentSelector::eInDAQ);
       iEye != fEvent.EyesEnd(ComponentSelector::eInDAQ);
       ++iEye) {

    if (iEye->GetStatus()==fevt::ComponentSelector::eDeSelected)
      continue;

    const fdet::Eye& eyeDet = theFDet.GetEye(*iEye);
    const TimeInterval& offsetFdSd = eyeDet.GetSDTimeOffset();

    bool eyeWithTrigger = false;

    TEyeEvent eyeEvent; // RU Fri May 13 11:22:22 CEST 2005
#if FDEVENTLIB_VERSION_CODE >= ModuleVersionCode(4,0,0)
    TEyeGeometryData* geoData = eyeEvent.GetGeometryData();
#endif

    for (fevt::Eye::TelescopeIterator iTel =  iEye->TelescopesBegin(ComponentSelector::eInDAQ);
         iTel != iEye->TelescopesEnd(ComponentSelector::eInDAQ);
         ++iTel) {

#if FDEVENTLIB_VERSION_CODE >= ModuleVersionCode(4,0,0)
      // Set up the telescope axis pointing for the T3.
      // The angles have to be in degrees in the coordinate system at the
      // telescope position with X pointing east and Z pointing up.
      // Thus the usual LocalCoordinateSystem gobbledigob.
      const fdet::Telescope& telDet = eyeDet.GetTelescope(iTel->GetId());
      const utl::Vector& axis = telDet.GetAxis();
      CoordinateSystemPtr localTelCS = fwk::LocalCoordinateSystem::Create(telDet.GetPosition());
      geoData->SetAxisDirection(iTel->GetId(),
                                axis.GetPhi(localTelCS) / degree,
                                (kPi/2.-axis.GetTheta(localTelCS)) / degree);
#endif

      if (!iTel->HasSimData())
        continue;

      bool foundFLT = FLTSim(*iTel);          // Check for channel FLT trigger

      if (foundFLT || fMinRequiredTriggerLevel>=eNoTrigger) {

        int timeT2_1000 = SLTSim(*iTel);

        if (timeT2_1000 || fMinRequiredTriggerLevel>=eFLT) {

          const int timeShift_100 = ShiftEventToSLT(*iTel, timeT2_1000);
          TMirrorEvent *mirrorEvent = MakeMirrorEvent(*iTel, event, offsetFdSd, timeShift_100);
          bool foundTLT = TLTSim(mirrorEvent, *iTel);

          if (foundTLT || fMinRequiredTriggerLevel>=eSLT) {

            eyeEvent.AddEvent(mirrorEvent);
            eyeWithTrigger = true;
          }

          delete mirrorEvent;

        }
      }
    }// end loop over Telescopes

    if (eyeWithTrigger) {

      if (fVerbosity >= 1)
        cout <<endl<< "**************** Event with SLTrigger on Eye "
             <<(*iEye).GetId()<<" ************************" << endl << endl;

      bool foundT3 = T3Sim(eyeEvent, *iEye);

#ifndef FDLIB_V3R3
      // for the case of no TLT/T3 simulation: put a valid timestamp into the T3 data
      // structure, in order to make the CentralTriggerSimulator work.
      const bool hasSimShower = event.HasSimShower(); // drum+photon simulation has no simshower
      const TimeStamp& coreTime = (hasSimShower ? event.GetSimShower().GetTimeStamp() : event.GetHeader().GetTime());
      const Point& eyePos = Detector::GetInstance().GetFDetector().GetEye(*iEye).GetPosition();
      const Point& corePos = (hasSimShower ? event.GetSimShower().GetPosition() : eyePos);
      const TimeStamp t3Time = coreTime + TimeInterval((eyePos-corePos).GetMag()/kSpeedOfLight);
      iEye->GetTriggerData().SetT3Time(t3Time); // time at core as seen by eye
#endif

      if (foundT3 || fMinRequiredTriggerLevel>=eTLT) {

        INFO("\n\n**************** Event with triggered eye ! ****************\n\n");
        AddEyeEvent(*iEye, event, eyeEvent);

      }
    }
  } // end loop over Eyes


  return eSuccess;

}// end of RunDefault





VModule::ResultFlag
FdTriggerSimulator::Finish()
{
#ifdef FDLIB_V3R3
  if (fMaxSimTriggerLevel<=eTLT && fTltLogFile ) {
    fclose(fTltLogFile);
  }
#endif

  return eSuccess;
} // end of Finish



bool
FdTriggerSimulator::FLTSim(fevt::Telescope& tel)
{
  const fdet::Telescope& detTel = Detector::GetInstance().GetFDetector().GetTelescope(tel);

  TelescopeSimData& telSim = tel.GetSimData();

  const unsigned int telId = tel.GetId();
  const unsigned int eyeId = tel.GetEyeId();

  // RESET FLT
  fChannelReadOutList.clear();
  fChannelFLT.clear();
  fFLT.clear();
  fMultiplicity.clear();

  unsigned int nPixShower = 0;
  unsigned int nPixFLT = 0;
  unsigned int nPixFLTShower = 0;
  unsigned int nPixFLTBg = 0;

  bool haspixel = false;

  // loop channels
  for (unsigned int channelId=1; channelId<=detTel.GetLastPixelId(); ++channelId) {

    const fdet::Channel& detChannel = detTel.GetChannel(channelId);
    unsigned int pixelId = detChannel.GetPixelId();                  // mapped pixel id

    if (!tel.HasChannel(channelId))
        continue;
    fevt::Channel& channel = tel.GetChannel(channelId);

    if (!channel.HasSimData())
        continue;
    fevt::ChannelSimData& channel_sim = channel.GetSimData();

    if (!channel_sim.HasFADCTrace (fevt::FdConstants::eTotal))
        continue;

    int threshold = 0;
    if (!tel.HasPixel(pixelId))
      continue;
    fevt::Pixel& pixel = tel.GetPixel(pixelId);

    if (pixel.GetStatus() == ComponentSelector::eDeSelected) {
      continue;
    }

    if (!pixel.HasSimData())
      continue;
    fevt::PixelSimData& sim = pixel.GetSimData();
    const bool hasShowerPhotons = sim.HasPhotonTrace(fevt::FdConstants::eTotal);

    if (fVerbosity >= 1 && !haspixel) {
      cout << endl << ">============== FIRST LEVEL TRIGGER on Mirror/Eye  "
           << telId << "/" << eyeId << " ===============< "<< endl;
      haspixel = true;
    }


    threshold = sim.GetThreshold();
    int nBox = sim.GetNumSamples();
    const TraceI& fadc_trace = channel_sim.GetFADCTrace(fevt::FdConstants::eTotal);


    if (hasShowerPhotons) {
      nPixShower++;
    }

    // Debug ==========================================================
    // Debug photon signal and ADC's
    // these are only shower photons !!!!
    if (hasShowerPhotons) {
      const utl::TraceD& phTrace = sim.GetPhotonTrace(fevt::FdConstants::eTotal);
      if (fVerbosity >= 2) {
        unsigned int duration = 0, t0_ph = 0;
        for (unsigned int t=0; t<phTrace.GetSize(); ++t) {
          if (phTrace[t] > 0) {
            if (duration == 0) {
              t0_ph = t;
            }
            duration++;
          }
        }
        if (duration>0) {
          cout << endl
               << " photons: pixel ID = " << pixelId
               << " signal duration : " << t0_ph << " -> " << t0_ph+duration-1
               << endl;
          for (unsigned int t=t0_ph; t<t0_ph + duration; ++t) {
            cout << "[" << setw(3)<<t << "] signal = " << setw(10)<<phTrace[t]
                 << " --> ADC = " << setw(4) << fadc_trace[t]
                 << endl;
          }
        }
      }
    }
    // END Debug =====================================================

    // the first possible t1 bin is around 10, IF there is no signal in bin 0 !
    if (T1Trigger(channelId, threshold, nBox, fadc_trace, sim, detChannel)) {

      nPixFLT++;
      if (hasShowerPhotons) {
        nPixFLTShower++;
      } else {
        nPixFLTBg++;
      }
    }

  } // loop pixels

  telSim.SetNumberOfFltPixels(nPixFLT);
  telSim.SetNumberOfFltPixelsFromShower(nPixFLTShower);
  telSim.SetNumberOfFltPixelsFromBackground(nPixFLTBg);
  telSim.SetNumberOfPixelsWithShowerPhotons(nPixShower);

  ostringstream msg;
  msg << "FLT simulation for eye=" << eyeId << " telescope=" << telId;
  msg << " number of FLT pixels=" << nPixFLT;
  INFO(msg);

  return (nPixFLT>0);

}// end of FLTSim



int
FdTriggerSimulator::SLTSim(fevt::Telescope& tel)
{
  if (fMaxSimTriggerLevel>eSLT) {
    if (fVerbosity>1) {
      INFO("skipping SLT simulation");
    }
    return 0;
  }

  const unsigned int telId = tel.GetId();
  const int eyeId = tel.GetEyeId();

  if (fVerbosity >= 1) {
    cout << endl<<">============== SECOND LEVEL TRIGGER on Mirror/Eye  "
         <<telId <<"/"<<eyeId<<" ==============< "<< endl;
  }

//#warning Where to get the length of the simulated traces here?
  // take trace length as photon trace length + full FADC trace length
  // this convention is used in FdElectronicsSimulatorOG
  const fdet::Telescope& detTel = Detector::GetInstance().GetFDetector().GetTelescope(tel);
  const fdet::Camera& detCamera = detTel.GetCamera();
  const unsigned int numFADCBin = detCamera.GetFADCTraceLength();
  const unsigned int lastFLTbin = tel.GetSimData().GetNumberOfPhotonBins() + numFADCBin;


  int timeT2_1000 = T2Trigger(detTel, lastFLTbin);    // T2 time in units of 100 ns
  fevt::TelescopeSimData& telSim = tel.GetSimData();
  //  telSim.SetSltTriggerTime(0);           // convert to ADC bin
  telSim.SetSltTriggerTime(timeT2_1000*10);           // convert to ADC bin

  ostringstream msg;
  msg << "SLT simulation for eye=" << eyeId << " telescope=" << telId;
  if(timeT2_1000) {
    msg << " SLT trigger time=" << timeT2_1000*10 << " [100ns]";
  } else {
    msg << " NO SLT";
  }
  INFO(msg);

  return timeT2_1000;
}



int
FdTriggerSimulator::ShiftEventToSLT(fevt::Telescope& tel, int timeT2_1000)
{
  const fdet::Camera& detCamera = Detector::GetInstance().GetFDetector().GetTelescope(tel).GetCamera();

  //  const int SLTbin = 300;
  //changes here too
  int SLTbin;
  if (pageCounter == 1) {
    SLTbin = 280;
  }else{
    SLTbin = 0;
   }
  
  const int adcTraceLength = detCamera.GetFADCTraceLength();
  cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" <<endl;
  cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" << SLTbin<< " " << adcTraceLength << " "<< timeT2_1000 << endl;
  cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" << endl;

  int timeT2_100 = timeT2_1000*10;

  // if event has no SLT, but we want to trigger anyway:
  // re-define event start bin
  if (fMinRequiredTriggerLevel>=eFLT &&  // -> FLT, all
      timeT2_100<=0) {
    timeT2_100 = SLTbin;
  }


  // if no SLT trigger time is available
  if (timeT2_100<0) {
    return -1;
  }

  // ++++ SLT trigger is only issued every 1000ns: ++++++=
  int timeShift_1000 = int(double(timeT2_100-SLTbin)/10);       // center at SLTbin
  int timeShift_100 = 10 * timeShift_1000;


  // +++++ FIRST check for event start time BEFORE start of trace simulation ++++++
  if (timeT2_100<SLTbin) {
    ostringstream err;
    err << " >>>>> Warning T2 = " << timeT2_100
        << " too low. SLTbin=" << SLTbin << ". Resetting Time Offset (was "
        << timeShift_100 << ")";
    ERROR(err);
    timeShift_100 = 0;
    timeShift_1000 = 0;
  }

//#warning Where to get the length of the simulated traces here?
  // take trace length as photon trace length + full FADC trace length
  // this convention is used in FdElectronicsSimulatorOG
  int numSimT1Bin = tel.GetSimData().GetNumberOfPhotonBins(); // ->ADC bins
  if (numSimT1Bin>0)
    numSimT1Bin += adcTraceLength;
  else
    numSimT1Bin = 2*adcTraceLength; // workaround for missing numberofphotonbins in old Offline files

  // ++++ SECOND check for event end time AFTER end of trace simulation +++++
  if (timeT2_100+adcTraceLength-SLTbin > numSimT1Bin) {
//if (timeT2_100+kActualT1Bins-fBinShift>kNumT1Bins) {
//                1000         300       2000
    ostringstream err;
    err << " >>>>> Warning timeT2_100 = " << timeT2_100
        << " too high: resetting Time Offset (was " << timeShift_100 << ")"
        << " adcTraceLength=" << adcTraceLength
        << " SLTbin="<< SLTbin << " numSimT1Bin=" << numSimT1Bin;
    ERROR(err);

    const double sltBinsPerFltBin = detCamera.GetFADCBinSize() / (100.*ns);
    timeShift_1000 = int( (-adcTraceLength+numSimT1Bin) * sltBinsPerFltBin) / 10;
    timeShift_100 = timeShift_1000 * 10;
    timeT2_1000 = (timeShift_100 + SLTbin) / 10;
    timeT2_100 = timeT2_1000 * 10;
  }

  // Debug ==========================================================
  if (fVerbosity>=1) {
    cout << " >>>>> Shifting traces to SLT trigger time: Time Offset is " << timeShift_100
         << " [100 ns] "
         << " ( = " << timeT2_100 << " - " << SLTbin << " )"
         << " <<<<<" << endl << endl;
  }
  // End Debug ==========================================================

  const fdet::FDetector& detFD = Detector::GetInstance().GetFDetector();
  const fdet::Telescope& detTel = detFD.GetTelescope(tel);

  DoPixelList(timeShift_100, detTel); // Prepare list of pixel to be read-out
  fevt::TelescopeSimData& telSim = tel.GetSimData();
  telSim.SetNumberOfReadOutPixels(fChannelReadOutList.size());
  telSim.SetSltTimeShift(timeShift_100);

  return timeShift_100;

} // end of ShiftEventToSLT




bool
FdTriggerSimulator::TLTSim(TMirrorEvent *tel, fevt::Telescope& evtTel)
{
#ifdef FDLIB_V3R3

  if (fMaxSimTriggerLevel>eTLT) {
    if (fVerbosity>1) {
      INFO("skipping TLT simulation");
    }
    return true;
  }

  const unsigned int telId = tel->GetEventHeader()->GetMirrorNo();
  const unsigned int eyeId = tel->GetEventHeader()->GetEyeNo();

  // Debug ==========================================================
  if (fVerbosity >= 1) {
    ostringstream info;
    info << endl << " >============== THIRD LEVEL TRIGGER on mirror "
         << telId << " eye " << eyeId << " ==============< "<< endl << endl;
    INFO(info);
  }
  // END Debug ==========================================================


  fTltProcessor->Clear();
  fTltProcessor->Decide(tel);

  bool result = fTltProcessor->IsEventAccepted();
  TMirrorEventHeader::EMiEventLabel tltLabel = fTltProcessor->GetEventLabel();

  // set FDAS data structure
  tel->GetEventHeader()->SetEventLabel(tltLabel);

  // set Offline data structure
  fevt::TelescopeTriggerData& telTrigger = evtTel.GetTriggerData();
  telTrigger.SetTLTAccepted(result);
  telTrigger.SetTLTLabel(tel->GetEventHeader()->GetVerboseEventLabel(tltLabel));

  ostringstream thisT3info;
  thisT3info<< " +++++++++++++++++++++ "
            << "TLT "
            << (result ? "accepted" : "rejected")
    //<< " type = " << tel->GetEventHeader()->GetVerboseEventType(miType)
            << " label = " << tel->GetEventHeader()->GetVerboseEventLabel(tltLabel)
            << " ++++++++++++ ";
  INFO(thisT3info);

  return result;


#else // if not using FDLIB_V3R3

  return true;

#endif

} // end of TLTmSim


bool
FdTriggerSimulator::T3Sim(TEyeEvent& eye, fevt::Eye& evtEye)
{
  // Init T3
  evtEye.MakeTriggerData();
  fevt::EyeTriggerData& eyeTrig = evtEye.GetTriggerData();

  eyeTrig.SetT3SDP(0,0,0); // theta and phi of SDP
  eyeTrig.SetT3NPixels(0);


#ifdef FDLIB_V3R3

  if (fMaxSimTriggerLevel>eT3) {
    if (fVerbosity>1) {
      INFO("skipping T3 simulation");
    }
    return true;
  }


  const unsigned int eyeId = eye.GetEventHeader()->GetEyeNo();

  // Debug ==========================================================
  if (fVerbosity >= 1) {
    ostringstream info;
    info << endl << " >============== T3 LEVEL TRIGGER on eye "
         << eyeId << " ==============< "<< endl << endl;
    INFO(info);
  }
  // END Debug ==========================================================


  /*
  TEyeEventHeader::EEventClass event_class = fT3->GetEventClass();

  header->SetEventClass(event_class);


  //  // only EventClass kShowerCandidate & kIsShower are sent as T3 to CDAS !!!
  //  // now also kCloseShower, sa 20040615
  //  if (    event_class == TEyeEventHeader::kShowerCandidate
  //      || event_class == TEyeEventHeader::kIsShower
  //       || event_class == TEyeEventHeader::kCloseShower) {

  */

  // -------- T3 --------
  TEyeEventHeader::EEventClass t3class = eye.GetEventHeader()->GetEventClass();

  ostringstream thisT3info;
  // only do T3 for TLT pre-selected events
  if (t3class!=TEyeEventHeader::kLargeEvent &&
      t3class!=TEyeEventHeader::kRejected &&
      t3class!=TEyeEventHeader::kIsMuon &&
      t3class!=TEyeEventHeader::kNoise) {

    // T3 parameters have to be set according to DAQ
    // if changed in the future, this has to be modified (last change 2016/11/28)
    TimeStamp T3time = TimeStamp(eye.GetEventHeader()->GetTimeStamp()->GetGPSSec());
    EyeEventClassifier::EventClassifierParameters p_STD = EVENTCLASSIFIER_PARAMETERS_STD;
    EyeEventClassifier::EventClassifierParameters p_OLD = EVENTCLASSIFIER_PARAMETERS_OLD;
    EyeEventClassifier::EventClassifierParameters p_HE2 = EVENTCLASSIFIER_PARAMETERS_HEAT_2; 
    if (eyeId < 5) {
      if (T3time < UTCDateTime(2004,5,10).GetTimeStamp()) {
        fT3->SetParameters( p_OLD, true );
      } else {
        fT3->SetParameters( p_STD, true );
      }
    } else {
      if (T3time < UTCDateTime(2012,2,15).GetTimeStamp()) {
	fT3->SetParameters( p_STD, true );
      } else {
	fT3->SetParameters( p_HE2, true );
      }
    }
    fT3->DoClassification (&eye);

    t3class = fT3->GetEventClass();
    eye.GetEventHeader()->SetEventClass(t3class);

    float t3Azimuth = fT3->GetAzimuthAtGround();
    /** Coordinates to send to CDAS.
      */
    int timeAtGround = fT3->GetTimeAtGround();
    /** Number of selected pixels. */
    int nPix = fT3->GetNPixels();
    /** Zenith angle of the SDP normal vector. */
    float sdpTheta =  fT3->GetSDPTheta();
    /** Azimuth angle of the SDP normal vector. */
    float sdpPhi = fT3->GetSDPPhi();
    /*  Total shower signal in ADC counts*/
    float t3TotalSignal =  fT3->GetTotalSignal();
    // T3PixelList* t3PixelList = fT3->GetT3PixelList();


    // fill FADS T3 data
    eye.GetT3Data()->SetNPixels(nPix);
    eye.GetT3Data()->SetSDPTheta(sdpTheta);
    eye.GetT3Data()->SetSDPPhi(sdpPhi);
    eye.GetT3Data()->SetAzimuthAtGround(t3Azimuth);
    eye.GetT3Data()->SetTimeAtGround(timeAtGround);
    eye.GetT3Data()->SetTotalSignal(t3TotalSignal);

    // fill Offline T3 data
    eyeTrig.SetT3Time(TimeStamp(eye.GetEventHeader()->GetTimeStamp()->GetGPSSec(),
                                timeAtGround));
    eyeTrig.SetT3SDP(sdpTheta, sdpPhi, t3Azimuth);
    eyeTrig.SetT3NPixels(nPix);
    eyeTrig.SetT3Class(eye.GetEventHeader()->GetVerboseEventClass(t3class));

    thisT3info << " +++++++++++ T3 resulting eventClass = "
         << eye.GetEventHeader()->GetVerboseEventClass(t3class)
         << " +++++++++++++++";

  } else {

    thisT3info << " +++++++++++ T3 not used due to eventClass = "
         << eye.GetEventHeader()->GetVerboseEventClass(t3class)
         << " +++++++++++++++";

  }

  INFO(thisT3info);

  switch (t3class) {

  case TEyeEventHeader::kCloseShower:
  case TEyeEventHeader::kShowerCandidate:
  case TEyeEventHeader::kHorizontalShower:
  case TEyeEventHeader::kIsShower:
    {
      thisT3info << " ++++++++++++ T3: accepted ";

      // ----------- do the EyeTriggerData ---------------
      // (maybe move this into the FDasToOfflineConverter)
      eyeTrig.SetT3Accepted(true);
    }

    //nT3 ++;
    //if (T3Type==TEyeEventHeader::kUnClassified)
    //T3Type = t3class;

    return true;
    break;

  default:
    thisT3info << " ++++++++++++ T3: rejected ";
    //T3Type = t3class;

    return false;
    break;
  };
  // -------------------

  INFO(thisT3info);

  return true;


#else // if not using FDLIB_V3R3

  return true;

#endif

} // end of T3Sim



/// Prepare list of pixel to be read-out
void
FdTriggerSimulator::DoPixelList(int sltTimeShift,
				const fdet::Telescope& detTel)
{
  vector<int> survivingFLT;    // just for debug output
  vector<int> timeOfStart;     // FLT pixel time shifted to SLT trigger (just for debug)
  vector<int> rejectedFLT;     // just for debug output
  vector<unsigned int> neighbourPixels; // just for debug output

  int nPixON = 0;       // number of FLT pixels in trace time window

  const double fltBinsPer100ns = detTel.GetChannel(1).GetFADCBinSize() / (100.*ns);

  const int sltTimeShiftInFltBins = int(sltTimeShift / fltBinsPer100ns);

  // loop over pixels
  for (list<int>::iterator iChannelFLT = fChannelFLT.begin();
       iChannelFLT != fChannelFLT.end(); ) {

    int channelId = *iChannelFLT;
    const fdet::Channel& detChannel = detTel.GetChannel(channelId);

    // Debug ==========================================================
    if (fVerbosity>=4)
      cout << " DoPixelList test pixel with FLT channelId=" << channelId;
    // END Debug ==========================================================

    // scan shifted trace for FLT
    bool fltInShiftedWindow = false;
    if (fFLT.count(channelId)) {

      // loop FLT trace and search for any FLT in trace time window
      for (FltTrace::const_iterator bin = fFLT[channelId].begin();
           bin != fFLT[channelId].end();
           ++bin) {

        // shift to SLT time
        const int tbin = bin->first;
        const int shiftBin = tbin - sltTimeShiftInFltBins;

        // check trace time window
        if (shiftBin>=0 && shiftBin<detChannel.GetFADCTraceLength()) {
          timeOfStart.push_back(tbin);
          survivingFLT.push_back(channelId);
          fltInShiftedWindow = true;
          break;
        }
      } // end loop bins

      if (fltInShiftedWindow) {

        // Debug ==========================================================
        if (fVerbosity>=4)
          cout << " : readout " << endl;
        // END Debug ==========================================================

        nPixON++;
        fChannelReadOutList.push_back(channelId);
        iChannelFLT++; // move to next channel

      } else {

        // Debug ==========================================================
        if(fVerbosity >=4)
          cout << " : has no FLT in trace (REMOVE FLT)!" << endl;
        // END Debug ==========================================================

        rejectedFLT.push_back(channelId);
        iChannelFLT = fChannelFLT.erase(iChannelFLT);       // turn off FLT

        // Debug ==========================================================
        if(fVerbosity>=2) // by sergio
          cout << " &&&& Pixel " << channelId << " is out of time (REMOVE FLT)!" << endl;
        // END Debug ==========================================================

      } // if else (fltInShiftedWindow)

    } // if (fFLT.count(channelId))
    else {

      iChannelFLT++; // move to next channel
      // Debug ==========================================================
      if (fVerbosity>=4)
        cout << " : ERROR fFLT not filled !" << endl;
      // END Debug ==========================================================

    }

  } // end loop channels



  // select pixels close to the FLT pixels (NEIGHBORS)
  for (list<int>::const_iterator iChannelFLT = fChannelFLT.begin();
       iChannelFLT != fChannelFLT.end();
       ++iChannelFLT) {

    int channelId = (*iChannelFLT);
    //int pixelId = detTel.GetChannel(channelId).GetPixelId();

    int col = ((channelId-1)/22) + 1;
    int row = ((channelId-1)%22) + 1;

    for (int c_neigh = -fColRO; c_neigh <= fColRO; ++c_neigh) {
      for (int r_neigh = -fRowRO; r_neigh <= fRowRO; ++r_neigh) {

        int c = col + c_neigh;
        int r = row + r_neigh;
        if (c>=1 && c<=20 &&
            r>=1 && r<=22) {

          int neigh_channel = 22 * (c - 1) + r;

          // flag non-FLT pixel in vicinity of FLT
          if (find(fChannelReadOutList.begin(), fChannelReadOutList.end(), neigh_channel)
              == fChannelReadOutList.end()) {
            fChannelReadOutList.push_back(neigh_channel);
            neighbourPixels.push_back(neigh_channel);
          }
        }
      } // loop rows
    } // loop cols
  } // end loop channels



  // now add the virtual channel
  for (list<int>::const_iterator iReadOutChannel = fChannelReadOutList.begin();
       iReadOutChannel != fChannelReadOutList.end();
       ++iReadOutChannel) {

    unsigned int channelId = (*iReadOutChannel);

    // check if already virtual channel and not pixel
    if (channelId>detTel.GetLastPixelId())
      continue;

    int virtual_channel = detTel.GetChannel(channelId).GetVirtualChannelId();

    if (find(fChannelReadOutList.begin(), fChannelReadOutList.end(), virtual_channel)
        == fChannelReadOutList.end()) {
      fChannelReadOutList.push_back(virtual_channel);
    }

  }
  // -- end algorithm --



  // DEBUG ==========================================================================
  // T1 summary
  if (fVerbosity>=1) {

    int tstartmin = 9999999;
    int tstartmax =-9999999;

    int maxmult = 0;
    int t_maxmult = 0;
    bool first = true;

    for (unsigned int k=0; k<timeOfStart.size(); ++k) {

      int t_bin = timeOfStart[k];

      int mult = 0;
      if (fMultiplicity.count(t_bin)) {
        mult = fMultiplicity[t_bin];
      }

      if (fVerbosity >=4) {
        cout << " surviving FLT in pixelId=" << setw(2) << survivingFLT[k]
             << " time of start " << setw(4) << timeOfStart[k]
             << " FLT multiplicity at bin " << setw(3) << mult
             << endl;
      } else if (fVerbosity >=2) { // by sergio
        cout << "flt pixel no " << setw(2) << k
             << " time of start " << setw(4) << timeOfStart[k]
             << " multiplicity " << setw(3) << mult
             << endl;
      }


      if (first) {
        first = false;
        tstartmin = timeOfStart[k];
        tstartmax = timeOfStart[k];
        maxmult = mult;
        t_maxmult = timeOfStart[k];
      } else {
        if (timeOfStart[k]<tstartmin) tstartmin = timeOfStart[k];
        if (timeOfStart[k]>tstartmax) tstartmax = timeOfStart[k];
        if (mult>maxmult) {
          maxmult = mult;
          t_maxmult = timeOfStart[k];
        }
      }
    }

    if (fVerbosity >=4) { // by Ralf
      for (unsigned int k=0; k<rejectedFLT.size(); ++k) {
        int channelId = rejectedFLT[k];
        cout << " rejected FLT in channelId=" << setw(2) << channelId;
        if (fVerbosity >=5) {
          cout << " FLTtrace: ";
          if (fFLT.count(channelId)) {
            bool firstInTrace = true;
            int previous = -100;
            for (FltTrace::const_iterator bin = fFLT[channelId].begin();
                 bin != fFLT[channelId].end();
                 ++bin) {
              if (previous!=(bin->first-1)) {
                if (firstInTrace) {
                  firstInTrace = false;
                  cout << bin->first;
                } else {
                  cout << "-" << previous << ", " << bin->first;
                }
              }
              previous = bin->first;
            } // loop tbin
            if (!firstInTrace) {
              cout << "-end";
            }
          } // if channel had FLT
        }
        cout << endl;

      }
    }

    cout << endl
         << " NpixON=" << nPixON
         << ",  NpixRO=" << fChannelReadOutList.size()
         << ", Max mult= " << maxmult
         << " @ t = " << t_maxmult
         << ", time of start: min " << tstartmin
         << ",  max " << tstartmax
         << endl;

    if (fVerbosity>=2) {
      cout << endl << " *** List of Read Out pixels ***" << endl;

      int n = 0;
      for (list<int>::const_iterator iReadOutPixel = fChannelReadOutList.begin();
           iReadOutPixel!=fChannelReadOutList.end();
           ++iReadOutPixel) {                      // loop read out pixels pixels

        unsigned int p = int(*iReadOutPixel);

        if (p <= detTel.GetLastPixelId()) {

          int col = ((p-1) / 22) + 1;
          int row = ((p-1) % 22) + 1;
          n++;

          if (find(neighbourPixels.begin(), neighbourPixels.end(), p)==neighbourPixels.end()) {

            cout << setw(4) << n
                 << " pixel " << setw(3) << p
                 << " col " << setw(3) << col
                 << " row " << setw(3) << row
                 << endl;

          } else {

            cout << setw(4) << n
                 << " pixel " << setw(3) << p
                 << " col " << setw(3) << col
                 << " row " << setw(3) << row
                 << " --> neighbour "
                 << endl;
          }

        } else {
          int virt = (p-1)%2 + 1;
          int col = ((p-1) - detTel.GetLastPixelId())/2 + 1;
          cout << setw(4) << n
               << " pixel " << setw(3) << p
               << " col " << setw(3) << col
               << " virt " << setw(3) << virt
               << " --> virtual"
               << endl;
        }
      }
    }
  }
  // ENd DEBUG ========================================================================

} // End of DoPixelList


bool
FdTriggerSimulator::T1Trigger(unsigned int channelId,
			      unsigned int pixelthresh,
			      unsigned int nSamp,
			      const TraceI& trace,
			      PixelSimData& pixelSim,
			      const fdet::Channel& detChannel)
{
  const unsigned int FLTprolongation = detChannel.GetFLTProlongation(); // [bins]

  unsigned int  tstart = 0;
  unsigned int  tend = 0;
  unsigned int  prolongationCounter = 0;
  bool pulse = false;
  bool hasFLT = false;

  unsigned int maxBoxcarSum  = 0;
  double meanBoxcarSum = 0;
  double rmsBoxcarSum  = 0;
  unsigned int ntslots = 0;
  double squareBoxcarSum = 0;

  vector<int> boxcar(nSamp, 0);

  for (unsigned int iBin=0; iBin<trace.GetSize(); ++iBin) { // loop over tslots

    // build boxcar
    for (unsigned int i=0; i<nSamp-1; ++i)
      boxcar[i] = boxcar[i+1];
    boxcar[nSamp-1] = trace [iBin];

    // sum up
    unsigned int boxcarsum = 0;
    for(unsigned int l=0; l<nSamp; ++l) {
      boxcarsum += boxcar[l];
    }

    // neglect the first partially filled boxcarsums
    if (iBin>nSamp-1) {
      ntslots += 1;
      meanBoxcarSum += boxcarsum;
      squareBoxcarSum += boxcarsum*boxcarsum;
    }

    if (boxcarsum>maxBoxcarSum)
      maxBoxcarSum = boxcarsum;

    // Test FLT on running sums wrt threshold
    bool bin_above_threshold = (boxcarsum > pixelthresh);
    if (bin_above_threshold) {

      hasFLT = true;
      prolongationCounter = 0;      // FLT prolongation counter

      if (!pulse) {                 // start of new FLT region (pulse)

        pulse = true;
        tstart = iBin;              // start of pulse/region

        // DEBUG ==========================================================================
        if (fVerbosity>=1){
          cout << "##### FLT start " << setw(5) << iBin
               << " pixel " << setw(3) << channelId
               << "; sum(" << nSamp << ") = " << setw(5) << boxcarsum
               << " threshold = " << pixelthresh
               << endl;
        }
        // ENDDEBUG ==========================================================================

      }

      // DEBUG ==========================================================================
      if (fVerbosity>=2) {
        int col = ((channelId - 1) / 22) + 1;
        int row = ((channelId - 1) % 22) + 1;
        cout << " > FLT Pixel # " << channelId << " col " << col << " row " << row
             << " ==> OVER threshold @ time " << iBin
             << " sum: " << boxcarsum << ", thres: " << pixelthresh
             << endl;
      }
      // END DEBUG ==========================================================================

    } else { // if (bin_above_threshold)

      if (pulse) { // if first bin below threshold after FLT pulse

        pulse = false;                // prepare next FLT pulse
        prolongationCounter = 0;      // FLT prolongation counter (start)
        tend = iBin-1;                // end of FLT pulse
      }

    } // if (bin_above_threshold) else


    if (hasFLT &&
        prolongationCounter<FLTprolongation) {

      if (bin_above_threshold) {

        if(!fMultiplicity.count(iBin)) {
          fMultiplicity[iBin] = 0;         // count FLTs in this time bin
        }
        fMultiplicity[iBin]++;             // count FLTs in this time bin

        fFLT[channelId][iBin] = 1;

      } else {

        fFLT[channelId][iBin] = 2;

      }

    }

    prolongationCounter++;

  } // end loop over tslots


  // DEBUG ==========================================================================
  if (fVerbosity>=4 && hasFLT) { // by ralf
    cout << "##### FLT summary channelId=" << channelId << " start=" << tstart << " length=" << tend-tstart+1
         << " end=" << tend
         << " FLTtrace(start-1|end+1): ";
    unsigned int i;
    for (i=tstart-1; i<=tend+1; i++) {
      if (fFLT[channelId].count(i)) cout << "1";
      else cout << "0";
    }
    for (i=tend; i<trace.GetSize(); i++) {
      if (!fFLT[channelId].count(i)) break;
    }
    cout << " first_0_bin=" << i
         << endl;
  } else  if (fVerbosity >= 2 && hasFLT) { // by sergio
    cout << "##### FLT start = " << tstart << " length = " << tend-tstart+1
         << ", end time = " << tend
         << endl;
  }
  // END DEBUG ==========================================================================


  if (ntslots) {
    meanBoxcarSum /= ntslots;
    rmsBoxcarSum = std::sqrt(double(squareBoxcarSum - ntslots * meanBoxcarSum * meanBoxcarSum)/(ntslots-1));
  }

  // remember all triggered pixels
  if (hasFLT) {
    fChannelFLT.push_back(channelId);
  }

  pixelSim.SetMaxBoxcarsum(maxBoxcarSum);
  pixelSim.SetFltTime(tstart);
  pixelSim.SetFltDuration(tend-tstart-1);
  pixelSim.SetMeanBoxcarsum(meanBoxcarSum);
  pixelSim.SetRmsBoxcarsum(rmsBoxcarSum);

  return hasFLT;
} // end of T1Trigger



int
FdTriggerSimulator::T2Trigger(const fdet::Telescope& detTel, unsigned int lastFLTbin)
{
  //const fdet::FDetector& detFD = Detector::GetInstance().GetFDetector();

  //if (fMultiplicity.size()<1) {
  //return 0;
  //}

  // Algorithm to select pixels

  int timeT2 = 0;

  const int pebit      = 0;
  const int triggerbit = 0;
  const int sparebit   = 1;

  // this is just for debugging
  map<int,int> sltMultiplicity;
  vector<int> usedSltPixels;


  // RESET SLT memory
  fSLT.clear();


  const double fltBinsPer50ns = detTel.GetChannel(1).GetFADCBinSize() / (50.*ns);
  int numT2bins = int(lastFLTbin * fltBinsPer50ns);


  // 50 ns SLT time slices : loop over time (2000*2 = 4000 bins)
  for (int tt2=0; tt2<numT2bins; ++tt2) {

    int tt2_100  = tt2/2;           // convert to 100ns binning (units of 100 ns (FLT))
    int flt_index = int(tt2/fltBinsPer50ns);
    int tt2_1000 = tt2/20;          // convert to 1000ns binning
    int tt2_1000_comp = tt2_1000;

    // define SLT window for pattern recognition
    const int sltWindow = 5;
    int col_cycle_load         = (tt2%20) + 1;     // define column of SLT memory to be refreshed
    int col_cycle_compare_low  = col_cycle_load - sltWindow + 1;    // start compare window
    if (col_cycle_compare_low<1) col_cycle_compare_low += 20;
    int col_cycle_compare_up   = col_cycle_compare_low + sltWindow - 1; // end compare window


    // column overflow for pattern recognition
    if (col_cycle_compare_up>20) {
      col_cycle_compare_up = 20;
      tt2_1000_comp -= 1;
      if (tt2_1000_comp<0) tt2_1000_comp = 0; // there are no FLTs anyway ...
    }





    // ++ Load SLT column col_cycle_load at time  tt2 ++

    int npixtot = 0;
    /* THIS DOES NOT TAKE INTO ACCOUNT PROLONGATED FLT TRACES !!!!!!!!!!
      if (fMultiplicity.count(tt2_100)) {
      npixtot = fMultiplicity[tt2_100];     // number of FLT pixels at tt2
      }
    */
    for (map<int,FltTrace>::const_iterator iFltTrace = fFLT.begin();
         iFltTrace != fFLT.end();
         ++iFltTrace) {
      const FltTrace& fltTrace = iFltTrace->second;
      if (fltTrace.count(flt_index)) npixtot++;
    }


    int rowMask = 0;


    // Debug =====================================================
    if (fVerbosity>=3) {
      cout << " SLT 50ns bin: " << setw(4) << tt2
           << ", /2: " << setw(4) << tt2_100
           << ", /20: " << setw(3) << tt2_1000
           << ", loadCol: " << setw(2) << col_cycle_load
           << ", compCol: " << setw(2) << col_cycle_compare_low
           << " - " << setw(2) << col_cycle_compare_up
           << " compBin1000: " << setw(3) << tt2_1000_comp
           << " nFLT: " << setw(3) << npixtot
           << " mask: ";
      if (npixtot==0) {
        if (fVerbosity>=4) { // by ralf
          cout << "0000000000000000000000.";
        } else { // by sergio
          cout << "0000000000000000000000";
        }
      }
    }

    // END Debug =====================================================



    if (npixtot>0) {

      // loop row
      for (int row=1; row<=22; ++row) {

        int channelMask = (1<<(row-1));
        int channelId = 22*(col_cycle_load-1) + row;

        if (fFLT.count(channelId) &&
            fFLT[channelId].count(flt_index)) {
          rowMask |= channelMask;
        }

        // Debug =====================================================
        if (fVerbosity>=3) {
          if (fFLT.count(channelId) &&
              fFLT[channelId].count(flt_index)) {
            if ((fVerbosity>=4) &&
                (fFLT[channelId][flt_index]==2)) cout << "2";
            else cout << "1";
          } else cout << "0";
        }
        // END Debug =====================================================


      }  // loop row

    } // if npixtot>0


    if (!fSLT.count(tt2_1000)) {
      fSLT[tt2_1000] = SLTData(detTel.GetLastColumn());
    }

    fSLT[tt2_1000].SetRowMask(col_cycle_load, rowMask);
    fSLT[tt2_1000].SetParityError(col_cycle_load, pebit);
    fSLT[tt2_1000].SetTrigger(col_cycle_load, triggerbit);
    fSLT[tt2_1000].SetSpare(col_cycle_load, sparebit);





    // ++ collect SLT channels in window ++

    set<int> pixInCycle;    // loop channel matrix 5x22
    for (int col=col_cycle_compare_low; col<=col_cycle_compare_up; col++) {
      for (int row=1; row<=22; ++row) {

        int channelId = 22*(col-1) + row;

        if (fSLT.count(tt2_1000_comp) &&
            fSLT[tt2_1000_comp].HasPixel(channelId)) {

          pixInCycle.insert(channelId); // retrieve channel within current cycle
                                        // (sorted by column!)

          // Debug ===================================================== by sergio
          if (fVerbosity ==2 && npixtot> 0)
            cout <<setw(2)<<pixInCycle.size()-1
                 <<". pixel "<<setw(3)<<channelId
                 <<"  col " <<setw(3)<<col<<"  row "<<setw(3)<<row<<endl;
          // END Debug =====================================================


        }  // if FLT
      }  // loop row
    } // loop col


    // Debug =====================================================
    if (fVerbosity>=3)
      cout << ", nPixCycle: " << setw(3) << pixInCycle.size()<<endl;
    if (fVerbosity >=2 && npixtot> 0)
      cout<<" @ SLT time "<<setw(3)<<tt2<<" col ["<<setw(2)<<col_cycle_compare_low
          <<"-"<<setw(2)<<col_cycle_compare_up<<"] -> there are "<<setw(2)
          <<npixtot<<" pixels of which "<<setw(2)<<pixInCycle.size()<<" within cycle"<<endl;
    if (fVerbosity>=200) {
      cout << " pixincycle list: ";
      for (set<int>::iterator cyclePixel = pixInCycle.begin();
           cyclePixel!=pixInCycle.end();
           ++cyclePixel) {
        cout << *cyclePixel << " ";
      }
    }
    // END Debug =====================================================



    // ####### SLT GEOMETRY ####### ++++= pattern search ++++++

    unsigned int  pattern = 0;
    int multipattern = 0;
    int colPrev = 0;

    // start loop over pixels within current SLT cycle
    for (set<int>::iterator cyclePixel = pixInCycle.begin();
         cyclePixel!=pixInCycle.end();
         ++cyclePixel) {

      int channelId = *cyclePixel;          // p[0] is the origin of the pattern
      int col = (channelId-1)/22 + 1;
      int row = (channelId-1)%22 + 1;;

      if (col != colPrev) multipattern = 0;
      colPrev = col;

      // SLT REQUIRES source pixel at 1st COLUMN
      if (col!=col_cycle_compare_low)
        continue;


      // SLT Pattern search
      unsigned int SLTpattern = TestPattern(pixInCycle, col, row,
                                            multipattern,
                                            usedSltPixels);
      if (SLTpattern>0) pattern = SLTpattern;

      // set SLT time
      if (pattern>0) {

        if (timeT2==0)
          timeT2 = tt2-(sltWindow-1);      // set SLT time T2 (start of SLT window!!)

        if (fSLT.count(tt2_1000_comp)) {
          fSLT[tt2_1000_comp].SetSLTPattern(col, pattern);
        } else {
          //ERROR ("MISSING SLT DATA ?????!!!!! ");
        }

        // Debug ===================================================== only if pixInCycle.size()>3
        //if (fVerbosity>=6)
        //cout << " -> pattern " << setw(3) << pattern;
        // END Debug =====================================================

      } // pattern>0

    } // end loop pixels in cycle

    // Debug =====================================================
    if (fVerbosity >= 1 && fVerbosity <3 && pattern > 0 && pixInCycle.size()>3)
      cout <<"=-=-=-=-> SecondLevelTrigger at time "<< setw(4) <<tt2
           <<" col "<< setw(2) <<col_cycle_compare_low<<" with pattern "<<setw(3) <<pattern<<endl;
    // END Debug =================================================

    // Debug bits ================================================
    if (fVerbosity >= 2 && pattern > 0 && pixInCycle.size()>3) {
      for (int col=col_cycle_compare_low; col<=col_cycle_compare_up; ++col) {
        cout << " ==> col " << setw(2) << col << " pixel bits: ";
        for(int b=0; b<22; ++b){
          unsigned int mask = (1<<b);
          if(fSLT[tt2_1000_comp].GetRowMask(col) & mask) cout << "1";
          else cout << "0";
        }
        if (col == col_cycle_compare_low) {
          cout << " -+-+-+- pattern: " << setw(3) << pattern << " = ";
          for(int b=23; b<30; ++b){
            unsigned int mask = (1<<b);
            if(fSLT[tt2_1000_comp].GetSLTDataWord(col) & mask) cout << "1";
            else cout << "0";
          }
        }
        cout << endl;
      } // end loop col
    } // END Debug bits =========================================

    // Debug =====================================================
    //if (fVerbosity>=5) {
    // cout << endl;
    //}
    // END Debug =====================================================

    if (col_cycle_load == 1) usedSltPixels.clear();
    else if (col_cycle_load == 20) sltMultiplicity[tt2_1000] = usedSltPixels.size();


  } // end loop over times tt2





  // Fill SLTWord at timeT2 (50ns)

  int timeT2_1000 = timeT2/20;
  int timeT2_100  = timeT2_1000*10;
  timeT2          = timeT2_1000*20;


  // SLT Debug =====================================================
  int SLT_T2[20];
  for (int tr_col=1; tr_col<=20; tr_col++) {
    if (fSLT.count(timeT2_1000)) {
      SLT_T2[tr_col-1] = fSLT[timeT2_1000].GetSLTDataWord(tr_col);
    }
  }
  if (fVerbosity>=1 && timeT2_100>0) {
    cout << endl
         << " **** SLT T2 time (units of 100 ns) is: " << timeT2_100
         << " (50ns: " << timeT2 << ", 1000ns: " << timeT2_1000 << ")"
         << " first SLT pattern was: "
         << endl
         << "                          ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ "
         << endl << endl;
    for ( int col=1; col<=20; ++col) {
      cout << " ==> col " << setw(2) << col << "  pixel bits: ";
      for(int b=0;b<22;++b){
        unsigned int mask = (1<<b);
        if(SLT_T2[col-1] & mask) cout << "1";
        else cout << "0";
      }
      unsigned int pattern =  SLT_T2[col-1] & (0x03f800000);
      pattern = pattern >> 23;
      cout << " -+-+-+- pattern: " << setw(3) << pattern << endl;
    }
    cout << endl
         << "-+-+- SLT multiplicity at each 1 mus ===>" << endl;
    for (int t=0; t<numT2bins/20; ++t) {
      int sltM = 0;
      if (sltMultiplicity.count(t)>0)
        sltM = sltMultiplicity[t];
      cout << setw(4) << sltM << ", ";
      if ((t+1)%20 ==0) cout << endl;
    }
    cout << endl;
  }
  // END SLT Debug =====================================================


  return timeT2_1000;

}// end of T2Trigger




unsigned int FdTriggerSimulator::TestPattern(const set<int>& pixInCycle, int col, int row,
                                             int& multipattern,
                                             vector<int>& usedInSlt) {

  if (pixInCycle.size()<=3) {
    return 0;
  }

  int pixelnumber;
  int pixels;
  unsigned int pattnum = 0;


  int   p[4]; // pixels numbers for this pattern
  short r[4]; // row number
  short c[4]; // col number
  p[0] = 22*(col-1) + row; //p[0] is the origin of the pattern
  c[0] = ((p[0] - 1) / 22) + 1;
  r[0] = ((p[0] - 1) % 22) + 1;

  // find candidate pixels of the pattern, given the starting pixel

  for (unsigned int patt=1; patt<109; ++patt) { // start loop over 108 patterns

    //usedInSlt.clear();
    pixels = 1;

    for( int m=1; m<4; ++m ) { //start loop over neighbours

      if( (row%2)==0){ // even row
        p[m] = p[0] + kSltDiffEven[patt-1][m-1]; //abs pixel number
        pixelnumber=p[m];
      }
      else {            // odd row
        p[m] = p[0] + kSltDiffOdd[patt-1][m-1]; // abs pixel number
        pixelnumber=p[m];
      }

      c[m] = ((p[m] - 1) / 22 ) + 1;  //col for this pixel
      r[m] = ((p[m] - 1) % 22 ) + 1;  //row for this pixel

      if (abs(r[0] - r[m])>4) pixelnumber = 0; // to avoid patterns across row 1->22

      if (pixInCycle.find(pixelnumber)!=pixInCycle.end()) {
        pixels++;
      }

    }// end loop over neighbours


    if(pixels==4) {

      pattnum = patt;
      multipattern++;
      for( int m=0; m<4; ++m ) {
        if (find(usedInSlt.begin(), usedInSlt.end(), p[m])==usedInSlt.end())
          usedInSlt.push_back(p[m]);
      }

      // Debug  ================================================
      if (fVerbosity>=6) { // by ralf
        cout << " \n       ### FOUND PATTERN no. " << setw(3) << patt
             << " = [";
        for (int dbgM=1; dbgM<4; dbgM++) {
          if( (row%2)==0) { // even row
            cout << kSltDiffEven[patt-1][dbgM-1]; //abs pixel number
          } else{            // odd row
            cout << kSltDiffOdd[patt-1][dbgM-1]; // abs pixel number
          }
          cout << ", ";
        }
        cout << "] = {" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << "}"
             << " from col " << c[0] << " row " << r[0] << " - mult = " << multipattern
             << endl;
      } else if (fVerbosity>=2) { // by sergio
        cout << " ### FOUND PATTERN n." << setw(3) << patt
             << " = ["<< setw(3) << p[0] << "," << setw(3)<< p[1] << ","
             << setw(3)<< p[2] << "," << setw(3)<< p[3] << "]"
             << " from col " << setw(2)<< c[0] << " row "
             << setw(2)<< r[0] << " - mult = " << multipattern << endl;
      }
      if (fVerbosity>=7) {
        for( int testCol=0; testCol<4; ++testCol) {
          cout << "        col=" << col+testCol << "   ";
          for( int testRow=1; testRow<=22; ++testRow) {
            int testPxlId = (col+testCol-1)*22 + testRow;
            if (pixInCycle.find(testPxlId)==pixInCycle.end()) cout << ".";
            else {
              bool testIsPat = false;
              for(int testM=1; testM<4; ++testM) { // start loop over neighbours
                int testPxlIdSLT = 0;
                if( (row%2)==0) { // even row
                  testPxlIdSLT = p[0] + kSltDiffEven[patt-1][testM-1]; //abs pixel number
                } else{            // odd row
                  testPxlIdSLT = p[0] + kSltDiffOdd[patt-1][testM-1]; // abs pixel number
                }
                if (testPxlId==testPxlIdSLT) testIsPat = true;
              }
              if (testIsPat || (col+testCol==col && testRow==row)) cout << "+";
              else cout << "1";
            }
          }
          cout << "\n";
        }
      }
      // END Debug  ================================================


    }


  } // end loop over patterns

  if (multipattern>1) { // if exist MULTIPATTERN
    pattnum = 127;
  } // end if MULTIPATTERN

  return pattnum;

}// End of TestPattern




TMirrorEvent*
FdTriggerSimulator::MakeMirrorEvent(fevt::Telescope& tel,
                                    const evt::Event& event,
                                    const TimeInterval& offsetFdSd,
                                    const int sltTimeShift)
{
  const fdet::Telescope& detTel = Detector::GetInstance().GetFDetector().GetTelescope(tel);
  const fdet::Camera& detCamera = detTel.GetCamera();
  //const int SLTbin = detCamera.GetSLTTriggerBin();
  int SLTbin;
  if (pageCounter == 1) {
    SLTbin = 280;
  }else{
    SLTbin = 0;
   }

  const double binsize = detCamera.GetFADCBinSize();
  const unsigned int tracelength = detCamera.GetFADCTraceLength();

  tel.SetStatus(fevt::ComponentSelector::eHasData);

  const fevt::FEvent& theFEvent = event.GetFEvent();
  const int eyeId = tel.GetEyeId();
  const int telId = tel.GetId();


  tel.MakeTriggerData();
  fevt::TelescopeTriggerData& telTrigger = tel.GetTriggerData();


  // RU Fri May 13 11:09:20 CEST 2005
  TMirrorEvent *MirrorEvent = new TMirrorEvent();

  // Objects to build the event
  TMirrorEventHeader    *header    = MirrorEvent->GetEventHeader();
  TMirrorPixelData      *pixelData = MirrorEvent->GetPixelData();
  TMirrorPixelList      *pixelList = MirrorEvent->GetPixelList();
  TMirrorFADCData       *fadcData  = MirrorEvent->GetFADCData();

  // EventHeader Info
  header->SetEventLabel( TMirrorEventHeader::kUnlabelled ); // classifcation is done by TLT (kShowerEvent)
  //header->SetEventLabel( TMirrorEventHeader::kShowerEvent );
  header->SetEventNo( theFEvent.GetHeader().GetId() );
  header->SetEventType( TMirrorEventHeader::kSimulatedEvent );
#if FDEVENTLIB_VERSION_CODE < ModuleVersionCode(4,0,0)
  header->SetTriggerNum( 3000 );
#else /* pre-v4 */
  header->SetTriggerNo( 3000 );
#endif /* FDEventLib at least as recent as v4 */
   /*
   * SLT page size is 1000 for old and new (HEAT) electronics.
   * TODO: Should be in FModelConfig?
   */
  const int onePage = 1000;

  const int nextPageDelay = onePage - SLTbin;


  // Time assigned in the LightAtDiaphragmSimulator as the begin of the
  // photon trace. This time represents the time of the first photons at the
  // apperture of the telescope
  // WARNING: In FdElectronicsSimulator the start time of the "prolonged" ADC
  //          traces (PhotonStartTime-(SLTbin*binsize))
  TimeStamp timeEvent = tel.GetSimData().GetPhotonsStartTime();
  const double eventStart = (sltTimeShift + (onePage-SLTbin)) * 100.*ns;  // used to be FADC bin size, but that's not correct! (This is SLT!)
  timeEvent += TimeInterval(eventStart);
  timeEvent -= offsetFdSd;


  ostringstream info;
  info << "Creating mirror event";
  info << " eye=" << eyeId << " tel=" << telId
       << " run=1 event=" << theFEvent.GetHeader().GetId()
       << " tStart=" << eventStart/ns << " ns";
  INFO(info);



  UInt_t nextPageNanoTime = timeEvent.GetGPSNanoSecond();
  UInt_t nextPageSec = timeEvent.GetGPSSecond();

  header->SetNextPageNanoTime (nextPageNanoTime);
  header->SetNextPageTime (nextPageSec);
  header->SetTriggerSource( TMirrorEventHeader::kInternalTrigger );
  header->SetTimeCorrectionParameters(100.0, 0.0);

  header->SetEyeNo( eyeId );
  header->SetMirrorNo( telId );
  header->SetNextPageDelay( nextPageDelay );

  //#warning ATTENTION! Run number in the FDAS output file is set to 1 (one).  Not suitable for blind analysis
  header->SetRunNo( 1 );

  header->SetT3Id( 4097 );
  header->SetTimeStampHigh( header->GetGPSTime() );  // TODO is this needed ?
  header->SetTimeStampLow( 0 );

  //MirrorEvent->GetEventHeader()->Dump();


  // Real signal

  // PixelList Info
  pixelList->SetNumPixels(0);

  Fd::PixelNumberRec pixelObj;
  Fd::PixelNumber pixel = &pixelObj;
  Fd::SetEyeNo(pixel, eyeId);
  Fd::SetMirrorNo(pixel, telId);

  for (list<int>::iterator iReadOutChannel = fChannelReadOutList.begin();
       iReadOutChannel != fChannelReadOutList.end();
       ++iReadOutChannel) {

    const unsigned int channelId = *iReadOutChannel;
    // Debug  ================================================
    if (fVerbosity>=11)
      cout << " adding trace eyeId=" << eyeId
           << " telId=" << telId
           << " pixelId=" << channelId
           << " (eyePixId=" << Fd::GetEyePixelNo( pixel ) << ")";
    // END Debug  ================================================

    const fdet::Channel& detChannel = detTel.GetChannel(channelId);
    if (!detChannel.IsVirtual()) {
      const unsigned int pixelId = detChannel.GetPixelId();
      if ( tel.HasPixel(pixelId, fevt::ComponentSelector::eDeSelected)
           && tel.GetPixel(pixelId, fevt::ComponentSelector::eDeSelected).GetStatus()
              == fevt::ComponentSelector::eDeSelected )
      {
        if (fVerbosity >= 11) {
          cout << " pixel is eDeSelected, status=";
          if (tel.HasPixel(pixelId, fevt::ComponentSelector::eExists))
            cout << tel.GetPixel(pixelId, fevt::ComponentSelector::eExists).GetStatus();
          else
            cout << "doesn't exist";
          cout << endl;
        }
        continue;
      }
    }


    // eye_no and mirror_no will not change
    Fd::SetPixelNo(pixel, channelId);
    pixelList->AddPixel(pixel);


    // FADCData info

    // union which contains eye, mirror and pixel number
    //Fd::PixelNumber pnumber = pixelList->GetPixel(i);
    // pnumber is used to select from internal 22*24 matrix
    TFADCData* data = fadcData->GetPixelFADCData(pixel);
    // access to the first bin, others are word[1] ... [999]
    TFADCData::FADCDataWord word = data->GetFADCTrace();
    data->SetPixelNumber(pixel);   // v2r0

    // The hardware trace starts at bin 0 and ends at bin 999 (1000 bins)
    // The number of bins readout can be limited around the peak of interest
    // as it is done for the calibration
    //const fdet::Channel& detChannel = detTel.GetChannel(channelId);

    data->SetTraceStartBin(0);
    data->SetTraceEndBin(tracelength-1);
    if (tracelength == 1000)
      data->SetTraceResolution(TFADCData::kFADC_RESOLUTION_1000_BINS);
    else if (tracelength == 2000)
      data->SetTraceResolution(TFADCData::kFADC_RESOLUTION_2000_BINS);
    // 4000 is not tested and thus not supported by FDEventLib
    //else if (tracelength == 4000) {
    //  data->SetTraceResolution(TFADCData::kFADC_RESOLUTION_4000_BINS);
    else {
      ostringstream errMsg;
      errMsg << "FADC trace length of " << tracelength
             << " has no correspondance in FDEventLib";
      ERROR(errMsg);
      throw OutOfBoundException(errMsg.str());
    }

    // set mean, variance and threshold
    unsigned int Mean = 0;
    unsigned int RMS = 0;
    unsigned int actualThr = 0;
    if (channelId<=detTel.GetLastPixelId()) {
      // mapped pixel Id
      const unsigned int pixelId = detTel.GetChannel(channelId).GetPixelId();
      if (tel.HasPixel(pixelId)) {
        const fevt::PixelSimData& pxsimdata = tel.GetPixel(pixelId).GetSimData();
        Mean = int(ceil(pxsimdata.GetMean()));
        RMS = int(ceil(pxsimdata.GetRMS()));
        actualThr = pxsimdata.GetThreshold();
      }
    }

    data->SetVariance(RMS);
    data->SetMean(Mean);
    data->GetPixelMonitorData()->fActualThreshold = actualThr;


    // Debug  ================================================
    if (fVerbosity>=11)
      cout << " rms=" << RMS
           << " mean=" << Mean
           << " threshold=" << actualThr;
    // END Debug  ================================================

    if (!tel.HasChannel(channelId)) { // skip non-sim pixels

      ostringstream err;
      err << " Channel/pixel with id=" << channelId << " was listed for FADC readout, "
          << " but there is no object like this available!"
          << endl;
      ERROR(err.str());

    } else {

      const fevt::Channel& ch = tel.GetChannel(channelId);
      const fevt::ChannelSimData& channel_sim = ch.GetSimData();
      const TraceI& fadc_trace = channel_sim.GetFADCTrace(fevt::FdConstants::eTotal);
      const int nFADCbins = fadc_trace.GetSize();

      // Debug  ================================================
      if (fVerbosity>=11)
        cout << " simTraceSize=" << nFADCbins << " trace:" << endl;
      // END Debug  ================================================

      const double sltBinsPerFltBin = detCamera.GetFADCBinSize() / (100.*ns);
      const int fltTimeShift = int(sltTimeShift / sltBinsPerFltBin);

      for(unsigned int k=0; k<tracelength; ++k ) {

        // Accessing via bits.fadcData is dangerous, names will change in v2r0
        // better use FADCDataWordSetData(&word[k],p0_data[k]) which is just a
        // macro and as fast as the direct access but portable.
        // I didn't provide C++ accessors for these elements simply for
        // speed reasons.
        FADCDataWordSetWord( &(word[k]), 0 );   // 16 bits of the union FDAS

        const int tbin = k+fltTimeShift;

        unsigned short int value = 0;
        if (tbin>=nFADCbins) {
          ostringstream err;
          err << "Bin=" << tbin << " out of simulated FADC trace min=" << 0
              << " max=" << nFADCbins-1 << " THIS SHOULD NEVER HAPPEN!"
              << endl;
          ERROR(err.str());
        }
        else value = fadc_trace[tbin];

        // Debug  ================================================
        if (fVerbosity>=12)
          cout << "    trace bin: " << setw(4) << k
               << " shifted " << setw(4) << tbin
               << " value=" << setw(5) << value;
        // END Debug  ================================================


        if (fFLT.count(channelId) &&
            fFLT[channelId].count(tbin)) {

          if (fFLT[channelId][tbin] == 1) {

            // trigger bit
            value |= (1<<15);

            // Debug  ================================================
            if (fVerbosity>=12)
              cout << "    FLT bit set";
            // END Debug  ================================================

          } else if (fFLT[channelId][tbin] == 2) {
            // Debug  ================================================
            if (fVerbosity>=12)
              cout << "    (FLT prolongated)";
            // END Debug  ================================================

          }

        }

        FADCDataWordSetWord( &(word[k]), value );   // FDAS

        // Debug  ================================================
        if (fVerbosity>=12)
          cout << endl;
        // END Debug  ================================================

      } // loop fadc bins
    } // has channel

    delete data; // RU: THIS IS NOT DONE AUTOMATICALLY !!!

  } // loop read out pixels

  // Filling TMirrorPixelData
  // I assume that all bits of the TMirrorPixelData are cleared
  // ---> ALSO DO SHIFT TO SLT TRIGGER TIME !!!
  vector<fevt::SLTData> sltData;
  const unsigned int nSltBins = detCamera.GetSLTTraceLength();
  for (unsigned int t2bin=0; t2bin<nSltBins; ++t2bin ) {     // 100 bins of 1mu s each

    const int t2binShift = t2bin + sltTimeShift/10;

    if (fVerbosity>101) {
      cout << " add slt t2bin=" << t2bin << " t2binShift=" << t2binShift << " \n";
    }

    fevt::SLTData sltDataCol(20);
    for (unsigned int col=1; col<=20; ++col) {
      TMirrorPixelData::PixelDataWord pword = pixelData->GetPixelData( col, t2bin );
      int slt_word = 0;
      if (fSLT.count(t2binShift)) {
        slt_word = fSLT[t2binShift].GetSLTDataWord(col);
      }
      PixelDataSetWord( pword, slt_word);        // FDAS
      sltDataCol.SetSLTDataWord(col, slt_word);  // Offline
    }
    sltData.push_back(sltDataCol);
  }
  telTrigger.SetSLTData(sltData);

  // Filling of multiplicity
  // Multiplicity is always 1000 words
  const double multiplicityBinsPerFltBin = (100.*ns) / detCamera.GetFADCBinSize();
  telTrigger.MakeMultiplicity(1000, binsize);
  TraceI& multiplicityTrace = telTrigger.GetMultiplicity();
  TMirrorPixelData::MultiplicityDataWord multiplicityData = pixelData->GetMultiplicityData();
  for (unsigned int tbin = 0; tbin < 1000; ++tbin) {

    // DKH: For HEAT, every second fMultiplicity bin is taken
    const int tbinShift = (tbin + sltTimeShift) * multiplicityBinsPerFltBin;

    int multiplicity = 0;
    // Note that fMultiplicity is in the same units as the FADC traces (50ns for HEAT)
    if (fMultiplicity.count(tbinShift))
      multiplicity = fMultiplicity[tbinShift];

    multiplicityData[tbin] = multiplicity;  // FDEventLib
    multiplicityTrace[tbin] = multiplicity; // Offline
  }

  return MirrorEvent;

} // end of AddMirrorEvent



bool
FdTriggerSimulator::AddEyeEvent(fevt::Eye& eye, evt::Event& event,
                                TEyeEvent& EyeEvent)
{
  eye.SetStatus(fevt::ComponentSelector::eHasData);

  eye.MakeHeader();

  fevt::FEvent& theFEvent = event.GetFEvent();
  fevt::EyeHeader& headerEye = eye.GetHeader();
  fevt::Header& headerFEvent = theFEvent.GetHeader();
  evt::Header& evtHeader = event.GetHeader();

  const int eyeId = eye.GetId();
  const int eventId = theFEvent.GetHeader().GetId();
  const int runId = 1;

  TEyeEventHeader *eyeheader = EyeEvent.GetEventHeader();

  eyeheader->SetEventNo(eventId);
  eyeheader->SetRunNo(runId);
  eyeheader->SetTimeCorrectionParameters(100, 0);
  //eyeheader->Set
    ////SetDeadTimeCounter(unsigned long long,UInt_t);
    //SetEyeDaqDeadTime(unsigned long long dead_time)
    //SetIsRejectedT3 CDAS veta
    //SetLidarIsTriggered()
    //SetMirrorDaqDeadTime(unsigned long long,unsigned int mirror_id);
    //SetMirrorNSkippedTriggers(unsigned int,unsigned int mirror_id);
    //void SetMirrorPresent(UInt_t);
    //SetNRejectedT3(unsigned int n_rejected) by CDAS
    //SetTotalT3DeadTime(float total_time)


  // Added by T. Porter - need to fill some more of the header parameters
  // based on current eye for HRec/HSim to work
  eyeheader->SetEyeNo(eyeId);

  // RawEvent part of the Event
  if (!event.HasRawEvent()) {
    INFO ("Making raw event for the current event");
    event.MakeRawEvent();
  }

  // this fills the AUGER raw event
  AugerEvent& rawEvent = event.GetRawEvent();
  rawEvent.PushEvent(EyeEvent);
  rawEvent.EventId = eventId;




  // Fill offline header
  unsigned int sec = eyeheader->GetTimeStamp()->GetGPSSec();
  unsigned int nsec = eyeheader->GetTimeStamp()->GetNanoSec();
  utl::TimeStamp tstamp(sec, nsec);

  //unsigned int runno = eyeheader->GetRunNo();
  //unsigned int evno  = eyeheader->GetEventNo();
  int eclass         = eyeheader->GetEventClass();
  int etype          = eyeheader->GetEventType();

  // fill evt::Header of Event
  const utl::TimeStamp& previous_time = evtHeader.GetTime();
  std::string previous_str = evtHeader.GetId();
  if (previous_time.GetGPSSecond() == 0 ||
      previous_time>tstamp) {
    event.GetHeader().SetTime(tstamp);
  }

  /*
    ostringstream id_str;
    if (previous_str!="")
    id_str << previous_str << "__";
    id_str << "eye" << eyeid
    << "_run" << runno
    << "_event" << evno;
    theEvent.GetHeader().SetId(id_str.str());
  */

  // Fill header of FEvent
  headerFEvent.SetNEyes(headerFEvent.GetNEyes()+1);

  // Fill header of eye
  headerEye.SetTimeStamp(tstamp);
  headerEye.SetEventNumber(eventId);
  headerEye.SetRunNumber(runId);
  headerEye.SetEventType(static_cast<fevt::EyeHeader::EventType>(etype));
  headerEye.SetEventClass(static_cast<fevt::EyeHeader::EventClass>(eclass));


  // set telescopes in DAQ bits
  for (fevt::Eye::ConstTelescopeIterator iTel = eye.TelescopesBegin(fevt::ComponentSelector::eInDAQ);
       iTel != eye.TelescopesEnd(fevt::ComponentSelector::eInDAQ);
       ++iTel) {
    eyeheader->SetMirrorPresent(iTel->GetId());
  }

  const fdet::FDetector& detFD = Detector::GetInstance().GetFDetector();
  const fdet::Eye& detEye = detFD.GetEye(eye);

  // copy the telescope time offset into the Offline data structure
  for (unsigned int iMirror = detEye.GetFirstTelescopeId();
       iMirror <= detEye.GetLastTelescopeId(); ++iMirror) {

    if (eyeheader->IsMirrorDataPresent(iMirror)) {

      if (!eye.HasTelescope(iMirror, ComponentSelector::eInDAQ))
        eye.MakeTelescope(iMirror, ComponentSelector::eInDAQ);

      fevt::Telescope& telescope = eye.GetTelescope(iMirror, ComponentSelector::eInDAQ);
      telescope.SetTimeOffset(eyeheader->GetMirrorTimeOffset(iMirror));

    } // if present in FDEventLib event
  }


  return true;

}// end of WriteEyeEvent


// Configure (x)emacs for this file ...
// Local Variables:
// mode:c++
// compile-command: "make -C .. -k"
// End:
