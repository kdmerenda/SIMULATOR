#ifndef _ELVESSimulatorNS_ELVESSimulator_h_
#define _ELVESSimulatorNS_ELVESSimulator_h_

/**
 * \file
 * \author Kevin-Druis Merenda
 * \date 10 Apr 2017
 */

#include <fwk/VModule.h>
#include <boost/utility.hpp>
#include <map>
#include <vector>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include <utl/CoordinateSystemPtr.h>
#include <utl/ShowerParticleIterator.h>
#include <utl/ShowerParticleList.h>
#include <utl/VParticleProperties.h>
#include <utl/Vector.h>
#include <utl/Point.h>
#include <utl/AugerUnits.h>
#include <utl/RandomEngine.h>
#include <utl/TimeStamp.h>

/*
 * Avoid using using namespace declarations in your headers,
 * doing so makes all symbols from each namespace visible
 * to the client which includes your header.
 */

namespace ELVESSimulatorNS {
  
  class ELVESSimulator : 
    public boost::noncopyable,
    public fwk::VModule {
      
      enum Status {
	eTransformCoordinates,
	eGeneratePhotons
      };
	
      struct ELVESSimData { 
	double timeEye1;
	double timeEye2;
	double timeEye3;
	double timeEye4;
	double ndensity;
	utl::Point positions;
      };
      
      struct by_timeEye1 { 
	bool operator()(ELVESSimData const &a, ELVESSimData const &b) { 
	  return a.timeEye1 < b.timeEye1;
	}
      };
      struct by_timeEye2 { 
	bool operator()(ELVESSimData const &a, ELVESSimData const &b) { 
	  return a.timeEye2 < b.timeEye2;
	}
      };
      struct by_timeEye3 { 
	bool operator()(ELVESSimData const &a, ELVESSimData const &b) { 
	  return a.timeEye3 < b.timeEye3;
	}
      };
      struct by_timeEye4 { 
	bool operator()(ELVESSimData const &a, ELVESSimData const &b) { 
	  return a.timeEye4 < b.timeEye4;
	}
      };

      struct by_ndensity { 
	bool operator()(ELVESSimData const &a, ELVESSimData const &b) { 
	  return a.ndensity < b.ndensity;
	}
      };
      
      struct  OutputTreeVariables{
	float PeakCurrent;                // VFloat storing the value used for peak current. 
	float AltSource;                  //Variable deinfing the Altitude of the source that was simulated
	float BRadius,BTheta,BPhi;
	int NStepsTime,NStepsRadius,NStepsTheta,NStepsPhi;
	float SizeStepTime,SizeStepRadiusLow,SizeStepRadiusHigh,SizeStepTheta,SizeStepPhi; 
	int TimeIntegratedSteps;          //to save memory, simulation done at 100ns, but data integrated to 500ns once simulated
	int Orientation;                  //will stay an integer describing the tilt of the lightining bolt. 0 for vertical and 1 for horizontal
	int LightningType;                // 0 for CG, 1 for IC, 2 for CID
	int GroundMethod; // 0 perf, 1 sibc, 2 real
	int GWave;                        //activated = 1, deactivated = 0
      }; 

  public:
      ELVESSimulator();
      ~ELVESSimulator();
      VModule::ResultFlag Init();
      VModule::ResultFlag Run(evt::Event& e);
      VModule::ResultFlag ELVESSimDataCreator();
      VModule::ResultFlag Finish();
  private:
      double fELVESCenterLat;
      double fELVESCenterLon;
      int fNumDiaGridPoints;
      int fPhotonDiscr;
      int fNumTreeEntries;
      std::string fELVESInputName;
      std::string fELVESTreeName;
      std::string fELVESParameterTreeName;
      TTree* fTree;
      TTree* fParameterTree;
      TFile* fIn;
      mutable utl::CoordinateSystemPtr fLocalCS;
      mutable utl::CoordinateSystemPtr fElvesCS;
      double fAzimuth;
      double fZenith;
      RandomEngine* fRandomEngine; 
      int  fNPhotonsToCreate;
      utl::TimeStamp fSimTime;
      /// Get the azimuth angle of the shower                                                                         
      OutputTreeVariables SimulationParameters;
  
      ///Display progress
      void displayProgress(int, int, int&);
      
      std::vector<ELVESSimData> ELVESData;
      std::vector<int> TimeCutIndices;
      
      int fLoop;
      int fLoopSelect;
      bool fInit;
      Status fStatus;
      
    // This goes at the end.
      REGISTER_MODULE("ELVESSimulator",ELVESSimulator);
  };
}

#endif // _ELVESSimulatorNS_ELVESSimulator_h_
