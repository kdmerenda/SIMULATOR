#ifndef _ELVESSimulationNS_ELVESSimulator_h_
#define _ELVESSimulationNS_ELVESSimulator_h_

/**
 * \file
 * \author Kevin-Druis Merenda
 * \date 06 Feb 2017
 */

#include <fwk/VModule.h>
#include <boost/utility.hpp>
#include <map>
#include <vector>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TCutG.h"
#include "TFile.h"
#include <utl/CoordinateSystemPtr.h>
#include <utl/ShowerParticleIterator.h>
#include <utl/ShowerParticleList.h>
#include <utl/VParticleProperties.h>
#include <utl/MultiTabulatedFunction.h>
#include <utl/Vector.h>
#include <utl/Point.h>
#include <utl/TimeStamp.h>
#include <utl/AugerUnits.h>
#include <utl/LameShadowPtr.h>
#include <utl/config.h>
#include <utl/RandomEngine.h>
#include <utl/ShadowPtr.h>
#include <fdet/Telescope.h>

/*
 * Avoid using using namespace declarations in your headers,
 * doing so makes all symbols from each namespace visible
 * to the client which includes your header.
 */




class ELVESSimulator : public boost::noncopyable, public fwk::VModule {

  enum Status {
    eTransformCoordinates,
    eGeneratePhotons
  };

  struct ELVESSimData { 
    float timeEye1;
    float timeEye2;
    float timeEye3;
    float timeEye4;
    float time;
    float nphotons;
    float nphotonsEye1;
    float nphotonsEye2;
    float nphotonsEye3;
    float nphotonsEye4;
    float nphotonsnormalized;
    float X;//in elves local cs
    float Y;
    float Z;
  };
  struct ELVESSimDataONE { 
    float timeEye1;
    float timeEye2;
    float timeEye3;
    float timeEye4;
    float X;//in elves local cs
    float Y;
    float Z;
    float nphotons;
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
      struct by_timeEye1ONE { 
	bool operator()(ELVESSimDataONE const &a, ELVESSimDataONE const &b) { 
	  return a.timeEye1 < b.timeEye1;
	}
      };
      struct by_timeEye2ONE { 
	bool operator()(ELVESSimDataONE const &a, ELVESSimDataONE const &b) { 
	  return a.timeEye2 < b.timeEye2;
	}
      };
      struct by_timeEye3ONE { 
	bool operator()(ELVESSimDataONE const &a, ELVESSimDataONE const &b) { 
	  return a.timeEye3 < b.timeEye3;
	}
      };
      struct by_timeEye4ONE { 
	bool operator()(ELVESSimDataONE const &a, ELVESSimDataONE const &b) { 
	  return a.timeEye4 < b.timeEye4;
	}
      };

      struct by_time { 
	bool operator()(ELVESSimData const &a, ELVESSimData const &b) { 
	  return a.time < b.time;
	}
      };
      
      struct by_nphotonsONE { 
	bool operator()(ELVESSimDataONE const &a, ELVESSimDataONE const &b) { 
	  return a.nphotons < b.nphotons;
	}
      };

      struct TimeCutIndices{
	int eye1;
	int eye2;
	int eye3;
	int eye4;
      };
      struct  OutputTreeVariables{
	float PeakCurrent;                // VFloat storing the value used for peak current. 
	float AltSource;                  //Variable deinfing the Altitude of the source that was simulated
	float BRadius,BTheta,BPhi;
	int NStepsTime,NStepsRadius,NStepsTheta,NStepsPhi;
	float SizeStepTime,SizeStepRadiusLow,SizeStepRadiusHigh,SizeStepTheta,SizeStepPhi; 
	int TimeIntegratedSteps;          //to save memory, simulation done at 500ns, but data integrated to 1000ns once simulated
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
  double fVODTot;
  int fNumTreeEntries;
  std::string fELVESInputName;
  std::string fELVESTreeName;
  std::string fELVESParameterTreeName;
  TTree* fTree;
  TTree* fTreeSim;
  TTree* fParameterTree;
  TFile* fIn;
  utl::ShadowPtr<utl::Point> fPosition;
  mutable utl::CoordinateSystemPtr fLocalCS;
  mutable utl::CoordinateSystemPtr fElvesCS;
  double fAzimuth;
  double fZenith;
  utl::RandomEngine* fRandomEngine; 
  int  nPhotonsToCreate;
  int  fNPhotonsToCreate;
  utl::TimeStamp fSimTime;
  /// Get the azimuth angle of the shower
  int fMaster = 0;
  OutputTreeVariables SimulationParameters;
  TTree* fSimDataTree;
  double fSimTimeStart;
  double fSimTimeEnd;
  int fEyeSelect;
  int fTelSelect;
  int fNPages;
  int fdogeomcorr;
  int fdoatmocorr;
  int fdoprecheck;
  int fdoradial;
  int fprodversion;
  double fNPhotonsMIN;
  double fNPhotonsMAX;
  double fNPhotonsMINEye1;
  double fNPhotonsMAXEye1;
  double fNPhotonsMINEye2;
  double fNPhotonsMAXEye2;
  double fNPhotonsMINEye3;
  double fNPhotonsMAXEye3;
  double fNPhotonsMINEye4;
  double fNPhotonsMAXEye4;
  ///Display progress
  void displayProgress(Int_t, Int_t, Int_t &);
  
  std::vector<ELVESSimData> ELVESData;
  std::vector<ELVESSimDataONE> ELVESDataONE;
  std::vector<TimeCutIndices> PhotonLoops100us;

  utl::CoordinateSystemPtr eye1CSSim;
  utl::CoordinateSystemPtr eye2CSSim; 
  utl::CoordinateSystemPtr eye3CSSim; 
  utl::CoordinateSystemPtr eye4CSSim; 
  utl::CoordinateSystemPtr eyeCSSim; 

  
  int fLoop;
  int fLoopSelect;
  bool fInit;
  Status fStatus;


  REGISTER_MODULE("ELVESSimulator",ELVESSimulator);
  
};


#endif // _ELVESSimulationNS_ELVESSimulator_h_
