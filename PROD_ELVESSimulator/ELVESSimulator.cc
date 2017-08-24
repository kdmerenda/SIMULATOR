#include "ELVESSimulator.h"
//#include "mpi.h"
#include <evt/Event.h>
#include <fevt/FEvent.h>
#include <fevt/Eye.h>
#include <fevt/TelescopeSimData.h>
#include <fevt/Telescope.h>
#include <fevt/PixelSimData.h>
#include <fevt/Pixel.h>
#include <fevt/FdConstants.h>
#include <det/Detector.h>
#include <fdet/FDetector.h>
#include <fdet/Eye.h>
#include <fevt/EyeHeader.h>
#include <fdet/Telescope.h>
#include <fdet/Camera.h>
#include <fdet/Pixel.h>
#include <fdet/Channel.h>
#include <fdet/Mirror.h>
#include <fdet/Filter.h>
#include <fdet/Corrector.h>

#include <utl/Point.h>
#include <utl/UTMPoint.h>
#include <utl/Particle.h>
#include <utl/Vector.h>
#include <utl/TimeStamp.h>
#include <utl/Particle.h>
#include <utl/ErrorLogger.h>
#include <utl/ReferenceEllipsoid.h>
#include <utl/RandomEngine.h>
#include <fwk/LocalCoordinateSystem.h>
#include <fwk/CoordinateSystemRegistry.h>
#include <fwk/CentralConfig.h>
#include <fwk/SVNGlobalRevision.h>
#include <fwk/RandomEngineRegistry.h>
#include <utl/Math.h>

#include <utl/MathConstants.h>
#include <utl/CoordinateSystemPtr.h>
#include <utl/AxialVector.h>
#include <utl/Vector.h>
#include <utl/PhysicalConstants.h>
#include <utl/AugerUnits.h>

#include <CLHEP/Random/Randomize.h>

#include <TMap.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include <TMath.h>
#include <TVector3.h>
#include <vector>

using namespace std;
using namespace utl;
using namespace fwk;
using namespace evt;
using namespace fevt;
using namespace fdet;
using namespace det;
using CLHEP::RandFlat;


ELVESSimulator::ELVESSimulator() 
{ }

ELVESSimulator::~ELVESSimulator()
{ }


VModule::ResultFlag 
ELVESSimulator::Init()
{
  INFO("ELVESSimulator::Init()");

  CentralConfig* cc = CentralConfig::GetInstance();
  Branch topB = cc->GetTopBranch("ELVESSimulator");

  
  if (!topB) {
    ostringstream err;
    err << "Missing configuration file for ELVESSimulator";
    ERROR(err);
    return eFailure;
  }

  fRandomEngine = &RandomEngineRegistry::GetInstance().Get(RandomEngineRegistry::eDetector);

  // Initialize the config file parameters. 
  topB.GetChild("ELVESInput").GetData(fELVESInputName);
  topB.GetChild("ELVESTreeName").GetData(fELVESTreeName);
  topB.GetChild("ELVESParameterTreeName").GetData(fELVESParameterTreeName);
  topB.GetChild("ELVESCenterLatitude").GetData(fELVESCenterLat);
  topB.GetChild("ELVESCenterLongitude").GetData(fELVESCenterLon);
  topB.GetChild("PhotonDiscretization").GetData(fPhotonDiscr);
  topB.GetChild("VODTot").GetData(fVODTot);
  topB.GetChild("FrameSelection").GetData(fLoopSelect);
  topB.GetChild("EyeSelection").GetData(fEyeSelect);
  topB.GetChild("FrameCount").GetData(fNPages);
  topB.GetChild("PreCheck").GetData(fdoprecheck);
  topB.GetChild("RadialAnalysis").GetData(fdoradial);
  topB.GetChild("ProdVersion").GetData(fprodversion);
  topB.GetChild("GeometricCorrection").GetData(fdogeomcorr);
  topB.GetChild("AtmosphericCorrection").GetData(fdoatmocorr);


  fIn = new TFile(fELVESInputName.data());
  fTree = (TTree*)fIn->Get(fELVESTreeName.data());
  //  if (fNumTreeEntries == 0) fNumTreeEntries = fTree->GetEntries();
  fNumTreeEntries = fTree->GetEntries();
  fSimTimeStart = fTree->GetMinimum("time");
  fSimTimeEnd = fTree->GetMaximum("time");
  
  fParameterTree = (TTree*)fIn->Get(fELVESParameterTreeName.data());
  fParameterTree->SetBranchAddress("parameters", &SimulationParameters);
  fParameterTree->GetEntry(0);


  ostringstream info;
  info << " Version: "
       << GetVersionInfo(VModule::eRevisionNumber) << "\n"
          " Parameters:\n"
          "                Input File: " << fELVESInputName << "\n"
          "                TTree Name: " << fELVESTreeName << "\n"
          "      Parameter TTree Name: " << fELVESParameterTreeName << "\n"
          "     Photon Discretization: " << fPhotonDiscr << "\n"
          "    Number of Tree Entries: " << fNumTreeEntries << "\n"
          "            Frame Selected: " << fLoopSelect << "\n"
          "                Page Count: " << fNPages << "\n"
          "               Do PreCheck: " << fdoprecheck << "\n"
          "        Do Radial Analysis: " << fdoradial << "\n"
          "        Production Version: " << fprodversion << "\n"
          " Do Atmospheric Correction: " << fdoatmocorr << "\n"
          "   Do Geometric Correction: " << fdogeomcorr << "\n"
          "                 Total VOD: " << fVODTot << "\n"
          "           Latitude Center: " << fELVESCenterLat/deg << "\n"
          "          Longitude Center: " << fELVESCenterLon/deg;
  INFO(info);



  fLoop = 0;
  fInit = false;
  fStatus = eTransformCoordinates;
  return eSuccess;
}

VModule::ResultFlag 
ELVESSimulator::Run(evt::Event& event)
{
  //INFO("ELVESSimulator::Run()");

  //this is to satisfy the module sequence design... fLoop = 0 for the coordinate transform, fLoop = fNPages for any pages saved. 
  if(fLoop > fNPages) return eBreakLoop;
  
  if(fStatus == eTransformCoordinates){
    INFO("Preparing ELVES Simulated Data Structure.");
    ELVESSimDataCreator();
    fIn->Close();    
  }
  
  if(fStatus == eGeneratePhotons){
    Detector& detector = Detector::GetInstance();
    const FDetector& detFD = detector.GetFDetector();
    
    if (!event.HasFEvent())
      event.MakeFEvent();
    if (!event.HasSimShower())
      event.MakeSimShower();
    
    
    fevt::FEvent& fEvent = event.GetFEvent();
    
    INFO("Creating FEvent");
    
    // set the FD event times
    fSimTime = detector.GetTime();
    event.GetHeader().SetTime(fSimTime);
         
    
    for (FDetector::EyeIterator iEye = detFD.EyesBegin();
	 iEye != detFD.EyesEnd() ; ++iEye) {
      const unsigned int eyeId = iEye->GetId();

      if(!fEvent.HasEye(eyeId))
	fEvent.MakeEye(eyeId, fevt::ComponentSelector::eInDAQ);
      fevt::Eye& eyeEvent = fEvent.GetEye(eyeId, fevt::ComponentSelector::eInDAQ);
      
      
      for (fdet::Eye::TelescopeIterator iTel = iEye->TelescopesBegin();
	   iTel != iEye->TelescopesEnd(); ++iTel) {
	const unsigned int telId = iTel->GetId();
	if(!eyeEvent.HasTelescope(telId))
	  eyeEvent.MakeTelescope(telId, fevt::ComponentSelector::eInDAQ);
	fevt::Telescope& telEvent = eyeEvent.GetTelescope(telId, fevt::ComponentSelector::eInDAQ);

	if(!telEvent.HasSimData())
	  telEvent.MakeSimData();
	fevt::TelescopeSimData& telSim = telEvent.GetSimData();
	
	//	telSim.SetPhotonsStartTime(fSimTime);
       	telSim.SetNumberOfPhotonBins(1000);
	
      }
    }  
  }

  if(fStatus == eGeneratePhotons){

    const ReferenceEllipsoid ellipsoid(ReferenceEllipsoid::Get(ReferenceEllipsoid::eWGS84));
    UTMPoint elvesWantedLocationUTM(fELVESCenterLat, fELVESCenterLon, 5.0*km, ellipsoid);
    const CoordinateSystemPtr WantedLocationCS =  fwk::LocalCoordinateSystem::Create(elvesWantedLocationUTM.GetPoint());
    const CoordinateSystemPtr referenceCS = Detector::GetInstance().GetReferenceCoordinateSystem();

    fevt::FEvent& fEvent = event.GetFEvent();
    Detector& detector = Detector::GetInstance();
    const FDetector& detFD = detector.GetFDetector();
    fEvent.GetHeader().SetTime(fSimTime+TimeInterval((fLoop-1)*100000));

    for (fevt::FEvent::EyeIterator iEye = fEvent.EyesBegin(fevt::ComponentSelector::eInDAQ),
	   end = fEvent.EyesEnd(fevt::ComponentSelector::eInDAQ);
	 iEye != end; ++iEye) {

      unsigned int eyeId = iEye->GetId();
      fevt::Eye& eyeEvent = *iEye;
      const CoordinateSystemPtr&  eyeCS = detFD.GetEye(eyeEvent).GetEyeCoordinateSystem();
      const CoordinateSystemPtr&  eyeCSLocal = detFD.GetEye(eyeEvent).GetLocalCoordinateSystem();
      // set the sim header id

      if(eyeId > 4)continue;
      ostringstream idStr;
      idStr << "eye" << eyeId <<"_run"<< "1" << "_event" <<  fLoop;
      cout << event.GetHeader().GetId() << endl;
      event.GetHeader().SetId(idStr.str());
      fEvent.GetHeader().SetId(fLoop);
      
      for (fevt::Eye::TelescopeIterator iTel = eyeEvent.TelescopesBegin(fevt::ComponentSelector::eInDAQ),
	     end = eyeEvent.TelescopesEnd(fevt::ComponentSelector::eInDAQ);
	   iTel != end; ++iTel) {
	
	const unsigned int telId = iTel->GetId();
	
	fevt::Telescope& telEvent = *iTel;
       	telEvent.SetTracesStartTime(fSimTime+TimeInterval((fLoop-1)*100000));
	const fdet::Telescope& detTel = detFD.GetTelescope(telEvent);
	
	const double rDia = detTel.GetDiaphragmRadius();

	fevt::TelescopeSimData& telSim = telEvent.GetSimData();
	telSim.SetPhotonsStartTime(fSimTime+TimeInterval((fLoop-1)*100000)); 

	const double normWavelength = detFD.GetReferenceLambda();
	//	const double fRDiaMin = rDia/10. ;//minimum radius on diaphragm, lets say a tenth of the radius
	const CoordinateSystemPtr& telCS = detTel.GetTelescopeCoordinateSystem();
	double detFOV = detTel.GetCamera().GetFieldOfView()/deg;//~23 deg
	detFOV += 10;//add 10 degrees to FOV so that more telescopes will be triggered. 

        double elvesTheta  = 90.-elvesWantedLocationUTM.GetPoint(referenceCS).GetTheta(referenceCS)/deg;
	double elvesPhi = elvesWantedLocationUTM.GetPoint(referenceCS).GetPhi(referenceCS)/deg;
	double elvesR = elvesWantedLocationUTM.GetPoint(referenceCS).GetR(referenceCS)/km;
	
	Vector telAxis = detTel.GetAxis();
	double telFOVLow =  telAxis.GetPhi(referenceCS)/deg - detFOV/2.;
	double telFOVHigh =  telAxis.GetPhi(referenceCS)/deg + detFOV/2.;
	
	//if we are outside detector field of view, go to next telescope, need to add elevation check too. 
	ostringstream info;  
	if (elvesPhi > telFOVLow && elvesPhi < telFOVHigh) {
	  if(!(fEyeSelect==0) && !(eyeId==fEyeSelect))   continue;
	  info << "Eye: " <<  eyeId << " Tel: " << telId << " - DETECTION!" << endl;
	  INFO(info);
	}else{
	  info << "Eye: " <<  eyeId << " Tel: " << telId << " - ELVES not in field of view";
	  INFO(info);
	  continue;//not in FOV so don't add photons for this telescope... 
	}
	
      //sort data for the given eye. and set the timing to be used below to the one of the right eye. Maybe try to find a faster way to do this. 
	INFO("Sorting Structure For Current Eye.");
  	if (eyeId == 1){
	  sort(ELVESData.begin(), ELVESData.end(), by_timeEye1());
          for (int iTime = 0; iTime < fNPhotonsToCreate; iTime ++){
	    ELVESData[iTime].time = ELVESData[iTime].timeEye1;
	    ELVESData[iTime].nphotonsnormalized = ELVESData[iTime].nphotonsEye1 ;   
	  }
	}
	if (eyeId == 2){
	  sort(ELVESData.begin(), ELVESData.end(), by_timeEye2());
	  for (int iTime = 0; iTime < fNPhotonsToCreate; iTime ++){
	    ELVESData[iTime].time = ELVESData[iTime].timeEye2;
	    ELVESData[iTime].nphotonsnormalized = ELVESData[iTime].nphotonsEye2 ;
	  }
	}
	if (eyeId == 3){
	  sort(ELVESData.begin(), ELVESData.end(), by_timeEye3());
	  for (int iTime = 0; iTime < fNPhotonsToCreate; iTime ++) {
	    ELVESData[iTime].time = ELVESData[iTime].timeEye3;
	    ELVESData[iTime].nphotonsnormalized = ELVESData[iTime].nphotonsEye3 ;
	  }
	}
	if (eyeId == 4){
	  sort(ELVESData.begin(), ELVESData.end(), by_timeEye4());
	  for (int iTime = 0; iTime < fNPhotonsToCreate; iTime ++){
	    ELVESData[iTime].time = ELVESData[iTime].timeEye4;
	    ELVESData[iTime].nphotonsnormalized = ELVESData[iTime].nphotonsEye4 ;
	  }
	}

	const unsigned int nBins = 1000;
	const double tracebin = 1e-7;
	telSim.ClearPhotons();
	if (!telSim.HasPhotonTrace(FdConstants::eFluorDirect))
	  telSim.MakePhotonTrace(FdConstants::eFluorDirect,1);
	if (!telSim.HasRayTracedPhotonTrace())
	  telSim.MakeRayTracedPhotonTrace(nBins, tracebin);
	telSim.ClearRayTracedPhotonTrace();


	//This is the creation of photons
	INFO("Creating Photons.");
	int photonCounter = 0;
	for(int iCell=0; iCell<fNPhotonsToCreate;iCell+=fPhotonDiscr){
	  
	  //check if we are in the window of interest.
	  if(!(ELVESData[iCell].time >= (fLoop-1)*0.0001 && ELVESData[iCell].time < (fLoop)*0.0001)) continue;

	  if (ELVESData[iCell].nphotonsnormalized < 1.0) continue;//is this fair?

	    const double diaPh = RandFlat::shoot(&fRandomEngine->GetEngine(), 0.0, kTwoPi);
	    const double diaR = rDia * sqrt(RandFlat::shoot(&fRandomEngine->GetEngine(), 0.2, 1.0));
	    const double xDia = diaR * cos(diaPh);
	    const double yDia = diaR * sin(diaPh);
	    Point pIn(xDia, yDia, 0.0, telCS);
	 	      
	    double photonTime=(ELVESData[iCell].time  - (fLoop-1)*100e-6)*1e9;
	    
	    Point pElves(ELVESData[iCell].X,ELVESData[iCell].Y,ELVESData[iCell].Z,WantedLocationCS);
	    Vector nIn(1.,(pElves-pIn).GetTheta(telCS),(pElves-pIn).GetPhi(telCS), telCS, Vector::kSpherical);
	    nIn.Normalize();
	    const float weight =  ELVESData[iCell].nphotonsnormalized;
	    utl::Photon photonIn(pIn, -nIn, normWavelength, weight);
	    photonIn.SetTime(TimeInterval(photonTime));
	    telSim.AddPhoton(photonIn);
	    photonCounter++;
	}//loop over entry poinits
	cout << "Initialized " << photonCounter << " photons to retrace. " << endl;
      } // End loop over Telescopes    
    } // End loop over Eyes
  }//fStatus
  
  
  
  fStatus = eGeneratePhotons;  
  cout << "=============================================== " << fLoop << endl;
  cout << "=============================================== " << fLoop << endl;
  cout << "=============================================== " << fLoop << endl;

  if (fLoopSelect == 0) {
    fLoop++;
  } else {fLoop = fLoopSelect;}

  return eSuccess;
}

VModule::ResultFlag 
ELVESSimulator::Finish() 
{
  INFO("ELVESSimulator::Finish()");
  return eSuccess;
}


VModule::ResultFlag
ELVESSimulator::ELVESSimDataCreator(){
  
  
  //initialize two different coordinate systems. One at the location where the ELVES was simualted in spherical coordinates. The other at the location where we want the event to be.
  const ReferenceEllipsoid ellipsoid(ReferenceEllipsoid::Get(ReferenceEllipsoid::eWGS84));

  UTMPoint elvesSimulatedLocationUTM(0.0, 0.0, SimulationParameters.AltSource+1400, ellipsoid);//altitude does matter. longitude doesn't matter. 
  const CoordinateSystemPtr SimulatedLocationCS =  fwk::LocalCoordinateSystem::Create(elvesSimulatedLocationUTM.GetPoint());
  
  UTMPoint elvesWantedLocationUTM(fELVESCenterLat, fELVESCenterLon, SimulationParameters.AltSource+1400, ellipsoid);
  const CoordinateSystemPtr WantedLocationCS =  fwk::LocalCoordinateSystem::Create(elvesWantedLocationUTM.GetPoint());
  
  //initializa the eye's CS to be able to calculate the arrival time with respect to each. 
  const fdet::Eye& eye1 = Detector::GetInstance().GetFDetector().GetEye(1);
  const CoordinateSystemPtr& eye1CS = eye1.GetEyeCoordinateSystem();
  const CoordinateSystemPtr& eye1CSLocal = eye1.GetLocalCoordinateSystem();
  const fdet::Eye& eye2 = Detector::GetInstance().GetFDetector().GetEye(2);
  const CoordinateSystemPtr& eye2CS = eye2.GetEyeCoordinateSystem();
  const CoordinateSystemPtr& eye2CSLocal = eye2.GetLocalCoordinateSystem();
  const fdet::Eye& eye3 = Detector::GetInstance().GetFDetector().GetEye(3);
  const CoordinateSystemPtr& eye3CS = eye3.GetEyeCoordinateSystem();
  const CoordinateSystemPtr& eye3CSLocal = eye3.GetLocalCoordinateSystem();
  const fdet::Eye& eye4 = Detector::GetInstance().GetFDetector().GetEye(4);
  const CoordinateSystemPtr& eye4CS = eye4.GetEyeCoordinateSystem();
  const CoordinateSystemPtr& eye4CSLocal = eye4.GetLocalCoordinateSystem();
  
  //bring detector to elves
  //initialize point at origin in each eye cs
  Point peye1eyeCS(0,0,0,eye1CSLocal);
  Point peye2eyeCS(0,0,0,eye2CSLocal);
  Point peye3eyeCS(0,0,0,eye3CSLocal);
  Point peye4eyeCS(0,0,0,eye4CSLocal);

  //initialize points in simulatedCS by using their location in WantedCS
  Point peye1SimulatedCS(peye1eyeCS.GetX(WantedLocationCS),peye1eyeCS.GetY(WantedLocationCS),peye1eyeCS.GetZ(WantedLocationCS),SimulatedLocationCS);
  Point peye2SimulatedCS(peye2eyeCS.GetX(WantedLocationCS),peye2eyeCS.GetY(WantedLocationCS),peye2eyeCS.GetZ(WantedLocationCS),SimulatedLocationCS);
  Point peye3SimulatedCS(peye3eyeCS.GetX(WantedLocationCS),peye3eyeCS.GetY(WantedLocationCS),peye3eyeCS.GetZ(WantedLocationCS),SimulatedLocationCS);
  Point peye4SimulatedCS(peye4eyeCS.GetX(WantedLocationCS),peye4eyeCS.GetY(WantedLocationCS),peye4eyeCS.GetZ(WantedLocationCS),SimulatedLocationCS);

  //initialize CS in new location now, note the locality to Auger, and not to eye itself!! This conserved GetR value in eyeCS for time calculations. 
  eye1CSSim = fwk::LocalCoordinateSystem::Create(peye1SimulatedCS);
  eye2CSSim = fwk::LocalCoordinateSystem::Create(peye2SimulatedCS);
  eye3CSSim = fwk::LocalCoordinateSystem::Create(peye3SimulatedCS);
  eye4CSSim = fwk::LocalCoordinateSystem::Create(peye4SimulatedCS);

  
  //initialize the variables for the input tree. 
  TVector3* position = 0;
  //    double radiusVar, phiVar, time, n, thetaVar;
  float radiusVar, phiVar, time, n, thetaVar;
  
  // return eSuccess;
  
  if(fprodversion == 0){
    fTree->SetBranchAddress("position", &position);
    fTree->SetBranchAddress("nN22P", &n);
    fTree->SetBranchAddress("time", &time);
  }
  if(fprodversion == 1){
    fTree->SetBranchAddress("radius", &radiusVar);
    fTree->SetBranchAddress("theta", &thetaVar);
    fTree->SetBranchAddress("phi", &phiVar);
    fTree->SetBranchAddress("nN22P", &n);
    fTree->SetBranchAddress("time", &time);
  }
  
  TGraph* hAtmo = new TGraph(fNumTreeEntries);
  TGraph* hGeom = new TGraph(fNumTreeEntries);
  bool fdocorrplots = false;

  ELVESDataONE.resize(fNumTreeEntries);
  //first pass to do a time check and select photns within range.
  INFO("First Pass - Read from Tree");
  int previousProgress1 = -1;
  for(int i=0; i<fNumTreeEntries;i+=fPhotonDiscr){
    fTree->GetEntry(i);  
    float thetaVarTMP, phiVarTMP, radiusVarTMP;
    if(fprodversion == 0){
      thetaVarTMP = position->Theta();
      phiVarTMP =  position->Phi();
      radiusVarTMP = position->Mag();
    }
    if(fprodversion == 1){
      thetaVarTMP  = thetaVar;
      phiVarTMP   =  phiVar;
      radiusVarTMP = radiusVar;
    } 

    Point posTMP(ellipsoid.LatitudeLongitudeHeightToPoint((90.*deg)-thetaVarTMP, phiVarTMP, radiusVarTMP-6370*km));

    ELVESSimDataONE ESDtmp;
    ESDtmp.timeEye1 = time+posTMP.GetR(eye1CSSim)/(kSpeedOfLight*1e9);
    ESDtmp.timeEye2 = time+posTMP.GetR(eye2CSSim)/(kSpeedOfLight*1e9);
    ESDtmp.timeEye3 = time+posTMP.GetR(eye3CSSim)/(kSpeedOfLight*1e9);
    ESDtmp.timeEye4 = time+posTMP.GetR(eye4CSSim)/(kSpeedOfLight*1e9);
    ESDtmp.X = posTMP.GetX(SimulatedLocationCS);
    ESDtmp.Y = posTMP.GetY(SimulatedLocationCS);
    ESDtmp.Z = posTMP.GetZ(SimulatedLocationCS);
    //the correction involves the einstein coefficient 2E7, the \Delta t, r, theta, phi, and the arc length is calculated at the 90 km altitude. This is to convert the number density to number of photons in individual cells. Also the output of the simulation integrates over XTimeIntegratedSteps, so a division needs to be done here. 
    ESDtmp.nphotons = (n/SimulationParameters.TimeIntegratedSteps)*2E7*SimulationParameters.SizeStepTime *
      SimulationParameters.SizeStepRadiusHigh * SimulationParameters.SizeStepPhi * 6460E3 * SimulationParameters.SizeStepTheta * 6460E3 ;

    ELVESDataONE[i] = ESDtmp;
    displayProgress(i,fNumTreeEntries,previousProgress1);    
  }
  cout << endl;

  
  std::vector<ELVESSimDataONE>::iterator maxn_it =  std::max_element(ELVESDataONE.begin(), ELVESDataONE.end(), by_nphotonsONE());
  Double_t nMAX = (*maxn_it).nphotons;
  cout << "MAX n: " << nMAX << endl;
  std::vector<ELVESSimDataONE>::iterator minTimeEye1_it =  std::min_element(ELVESDataONE.begin(), ELVESDataONE.end(), by_timeEye1ONE());
  Double_t timeEye1MIN = (*minTimeEye1_it).timeEye1;
  cout << "Min Time E1: " << timeEye1MIN << endl;
  std::vector<ELVESSimDataONE>::iterator minTimeEye2_it =  std::min_element(ELVESDataONE.begin(), ELVESDataONE.end(), by_timeEye2ONE());
  Double_t timeEye2MIN = (*minTimeEye2_it).timeEye2;
  cout << "Min Time E2: " << timeEye2MIN << endl;
  std::vector<ELVESSimDataONE>::iterator minTimeEye3_it =  std::min_element(ELVESDataONE.begin(), ELVESDataONE.end(), by_timeEye3ONE());
  Double_t timeEye3MIN = (*minTimeEye3_it).timeEye3;
  cout << "Min Time E3: " << timeEye3MIN << endl;
  std::vector<ELVESSimDataONE>::iterator minTimeEye4_it =  std::min_element(ELVESDataONE.begin(), ELVESDataONE.end(), by_timeEye4ONE());
  Double_t timeEye4MIN = (*minTimeEye4_it).timeEye4;
  cout << "Min Time E4: " << timeEye4MIN << endl;

  INFO("Second Pass - Calculate Arrival Time at Each Eye and Cut on Time to save RAM");
  int previousProgress2 = -1;
  //second pass: need  to remove the offset of each telescope to 0 the arrival time of the first photon and delete things too far away. 
  int finalDataNEntries = fNumTreeEntries;
  for(int i=0; i<finalDataNEntries;i++){
    ELVESDataONE[i].timeEye1 -= timeEye1MIN;
    ELVESDataONE[i].timeEye2 -= timeEye2MIN;
    ELVESDataONE[i].timeEye3 -= timeEye3MIN;
    ELVESDataONE[i].timeEye4 -= timeEye4MIN;
    float timeTMP = TMath::Max(TMath::Max(ELVESDataONE[i].timeEye1,ELVESDataONE[i].timeEye2),TMath::Max(ELVESDataONE[i].timeEye4,ELVESDataONE[i].timeEye3));
    //only if within window of interest plus a 1 page margin...
    if(timeTMP >= fNPages*1e-4+1e-4){
      std::iter_swap(ELVESDataONE.begin()+i,ELVESDataONE.end()-1);
      ELVESDataONE.pop_back();
      i--;
      finalDataNEntries--;
    }
    displayProgress(i,finalDataNEntries,previousProgress2);
  }
  ELVESDataONE.shrink_to_fit();
  fNPhotonsToCreate = ELVESDataONE.size();
  ELVESData.resize(fNPhotonsToCreate);
  cout << endl <<  "Number of Simulated Source Points within Selection: " << fNPhotonsToCreate << endl;
  

  INFO("Third Pass - Apply Geometric and Atmospheric Factors");
  //Third pass
  int previousProgress3 = -1;//for progress bar
  float attenuationFactor = 10.;
  for(int i=0; i<fNPhotonsToCreate;i++){

    //merge the structs, popback on old and reduce size. 
    ELVESData[i].nphotons = ELVESDataONE[i].nphotons/attenuationFactor;
    ELVESData[i].timeEye1 = ELVESDataONE[i].timeEye1;
    ELVESData[i].timeEye2 = ELVESDataONE[i].timeEye2;
    ELVESData[i].timeEye3 = ELVESDataONE[i].timeEye3;
    ELVESData[i].timeEye4 = ELVESDataONE[i].timeEye4;
    ELVESData[i].X = ELVESDataONE[i].X;
    ELVESData[i].Y = ELVESDataONE[i].Y;
    ELVESData[i].Z = ELVESDataONE[i].Z;
    
    Point posTMP(ELVESData[i].X,ELVESData[i].Y,ELVESData[i].Z,WantedLocationCS);
    
    //this is applying the geometric correction and atmoaspheric correction independently for the individual eyes, wrt to what the elves looks like in their corrd. sys.
    double atmocorrEye1, atmocorrEye2, atmocorrEye3, atmocorrEye4;
    double KYParam_a = 0.50572, KYParam_b = 6.07995*deg, KYParam_c = 1.6364; //Kasten and Young 1989.. careful all defined in degrees, while points are in rad
    double VODTot = fVODTot;
    if(fdoatmocorr){    
      atmocorrEye1=TMath::Exp(-(VODTot)/(TMath::Sin((TMath::Pi()/2.)-posTMP.GetTheta(eye1CS))+KYParam_a*pow(((TMath::Pi()/2.)-posTMP.GetTheta(eye1CS) + KYParam_b)/deg,-KYParam_c)));
      if(fdocorrplots)hAtmo->SetPoint(i,90.-posTMP.GetTheta(eye1CS)/deg, atmocorrEye1);
      atmocorrEye2=TMath::Exp(-(VODTot)/(TMath::Sin((TMath::Pi()/2.)-posTMP.GetTheta(eye2CS))+KYParam_a*pow(((TMath::Pi()/2.)-posTMP.GetTheta(eye2CS) + KYParam_b)/deg,-KYParam_c)));
      atmocorrEye3=TMath::Exp(-(VODTot)/(TMath::Sin((TMath::Pi()/2.)-posTMP.GetTheta(eye3CS))+KYParam_a*pow(((TMath::Pi()/2.)-posTMP.GetTheta(eye3CS) + KYParam_b)/deg,-KYParam_c)));
      atmocorrEye4=TMath::Exp(-(VODTot)/(TMath::Sin((TMath::Pi()/2.)-posTMP.GetTheta(eye4CS))+KYParam_a*pow(((TMath::Pi()/2.)-posTMP.GetTheta(eye4CS) + KYParam_b)/deg,-KYParam_c)));
    }else{
      atmocorrEye1=1; atmocorrEye2=1; atmocorrEye3=1;atmocorrEye4=1;
    }

    double geomcorrEye1, geomcorrEye2, geomcorrEye3,geomcorrEye4;
    if(fdogeomcorr){
      geomcorrEye1 = ((1.8*1.8*kPi)/(4*kPi*posTMP.GetR(eye1CS)*posTMP.GetR(eye1CS)));
      if(fdocorrplots)hGeom->SetPoint(i,posTMP.GetR(eye1CS), geomcorrEye1);
      geomcorrEye2 = ((1.8*1.8*kPi)/(4*kPi*posTMP.GetR(eye2CS)*posTMP.GetR(eye2CS)));
      geomcorrEye3 = ((1.8*1.8*kPi)/(4*kPi*posTMP.GetR(eye3CS)*posTMP.GetR(eye3CS)));
      geomcorrEye4 = ((1.8*1.8*kPi)/(4*kPi*posTMP.GetR(eye4CS)*posTMP.GetR(eye4CS)));      
    }else{
      geomcorrEye1=1; geomcorrEye2=1; geomcorrEye3=1;geomcorrEye4=1;
    }

    ELVESData[i].nphotonsEye1 = ELVESData[i].nphotons*atmocorrEye1*geomcorrEye1;
    ELVESData[i].nphotonsEye2 = ELVESData[i].nphotons*atmocorrEye2*geomcorrEye2;
    ELVESData[i].nphotonsEye3 = ELVESData[i].nphotons*atmocorrEye3*geomcorrEye3;
    ELVESData[i].nphotonsEye4 = ELVESData[i].nphotons*atmocorrEye4*geomcorrEye4;
    displayProgress(i,fNPhotonsToCreate,previousProgress3);
  }
  cout << endl;
  ELVESDataONE.clear();//remove temporary storage of tree data to save RAM. 

  //plot the corrections
  if(fdoatmocorr && fdocorrplots){
    TCanvas *cAtmo = new TCanvas("cAtmo","cAtmo",800,600);
    hAtmo->SetTitle("Atmospheric Attenuation for Individual Grid Cells");
    hAtmo->GetXaxis()->SetTitle("Elevation Angle (Degrees)");
    hAtmo->GetYaxis()->SetTitle("Correction Factor");
    hAtmo->Draw("A*");
    cAtmo->SetLogy();
    cAtmo->SaveAs("outputs/AtmoCorr.png");    
  }
  if(fdogeomcorr && fdocorrplots){
    TCanvas *cGeom = new TCanvas("cGeom","cGeom",800,600);
    hGeom->SetTitle("Geometric Correction for Individual Grid Cells");
    hGeom->GetXaxis()->SetTitle("Distance from Individual Grid Cells (m)");
    hGeom->GetYaxis()->SetTitle("Correction Factor");
    hGeom->Draw("A*");
    cGeom->SetLogy();
    cGeom->SaveAs("outputs/GeomCorr.png");    
  }
  
  /*
    Alright, now the data should be in a format that can be used easily to create the photons in the telescope simulator. A struct was created to contain each of the grid cells simulated with the number density changing through time wrt to each eyes. The next step is to loop through the telescopes and start adding events. ... checks of the data have been done and where presented in Malargue March 2017.  
     */
  
  fInit=true;
  return eSuccess;
}


// This function prints out a progress bar based on how far through the data the program is.                                  
void ELVESSimulator::displayProgress (Int_t currentLoop, Int_t totalLoops, Int_t &previousProgress) {
  float progress = (float)currentLoop / (float)totalLoops;
  if (int(progress*100.0) > previousProgress) {
    previousProgress = int(progress*100.0);
    int barWidth = 70;
    cout << "[";
    int pos = barWidth * progress;
    for (int j = 0; j < barWidth; ++j) {
      if (j < pos) cout << "=";
      else if (j == pos) cout << ">";
      else cout << " ";
    }
    cout << "] " << int(progress * 100.0) << " %\r";
    cout.flush();
  }
  return;
}


