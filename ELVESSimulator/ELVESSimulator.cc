#include "ELVESSimulator.h"
#include <iostream>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TFile.h>
#include <TVector3.h>
#include <vector>


#include <utl/Point.h>
#include <utl/UTMPoint.h>
#include <utl/Branch.h>
#include <utl/Particle.h>
#include <utl/Vector.h>
#include <utl/TimeStamp.h>
#include <utl/Particle.h>
#include <utl/ErrorLogger.h>
#include <utl/ReferenceEllipsoid.h>
#include <utl/RandomEngine.h>
#include <utl/Math.h>
#include <utl/MathConstants.h>
#include <utl/CoordinateSystemPtr.h>
#include <utl/AxialVector.h>
#include <utl/Vector.h>
#include <utl/PhysicalConstants.h>
#include <utl/AugerUnits.h>

#include <fwk/LocalCoordinateSystem.h>
#include <fwk/CoordinateSystemRegistry.h>
#include <fwk/CentralConfig.h>
#include <fwk/SVNGlobalRevision.h>
#include <fwk/RandomEngineRegistry.h>

#include <det/Detector.h>

#include <fdet/FDetector.h>
#include <fdet/Eye.h>

#include <fevt/FEvent.h>
#include <evt/Event.h>
// #include <fevt/Eye.h>
// #include <fevt/TelescopeSimData.h>
// #include <fevt/Telescope.h>
// #include <fevt/PixelSimData.h>
// #include <fevt/Pixel.h>
// #include <fevt/FdConstants.h>

#include <CLHEP/Random/Randomize.h>

using namespace std;
using namespace utl;
using namespace fwk;
using namespace det;
using namespace fdet;
using namespace evt;
using namespace fevt;


namespace ELVESSimulatorNS {

  ELVESSimulator::ELVESSimulator() 
  {
  }
  
  ELVESSimulator::~ELVESSimulator()
  {
  }
  
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
    topB.GetChild("NumberOfDiaphragmEntryPoints").GetData(fNumDiaGridPoints);
    topB.GetChild("PhotonDiscretization").GetData(fPhotonDiscr);
    topB.GetChild("NumberOfTreeEntries").GetData(fNumTreeEntries);
    topB.GetChild("FrameSelection").GetData(fLoopSelect);
    
    
    fIn = new TFile(fELVESInputName.data());
    fTree = (TTree*)fIn->Get(fELVESTreeName.data());
    if (fNumTreeEntries == 0) fNumTreeEntries = fTree->GetEntries();
    //(if set to zero, run over all entries)
    
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
      "Number of Diaphragm Points: " << fNumDiaGridPoints << "\n"
          "     Photon Discretization: " << fPhotonDiscr << "\n"
      "    Number of Tree Entries: " << fNumTreeEntries << "\n"
      "            Frame Selected: " << fLoopSelect << "\n"
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
    
    INFO("ELVESSimulator::Run()");

    //do coordinate transform and time propagauation for each eyes. 
    if(fStatus == eTransformCoordinates){
      ELVESSimDataCreator();
      
      //initialize the event 
      Detector& detector = Detector::GetInstance();

      const fdet::FDetector& detFD = detector.GetFDetector();
      
      if (!event.HasFEvent())
	event.MakeFEvent();
      
      const FEvent& fEvent = event.GetFEvent();
      
      INFO("Creating FEvent");
      
      // set the FD event times
      fSimTime = detector.GetTime();
      fEvent.GetHeader().SetTime(fSimTime);
      //      event.GetHeader().SetTime(fSimTime);
      
      // set the sim header id
      ostringstream idStr;
      idStr << "ELVESSimulation_0" ;
      event.GetHeader().SetId(idStr.str());
      
      for (FDetector::EyeIterator iEye = detFD.EyesBegin();
	   iEye != detFD.EyesEnd() ; ++iEye) {
	const unsigned int eyeId = iEye->GetId();
	fEvent.MakeEye(eyeId, fevt::ComponentSelector::eInDAQ);
	fevt::Eye& eyeEvent = fEvent.GetEye(eyeId, fevt::ComponentSelector::eInDAQ);

	for (fdet::Eye::TelescopeIterator iTel = iEye->TelescopesBegin();
	     iTel != iEye->TelescopesEnd(); ++iTel) {
	  const unsigned int telId = iTel->GetId();
	  eyeEvent.MakeTelescope(telId, fevt::ComponentSelector::eInDAQ);
	  fevt::Telescope& telEvent = eyeEvent.GetTelescope(telId, fevt::ComponentSelector::eInDAQ);
	  telEvent.MakeSimData();
	  fevt::TelescopeSimData& telSim = telEvent.GetSimData();	  
	  telSim.SetNumberOfPhotonBins(700);	  
	}
      }            
    }
    /*    if(fStatus == eGeneratePhotons){
      fevt::FEvent& fEvent = event.GetFEvent();
      Detector& detector = Detector::GetInstance();
      const FDetector& detFD = detector.GetFDetector();
      fEvent.GetHeader().SetTime(fSimTime+TimeInterval((fLoop-1)*100000));
      fEvent.SetId(fLoop);
      for (fevt::FEvent::EyeIterator iEye = fEvent.EyesBegin(fevt::ComponentSelector::eInDAQ),
	     end = fEvent.EyesEnd(fevt::ComponentSelector::eInDAQ);
	   iEye != end; ++iEye) {

	unsigned int eyeId = iEye->GetId();
	fevt::Eye& eyeEvent = *iEye;
	
	for (fevt::Eye::TelescopeIterator iTel = eyeEvent.TelescopesBegin(fevt::ComponentSelector::eInDAQ),
	       end = eyeEvent.TelescopesEnd(fevt::ComponentSelector::eInDAQ);
	     iTel != end; ++iTel) {
	  
	    const unsigned int telId = iTel->GetId();
	  
	  fevt::Telescope& telEvent = *iTel;
	  telEvent.SetTracesStartTime(fSimTime+TimeInterval((fLoop-1)*100000));
	  //       	telEvent.SetTracesStartTime(TimeStamp(0));
	  const fdet::Telescope& detTel = detFD.GetTelescope(telEvent);
	  
	  const double rDia = detTel.GetDiaphragmRadius();
	  const double diaphragmArea = detTel.GetDiaphragmArea();
	  
	  
	  fevt::TelescopeSimData& telSim = telEvent.GetSimData();
	  telSim.SetPhotonsStartTime(fSimTime+TimeInterval((fLoop-1)*100000)); 

	  //	fevt::TelescopeSimData& telSim = telEvent.GetSimData();
	  const double binWidth = detTel.GetCamera().GetFADCBinSize();
	  const double normWavelength = detFD.GetReferenceLambda();
	  const double fRDiaMin = rDia/10. ;//minimum radius on diaphragm, lets say a tenth of the radius
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
	    
	    //quick cut to only do 1 telescope
	    if (telId != 4) continue;
	    
	    info << "Eye: " <<  eyeId << " Tel: " << telId << endl;
	    info << "Tel Axis wrt refCS (Az,El): " << telAxis.GetPhi(referenceCS)/deg << " "  <<   90.-telAxis.GetTheta(referenceCS)/deg <<  endl;
	    info << "FOV: " << detFOV << endl;
	    info << "ELVES Coord wrt refCS (r, Az, El): " << elvesR << "  " << elvesPhi  << " " << elvesTheta << endl;
	    info << "\n\nAdding Photons to TelescopeSimData... ";
	    INFO(info);
	  }else{
	    info << "Eye: " <<  eyeId << " Tel: " << telId << " - ELVES not in field of view";
	    INFO(info);
	    continue;//not in FOV so don't add photons. 
	  }
	  
	  
	  cout << fSimTime.GetGPSSecond() << " " << fSimTime.GetGPSNanoSecond()	<< " "<< TimeInterval((fLoop-1)*100000) << " " << TimeInterval((fLoop-1)*100000).GetSecond() << " "  << TimeInterval((fLoop-1)*100000).GetNanoSecond() << endl;	
	  
	  
	  
	  unsigned int countRTPhotons = 0;
	  double totalWeightPhotons = 0;
	  
	  
	  //	const CoordinateSystemPtr& telCS = detTel.GetTelescopeCoordinateSystem();
	  const Point telescopePos(0.0, 0.0, 0.0, telCS);
	  
	  //THIS NEEDS TO BE REDEFINED
	  // definition of the light trace at the diaphragm
	  // const unsigned int nBins = distanceTrace.GetSize();
	  // const double tracebin = distanceTrace.GetBinning();
	  const unsigned int nBins = 1000;
	  const double tracebin = 1e-7;
	
	  
	  telSim.ClearPhotons();
	  if (!telSim.HasPhotonTrace(FdConstants::eFluorDirect))
	    telSim.MakePhotonTrace(FdConstants::eFluorDirect,1);
	  if (!telSim.HasRayTracedPhotonTrace())
	    telSim.MakeRayTracedPhotonTrace(nBins, tracebin);
	  telSim.ClearRayTracedPhotonTrace();
	  TraceI& rayTracedPhotonTrace = telSim.GetRayTracedPhotonTrace();
	  
	  const fdet::Camera& detCamera = detTel.GetCamera();
	  //	const double normWavelength = detFD.GetReferenceLambda();
	  const double telTraceBinWidth = detCamera.GetFADCBinSize();
	  const double totalTraceDuration = tracebin * nBins;
	  const unsigned int telTraceNBins = int(0.5 + totalTraceDuration / telTraceBinWidth);
	  //	telSim.SetNumberOfPhotonBins(telTraceNBins);
	  
	  cout << " photonTraceSize=" << nBins
	       << " tracebin=" << tracebin
	       << " telTraceBinWidth=" << telTraceBinWidth
	       << " telTraceNBins=" << telTraceNBins << endl;
	  
	  
	  //      for (unsigned int iTrace = 0; iTrace < nBins; ++iTrace) {
	  
	  unsigned int totalNRayTracedPhotonsPerBin = 0;
	  const unsigned int nRayTrace = 30;
	  
	  totalNRayTracedPhotonsPerBin += nRayTrace;
	  countRTPhotons += nRayTrace;
	  
	  //        for (unsigned int iSample = 0; iSample < nRayTrace; ++iSample) {
	  //	for(int iCell=TimeCutIndices[fLoop-1]; iCell<TimeCutIndices[fLoop];iCell+=5){
	  cout << TimeCutIndices[fLoop-1] << " " <<  TimeCutIndices[fLoop]<< " " << TimeCutIndices[fLoop]-TimeCutIndices[fLoop-1]  << endl;
	  
	  for(int iCell=TimeCutIndices[fLoop-1]; iCell<TimeCutIndices[fLoop];iCell+=fPhotonDiscr){
	    for(int iDiaPoint=0; iDiaPoint<fNumDiaGridPoints;iDiaPoint++){
	      // Generating random point on the diaphragm
	      // Uniform phi
	      
	      const double diaPh = RandFlat::shoot(&fRandomEngine->GetEngine(), 0.0, kTwoPi);
	      const double diaR = rDia * sqrt(RandFlat::shoot(&fRandomEngine->GetEngine(), 0.0, 1.0));
	      const double xDia = diaR * cos(diaPh);
	      const double yDia = diaR * sin(diaPh);
	      const Point pIn(xDia, yDia, 0.0, telCS);
	      
	      
	      
	      //const double distancePhoton = distanceBin + tracebin*RandFlat::shoot(&fRandomEngine->GetEngine(), -0.5, 0.5);
	      
	      //	  Point pointOnShower = showerCore + showerAxis * distancePhoton;
	      //	  Vector photonDir = pointOnShower-telescopePos;
	      //	  double timeShift = 0; // time shift with respect to 1D flight path
	      
	      double photonTime=(ELVESData[iCell].time-ELVESData[0].time - (fLoop-1)*0.0001)*1e9;//need floop to reset time at each event to 0. Also need nanoseconds?
	      
	      //	    photonTime = round(photonTime/100)*100;
	      photonTime = round(photonTime/100)*100;
	      
	      Vector nIn (1.,(ELVESData[iCell].positions-pIn).GetTheta(telCS) ,(ELVESData[iCell].positions-pIn).GetPhi(telCS), telCS, Vector::kSpherical);
	      
	      nIn.Normalize();
	      const double projectedDiaphragmArea = diaphragmArea * nIn.GetCosTheta(telCS);
	      //          const double weight = nBunchPhotons * projectedDiaphragmArea;
	      //	    const double weight = ELVESData[iCell].ndensity * 2.0e1 * 500e-9;//TBD	  
	    //	    const double weight = 1000*ELVESData[iCell].ndensitynormalized;//TBD	  
	      const double weight = 100*ELVESData[iCell].ndensity;//TBD	  
	      //	  const double weight = ELVESData[iCell].ndensity * 1e12;//TBD	  
	    //	  const double weight =  1;
	    //	  totalWeightPhotons += weight;
	      utl::Photon photonIn(pIn, -nIn, normWavelength, weight);
	      //	  photonIn.SetTime(TimeInterval(iTrace*tracebin + timeShift));
	      photonIn.SetTime(TimeInterval(photonTime));
	      telSim.AddPhoton(photonIn);
	      
	      // cout << photonTime << " " << weight <<  endl;
	      
	    } // loop over all simulated entries. 
	    
	    //        rayTracedPhotonTrace[iTrace] += totalNRayTracedPhotonsPerBin;
	    //      } // loop over time bin
	  }      //lloop over entry poinits
	  // info << ", generated " << setw(6) << countRTPhotons << " photons"
	  //      << ", total weight=" << totalWeightPhotons;
	  // INFO(info);
	  
	}	
      }
    }
*/    
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
  ELVESSimulator::ELVESSimDataCreator(){


    //initialize two different coordinate systems. One at the location where the ELVES was simualted in spherical coordinates. The other at the location where we want the event to be.
    const ReferenceEllipsoid ellipsoid(ReferenceEllipsoid::Get(ReferenceEllipsoid::eWGS84));
    UTMPoint elvesSimulatedLocationUTM(0.0, 0.0, 5.0*km, ellipsoid);//altitude does matter. This is just for reference
    UTMPoint elvesWantedLocationUTM(fELVESCenterLat, fELVESCenterLon, 5.0*km, ellipsoid);
    const CoordinateSystemPtr WantedLocationCS =  fwk::LocalCoordinateSystem::Create(elvesWantedLocationUTM.GetPoint());
    const CoordinateSystemPtr SimulatedLocationCS =  fwk::LocalCoordinateSystem::Create(elvesSimulatedLocationUTM.GetPoint());

    //initializa the eye's CS to be able to calculate the arrival time with respect to each. 
    const fdet::Eye& eye1 = Detector::GetInstance().GetFDetector().GetEye(1);
    const CoordinateSystemPtr& eye1CS = eye1.GetEyeCoordinateSystem();
    const fdet::Eye& eye2 = Detector::GetInstance().GetFDetector().GetEye(2);
    const CoordinateSystemPtr& eye2CS = eye2.GetEyeCoordinateSystem();
    const fdet::Eye& eye3 = Detector::GetInstance().GetFDetector().GetEye(3);
    const CoordinateSystemPtr& eye3CS = eye3.GetEyeCoordinateSystem();
    const fdet::Eye& eye4 = Detector::GetInstance().GetFDetector().GetEye(4);
    const CoordinateSystemPtr& eye4CS = eye4.GetEyeCoordinateSystem();


    //initialize the variables for the input tree. 
    TVector3* position = 0;
    Double_t time;
    Double_t n;
    
    //TBD: Final structure of Input files
    fTree->SetBranchAddress("position", &position);
    fTree->SetBranchAddress("nN22P", &n);
    fTree->SetBranchAddress("time", &time);
        
    
    int previousProgress = -1;//for progress bar
    for(int i=0; i<fNumTreeEntries;i++){
      fTree->GetEntry(i);    
      
      UTMPoint locationCellUTM(ellipsoid.LatitudeLongitudeHeightToPoint((90.*deg)-position->Theta(), position->Phi(), position->Mag()-6370*km),ellipsoid);
      
      Point posTMP =   Point(locationCellUTM.GetPoint(SimulatedLocationCS).GetX(SimulatedLocationCS),
			     locationCellUTM.GetPoint(SimulatedLocationCS).GetY(SimulatedLocationCS),
			     locationCellUTM.GetPoint(SimulatedLocationCS).GetZ(SimulatedLocationCS),
			     WantedLocationCS);
      
      ELVESSimData ESDtmp;
      ESDtmp.timeEye1 = time+posTMP.GetR(eye1CS)/(kSpeedOfLight*pow(10,9));
      ESDtmp.timeEye2 = time+posTMP.GetR(eye2CS)/(kSpeedOfLight*pow(10,9));
      ESDtmp.timeEye3 = time+posTMP.GetR(eye3CS)/(kSpeedOfLight*pow(10,9));
      ESDtmp.timeEye4 = time+posTMP.GetR(eye4CS)/(kSpeedOfLight*pow(10,9));
      ESDtmp.ndensity = n;
      ESDtmp.positions = posTMP;
      ELVESData.push_back(ESDtmp);
      
      displayProgress(i,fNumTreeEntries,previousProgress);
    }
    
    fNPhotonsToCreate = ELVESData.size();
    cout << endl <<  "Number of Simulated Source Points: " << fNPhotonsToCreate << endl;

    //sorting the struct to get the ranges in time for each eye and the num density
    sort(ELVESData.begin(), ELVESData.end(), by_ndensity());
    cout << "nN22P Range: " << ELVESData[0].ndensity  << " " << ELVESData[ELVESData.size()-1].ndensity << endl;
    sort(ELVESData.begin(), ELVESData.end(), by_timeEye1());
    cout << "Eye1 Time Range: " << ELVESData[0].timeEye1  << " " << ELVESData[ELVESData.size()-1].timeEye1 << endl; Double_t timeEye1MIN = ELVESData[0].timeEye1;
    sort(ELVESData.begin(), ELVESData.end(), by_timeEye2());
    cout << "Eye2 Time Range: " << ELVESData[0].timeEye2  << " " << ELVESData[ELVESData.size()-1].timeEye2 << endl; Double_t timeEye2MIN = ELVESData[0].timeEye2;
    sort(ELVESData.begin(), ELVESData.end(), by_timeEye3());
    cout << "Eye3 Time Range: " << ELVESData[0].timeEye3  << " " << ELVESData[ELVESData.size()-1].timeEye3 << endl; Double_t timeEye3MIN = ELVESData[0].timeEye3;
    sort(ELVESData.begin(), ELVESData.end(), by_timeEye4());
    cout << "Eye4 Time Range: " << ELVESData[0].timeEye4  << " " << ELVESData[ELVESData.size()-1].timeEye4 << endl; Double_t timeEye4MIN = ELVESData[0].timeEye4;
    //don't forget that whenever an eye is reached, the sorting of the struct needs to be done again. right now we are sorted wrt eye4. 

    //need  to remove the offset of each telescope to 0 the arrival time of the first photon. 
    for(int i=0; i<fNPhotonsToCreate;i++){
      ELVESData[i].timeEye1 -= timeEye1MIN;
      ELVESData[i].timeEye2 -= timeEye2MIN;
      ELVESData[i].timeEye3 -= timeEye3MIN;
      ELVESData[i].timeEye4 -= timeEye4MIN;
    }

    /*
      Alright, now the data should be in a format that can be used easily to create the photons in the telescope simulator. A struct was created to contain each of the grid cells simulated with the number density changing through time wrt to each eyes. The next step is to loop through the telescopes and start adding events. ... checks of the data have been done and where presented in Malargue March 2017.  
     */
    
    fInit=true;
    return eSuccess;
  }

  VModule::ResultFlag 
  ELVESSimulator::Finish() 
  {
    return eSuccess;
  }
  
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

}

// For special applications.
void AugerOfflineUser()
{
}
