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
#include <TLatex.h>
#include <TH1D.h>
#include <TFile.h>
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

  if(fStatus == eTransformCoordinates){

    INFO("Preparing ELVES Simulated Data...");
    ELVESSimDataCreator();
    Detector& detector = Detector::GetInstance();
    const FDetector& detFD = detector.GetFDetector();
    
    // if (event.HasFEvent())
    //   delete &event.GetFEvent();
    if (!event.HasFEvent())
      event.MakeFEvent();

    if (!event.HasSimShower())
      event.MakeSimShower();
      
    fevt::FEvent& fEvent = event.GetFEvent();
    
    INFO("Creating FEvent");
    
    // set the FD event times
    fSimTime = detector.GetTime();
    event.GetHeader().SetTime(fSimTime);
    //    fEvent.GetHeader().SetTime(fSimTime);
         
    
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
	//	fevt::TelescopeSimData& telSim = telEvent.GetSimData();
	
	//	telSim.SetPhotonsStartTime(fSimTime);
	//	telSim.SetNumberOfPhotonBins(1000);
	
      }
    }  
  }

  if(fStatus == eGeneratePhotons){
    const ReferenceEllipsoid ellipsoid(ReferenceEllipsoid::Get(ReferenceEllipsoid::eWGS84));
    UTMPoint elvesWantedLocationUTM(fELVESCenterLat, fELVESCenterLon, 5.0*km, ellipsoid);
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
      // set the sim header id

      if(eyeId > 4)continue;
      ostringstream idStr;
      idStr << "eye" << eyeId <<"_run"<< "1" << "_event" <<  fLoop;
      event.GetHeader().SetId(idStr.str());
      fEvent.GetHeader().SetId(fLoop);
      //      if(!eyeEvent.HasHeader()) eyeEvent.MakeHeader();
      //      EyeHeader& eEvent = eyeEvent.GetHeader();
      // eEvent.SetRunNumber(1);
      // eEvent.SetEventNumber(fLoop);
      //cout << endl<<endl<<endl<< eyeId<<endl<<endl<<fEvent.GetHeader().GetId() << endl;
      
      for (fevt::Eye::TelescopeIterator iTel = eyeEvent.TelescopesBegin(fevt::ComponentSelector::eInDAQ),
	     end = eyeEvent.TelescopesEnd(fevt::ComponentSelector::eInDAQ);
	   iTel != end; ++iTel) {
	
	const unsigned int telId = iTel->GetId();
	
	fevt::Telescope& telEvent = *iTel;
       	telEvent.SetTracesStartTime(fSimTime+TimeInterval((fLoop-1)*100000));
	//	telEvent.SetTracesStartTime(TimeStamp(0));
	const fdet::Telescope& detTel = detFD.GetTelescope(telEvent);
	
	const double rDia = detTel.GetDiaphragmRadius();
	//	const double diaphragmArea = detTel.GetDiaphragmArea();

	fevt::TelescopeSimData& telSim = telEvent.GetSimData();
	telSim.SetPhotonsStartTime(fSimTime+TimeInterval((fLoop-1)*100000)); 

	//	fevt::TelescopeSimData& telSim = telEvent.GetSimData();
	//	const double binWidth = detTel.GetCamera().GetFADCBinSize();
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
	
	

      //sort data for the given eye. and set the timing to be used below to the one of the right eye. Maybe try to find a faster way to do this. 
	//	int normFactor = 1000;
	if (eyeId == 1){
	  sort(ELVESData.begin(), ELVESData.end(), by_timeEye1());
          for (int iTime = 0; iTime < fNPhotonsToCreate; iTime ++){
	    ELVESData[iTime].time = ELVESData[iTime].timeEye1;
	    //	    ELVESData[iTime].nphotonsnormalized = normFactor*(ELVESData[iTime].nphotonsEye1 - fNPhotonsMINEye1)/(fNPhotonsMAXEye1-fNPhotonsMINEye1);    
	    ELVESData[iTime].nphotonsnormalized = ELVESData[iTime].nphotonsEye1 ;    
	  }
	}
	if (eyeId == 2){
	  sort(ELVESData.begin(), ELVESData.end(), by_timeEye2());
	  for (int iTime = 0; iTime < fNPhotonsToCreate; iTime ++){
	    ELVESData[iTime].time = ELVESData[iTime].timeEye2;
	    //	    ELVESData[iTime].nphotonsnormalized = normFactor*(ELVESData[iTime].nphotonsEye2 - fNPhotonsMINEye2)/(fNPhotonsMAXEye2-fNPhotonsMINEye2);
	    ELVESData[iTime].nphotonsnormalized = ELVESData[iTime].nphotonsEye2 ;    

	  }
	}
	if (eyeId == 3){
	  sort(ELVESData.begin(), ELVESData.end(), by_timeEye3());
	  for (int iTime = 0; iTime < fNPhotonsToCreate; iTime ++) {
	    ELVESData[iTime].time = ELVESData[iTime].timeEye3;
	    ELVESData[iTime].nphotonsnormalized = ELVESData[iTime].nphotonsEye3 ;    
	    //	    ELVESData[iTime].nphotonsnormalized = normFactor*(ELVESData[iTime].nphotonsEye3 - fNPhotonsMINEye3)/(fNPhotonsMAXEye3-fNPhotonsMINEye3);    
	  }
	}
	if (eyeId == 4){
	  sort(ELVESData.begin(), ELVESData.end(), by_timeEye4());
	  for (int iTime = 0; iTime < fNPhotonsToCreate; iTime ++){
	    ELVESData[iTime].time = ELVESData[iTime].timeEye4;
	    //	    ELVESData[iTime].nphotonsnormalized = normFactor*(ELVESData[iTime].nphotonsEye4 - fNPhotonsMINEye4)/(fNPhotonsMAXEye4-fNPhotonsMINEye4);
	    ELVESData[iTime].nphotonsnormalized = ELVESData[iTime].nphotonsEye4 ;    

	  }
	}
	
	// const unsigned int nBins = 1000;
	// const double tracebin = 1e-7;
	
	telSim.SetNumberOfPhotonBins(1000);
	//	telSim.SetNumberOfPhotonBins(700);
	
	telSim.ClearPhotons();
	// if (!telSim.HasPhotonTrace(FdConstants::eFluorDirect))
	//   telSim.MakePhotonTrace(FdConstants::eFluorDirect,1);
	// if (!telSim.HasRayTracedPhotonTrace())
	//   telSim.MakeRayTracedPhotonTrace(nBins, tracebin);
	// telSim.ClearRayTracedPhotonTrace();
	//	TraceI& rayTracedPhotonTrace = telSim.GetRayTracedPhotonTrace();
 
	//	const fdet::Camera& detCamera = detTel.GetCamera();
	//	const double normWavelength = detFD.GetReferenceLambda();
	// const double telTraceBinWidth = detCamera.GetFADCBinSize();
	// const double totalTraceDuration = tracebin * nBins;
	// const unsigned int telTraceNBins = int(0.5 + totalTraceDuration / telTraceBinWidth);
	//telSim.SetNumberOfPhotonBins(telTraceNBins);
	
	
	
	TString hTelName("hTelCS"); hTelName+=telId;
	TString hPINName("hPIN"); hPINName+=telId;
	TString hAINName("hAIN"); hAINName+=telId;
	TH2F* hTelCS = new TH2F(hTelName, hTelName,60,-15,15,60,0,30);
	TH2F* hPIN = new TH2F(hPINName, hPINName,60,-1.5,1.5,60,-1.5,1.5);
	TH2F* hAIN = new TH2F(hAINName,hAINName,60,-100,100,60,-1,10);
	//TH2F* hTelCS = new TH2F(hTelName, hTelName,100,-180,180,40,-20,20);
	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gStyle->SetOptStat(0);
	
	for(int iCell=0; iCell<fNPhotonsToCreate;iCell+=fPhotonDiscr){
	  //check if we are in the window of interest.
	  if(!(ELVESData[iCell].time >= (fLoop-1)*0.0001 && ELVESData[iCell].time < (fLoop)*0.0001)) continue;
	  //	  if(!(ELVESData[iCell].time >= (fLoop-1)*0.00007 && ELVESData[iCell].time < (fLoop)*0.00007)) continue;
	  //if (ELVESData[iCell].nphotonsnormalized < 1.0) continue;

	  //	  double photonNumber = ELVESData[iCell].nphotonsnormalized;
	  for(int iPhoton=0; iPhoton<fNumDiaGridPoints; iPhoton++){

	    // Generating random point on the diaphragm
	    // Uniform phi
	    
	    const double diaPh = RandFlat::shoot(&fRandomEngine->GetEngine(), 0.0, kTwoPi);
	    const double diaR = rDia * sqrt(RandFlat::shoot(&fRandomEngine->GetEngine(), 0.2, 1.0));
	    const double xDia = diaR * cos(diaPh);
	    const double yDia = diaR * sin(diaPh);
	    const Point pIn(xDia, yDia, 0.0, telCS);
	    
	    
	    
	    //	    double photonTime=(ELVESData[iCell].time  - (fLoop-1)*0.00007)*1e9;//need floop to reset time at each event to 0. Also need nanoseconds? and need to shift by 30 us to stay away form pedestals... attempt
	    double photonTime=(ELVESData[iCell].time  - (fLoop-1)*0.0001)*1e9;
	    // photonTime = round(photonTime/100)*100+30000;
	    photonTime = round(photonTime/100)*100;
	    //cout << photonTime << endl;
	    
	    Vector nIn (1.,(ELVESData[iCell].positions-pIn).GetTheta(telCS),(ELVESData[iCell].positions-pIn).GetPhi(telCS), telCS, Vector::kSpherical);
	    
	    nIn.Normalize();
	    const double weight =  ELVESData[iCell].nphotonsnormalized / fNumDiaGridPoints;
	    utl::Photon photonIn(pIn, -nIn, normWavelength, weight);
	    photonIn.SetTime(TimeInterval(photonTime));
	    telSim.AddPhoton(photonIn);
	    
	    double elevationCell = 90.-ELVESData[iCell].positions.GetTheta(eyeCS)/deg;
	    double azimuthCell = (ELVESData[iCell].positions.GetPhi(eyeCS)-telAxis.GetPhi(eyeCS))/deg;
	    hTelCS->Fill(azimuthCell,elevationCell,weight);
	    hPIN->Fill(xDia,yDia);
	    hAIN->Fill((ELVESData[iCell].positions-pIn).GetPhi(telCS)/deg,(ELVESData[iCell].positions-pIn).GetTheta(telCS)/deg,weight);
	  }//lloop over entry poinits
	} // loop over all simulated entries. 
	c1->SetRightMargin(0.16);
	hTelCS->Draw("colz");
	c1->SaveAs((TString)hTelCS->GetName()+".png");
	hPIN->Draw("colz");
	c1->SaveAs((TString)hPIN->GetName()+".png");
	hAIN->Draw("colz");
	c1->SaveAs((TString)hAIN->GetName()+".png");
	
      } // End loop over Telescopes
      
    } // End loop over Eyes
  }
  //Notes: the weights are working
  // however, I would like to find the way to start the traces earlier, to take advantage of full trace
  // also, in CompStudy, i am trying to read out the files properly so that i can overlay traces from both data and simualtion.
  // the header names and the a lot of other stuff still cause the codes to complain.
  // the profiles still don;t seem too right. Maybe it is because we are lloking at the traces with a weird binning.,
  //AH still, it seems the coverage we were getting before, is not the same on the wole FOV. Not enough pixels are covered.
  //CID seem to be the cause of double elves. I would like to look into them more. as well as magnetic effects!
  
  // fInit = true; //initialization actually done at first round now it's time to generate photons
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

    //the correction involves the einstein coefficient 2E7, the \Delta t, r, theta, phi, and the arc length is calculated at the 90 km altitude. This is to convert the number density to number of photons in individual cells. 
    ESDtmp.nphotons = n*2E7*SimulationParameters.SizeStepTime *
      SimulationParameters.SizeStepRadiusHigh * SimulationParameters.SizeStepPhi * 6460E3 * SimulationParameters.SizeStepTheta * 6460E3 ;
    //this is applying the geometric correction independently for the individual eyes. 
    ESDtmp.nphotonsEye1 = ESDtmp.nphotons *((3.6*3.6*kPi)/(4*kPi*posTMP.GetR(eye1CS)*posTMP.GetR(eye1CS)));
    ESDtmp.nphotonsEye2 = ESDtmp.nphotons *((3.6*3.6*kPi)/(4*kPi*posTMP.GetR(eye2CS)*posTMP.GetR(eye2CS)));
    ESDtmp.nphotonsEye3 = ESDtmp.nphotons *((3.6*3.6*kPi)/(4*kPi*posTMP.GetR(eye3CS)*posTMP.GetR(eye3CS)));
    ESDtmp.nphotonsEye4 = ESDtmp.nphotons *((3.6*3.6*kPi)/(4*kPi*posTMP.GetR(eye4CS)*posTMP.GetR(eye4CS)));
    ESDtmp.positions = posTMP;
    ELVESData.push_back(ESDtmp);
    displayProgress(i,fNumTreeEntries,previousProgress);
  }
    
  fNPhotonsToCreate = ELVESData.size();
  cout << endl <<  "Number of Simulated Source Points: " << fNPhotonsToCreate << endl;
  
  //sorting the struct to get the ranges in time for each eye and the num density
  
  sort(ELVESData.begin(), ELVESData.end(), by_nphotons());
  fNPhotonsMIN = ELVESData[0].nphotons;
  fNPhotonsMAX = ELVESData[ELVESData.size()-1].nphotons;
  cout << "nPhotons Range NO GEOM: " << ELVESData[0].nphotons  << " " << ELVESData[ELVESData.size()-1].nphotons << endl;
  sort(ELVESData.begin(), ELVESData.end(), by_nphotonsEye1());
  fNPhotonsMINEye1 = ELVESData[0].nphotonsEye1;
  fNPhotonsMAXEye1 = ELVESData[ELVESData.size()-1].nphotonsEye1;
  cout << "nPhotons Range E1 GEOM: " << ELVESData[0].nphotonsEye1  << " " << ELVESData[ELVESData.size()-1].nphotonsEye1 << endl;
  sort(ELVESData.begin(), ELVESData.end(), by_nphotonsEye2());
  fNPhotonsMINEye2 = ELVESData[0].nphotonsEye2;
  fNPhotonsMAXEye2 = ELVESData[ELVESData.size()-1].nphotonsEye2;
  cout << "nPhotons Range E2 GEOM: " << ELVESData[0].nphotonsEye2  << " " << ELVESData[ELVESData.size()-1].nphotonsEye2 << endl;
  sort(ELVESData.begin(), ELVESData.end(), by_nphotonsEye3());
  fNPhotonsMINEye3 = ELVESData[0].nphotonsEye3;
  fNPhotonsMAXEye3 = ELVESData[ELVESData.size()-1].nphotonsEye3;
  cout << "nPhotons Range E3 GEOM: " << ELVESData[0].nphotonsEye3  << " " << ELVESData[ELVESData.size()-1].nphotonsEye3 << endl;
  sort(ELVESData.begin(), ELVESData.end(), by_nphotonsEye4());
  fNPhotonsMINEye4 = ELVESData[0].nphotonsEye4;
  fNPhotonsMAXEye4 = ELVESData[ELVESData.size()-1].nphotonsEye4;
  cout << "nPhotons Range E4 GEOM: " << ELVESData[0].nphotonsEye4  << " " << ELVESData[ELVESData.size()-1].nphotonsEye4 << endl;
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


VModule::ResultFlag 
ELVESSimulator::Finish() 
{
  INFO("ELVESSimulator::Finish()");
  return eSuccess;
}

