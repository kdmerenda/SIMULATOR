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
  topB.GetChild("NumberOfDiaphragmEntryPoints").GetData(fNumDiaGridPoints);
  topB.GetChild("PhotonDiscretization").GetData(fPhotonDiscr);
  topB.GetChild("NumberOfTreeEntries").GetData(fNumTreeEntries);
  topB.GetChild("FrameSelection").GetData(fLoopSelect);
  topB.GetChild("PreCheck").GetData(fdoprecheck);
  topB.GetChild("ProdVersion").GetData(fprodversion);
  topB.GetChild("GeometricCorrection").GetData(fdogeomcorr);
  topB.GetChild("AtmosphericCorrection").GetData(fdoatmocorr);


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
          "               Do PreCheck: " << fdoprecheck << "\n"
          "        Production Version: " << fprodversion << "\n"
          " Do Atmospheric Correction: " << fdoatmocorr << "\n"
          "   Do Geometric Correction: " << fdogeomcorr << "\n"
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
    //    fEvent.GetHeader().SetTime(fSimTime);
         
    
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
      cout << event.GetHeader().GetId() << endl;
      event.GetHeader().SetId(idStr.str());
      fEvent.GetHeader().SetId(fLoop);
       if(!eyeEvent.HasHeader()) eyeEvent.MakeHeader();
      EyeHeader& eEvent = eyeEvent.GetHeader();
      eEvent.SetRunNumber(1);
      eEvent.SetEventNumber(fLoop);
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
	  //	  if (telId != 1) continue;
	  
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
	int normFactor = 4000;
	if (eyeId == 1){
	  sort(ELVESData.begin(), ELVESData.end(), by_timeEye1());
          for (int iTime = 0; iTime < fNPhotonsToCreate; iTime ++){
	    ELVESData[iTime].time = ELVESData[iTime].timeEye1;
	    //ELVESData[iTime].nphotonsnormalized = normFactor*(ELVESData[iTime].nphotonsEye1 - fNPhotonsMINEye1)/(fNPhotonsMAXEye1-fNPhotonsMINEye1);    
	    ELVESData[iTime].nphotonsnormalized = ELVESData[iTime].nphotonsEye1 ;    
	  }
	}
	if (eyeId == 2){
	  sort(ELVESData.begin(), ELVESData.end(), by_timeEye2());
	  for (int iTime = 0; iTime < fNPhotonsToCreate; iTime ++){
	    ELVESData[iTime].time = ELVESData[iTime].timeEye2;
	    //ELVESData[iTime].nphotonsnormalized = normFactor*(ELVESData[iTime].nphotonsEye2 - fNPhotonsMINEye2)/(fNPhotonsMAXEye2-fNPhotonsMINEye2);
	    ELVESData[iTime].nphotonsnormalized = ELVESData[iTime].nphotonsEye2 ;    

	  }
	}
	if (eyeId == 3){
	  sort(ELVESData.begin(), ELVESData.end(), by_timeEye3());
	  for (int iTime = 0; iTime < fNPhotonsToCreate; iTime ++) {
	    ELVESData[iTime].time = ELVESData[iTime].timeEye3;
	    ELVESData[iTime].nphotonsnormalized = ELVESData[iTime].nphotonsEye3 ;    
	    //ELVESData[iTime].nphotonsnormalized = normFactor*(ELVESData[iTime].nphotonsEye3 - fNPhotonsMINEye3)/(fNPhotonsMAXEye3-fNPhotonsMINEye3);    
	  }
	}
	if (eyeId == 4){
	  sort(ELVESData.begin(), ELVESData.end(), by_timeEye4());
	  for (int iTime = 0; iTime < fNPhotonsToCreate; iTime ++){
	    ELVESData[iTime].time = ELVESData[iTime].timeEye4;
	    //ELVESData[iTime].nphotonsnormalized = normFactor*(ELVESData[iTime].nphotonsEye4 - fNPhotonsMINEye4)/(fNPhotonsMAXEye4-fNPhotonsMINEye4);
	    ELVESData[iTime].nphotonsnormalized = ELVESData[iTime].nphotonsEye4 ;    
	    //try w/o normalization but by removing any values below1. putting 1 as min caused values dipping below 0... not good, prob not correct thing to do. 
	  }
	}
	
	// const unsigned int nBins = 1000;
	// const double tracebin = 1e-7;
	
	//	telSim.SetNumberOfPhotonBins(8000);
	//	telSim.SetNumberOfPhotonBins(8000);
	//	telSim.SetNumberOfPhotonBins(700);
	
	 telSim.ClearPhotons();
	// if (!telSim.HasPhotonTrace(FdConstants::eFluorDirect))
	//   telSim.MakePhotonTrace(FdConstants::eFluorDirect,1);
	// if (!telSim.HasRayTracedPhotonTrace())
	//   telSim.MakeRayTracedPhotonTrace(nBins, tracebin);
	// telSim.ClearRayTracedPhotonTrace();	
	

	//ROOT OUTPUT TO CHECK THE SIGNAL BEFORE THE PHOTON CREATION
       	if (fdoprecheck){
	  
	  //mitchell command to make TCutG of hexagonal pixel. 
	  makePixels(eyeId,telId);
	  //return a tree on which to apply the graphical cut 
	  fSimDataTree = makeTelDataTree(telId, telAxis, eyeId, eyeCS);
	  fSimDataTree->Print();

	  
	  //init canvas for visual of data in tel cs with grid
	  TCanvas* c1 = new TCanvas("c1","c1",600,600);
	  gStyle->SetOptStat(0);
	  //check pixel grid overlay
	  TString hTelName("hTelCS"); hTelName+=eyeId; hTelName+=telId;
	 
	  TH2F* hTelCS = new TH2F(hTelName, hTelName,60,-17,17,60,0,32);
	  c1->SetRightMargin(0.16);
	  fSimDataTree->Project(hTelName,"elevation:azimuth","nphotons*(time<72e-6)");
	  hTelCS->Draw("colz");
	  for(int ipix = 0; ipix < pixelCuts.size(); ipix++){
	    if(detTel.GetPixel(ipix+1).GetColumn()==10) pixelCuts[ipix]->SetLineWidth(2);
	    if(detTel.GetPixel(ipix+1).GetRow()==10) pixelCuts[ipix]->SetLineWidth(2);
	    pixelCuts[ipix]->Draw("L");	    
	  }
	  c1->SaveAs((TString)hTelCS->GetName()+".png");
	 
	  //20columns 22rows. 
	  //output traces for column
	  MakeTraces("Column",10,72e-6,1e-6,detTel, telId, eyeId);
	  //output traces for row 
	  //	  MakeTraces("Row",10,300e-6,0.5e-6,detTel, telId, eyeId);
	  
	}//doprecheck
	

	//This is the creation of photons
	for(int iCell=0; iCell<fNPhotonsToCreate;iCell+=fPhotonDiscr){
	  
	  //check if we are in the window of interest.

	  if(!(ELVESData[iCell].time >= (fLoop-1)*0.0001 && ELVESData[iCell].time < (fLoop)*0.0001)) continue;
	  //	 if(fLoop == 1 && !(ELVESData[iCell].time >= 0. && ELVESData[iCell].time < 0.00072)) continue;

	 if (ELVESData[iCell].nphotonsnormalized < 1.0) continue;//is this fair?

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
	    double photonTime=(ELVESData[iCell].time  - (fLoop-1)*100e-6)*1e9;
	    //photonTime = round(photonTime/100)*100;
	    //photonTime = round(photonTime/1000)*1000;//should the binning be done here or in post processing?
	    //cout << photonTime << endl;
	    
	    Vector nIn (1.,(ELVESData[iCell].positions-pIn).GetTheta(telCS),(ELVESData[iCell].positions-pIn).GetPhi(telCS), telCS, Vector::kSpherical);
	    
	    nIn.Normalize();
	    const double weight =  ELVESData[iCell].nphotonsnormalized / fNumDiaGridPoints;
	    utl::Photon photonIn(pIn, -nIn, normWavelength, weight);
	    photonIn.SetTime(TimeInterval(photonTime));
	    telSim.AddPhoton(photonIn);
	    
	  }//loop over entry poinits
	} // loop over all simulated entries. 
      } // End loop over Telescopes    
    } // End loop over Eyes
  }//fStatus
  
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
ELVESSimulator::Finish() 
{
  INFO("ELVESSimulator::Finish()");
  return eSuccess;
}

  
VModule::ResultFlag 
ELVESSimulator::MakeTraces(TString rctag,int rcnum,double maxtime,double timebinsize, const fdet::Telescope& detTel,int telId, int eyeId){


  int numColors = 11;
  Color_t colorList[numColors] = {kPink,kRed,kOrange,kSpring+8,kTeal+8,kTeal,kCyan,kAzure-4,kBlue,kViolet+10,kBlack};//watch out, 9 is hardcoded in this function
  
  TString fNameTraces("PreCheckTraces");fNameTraces+=eyeId;fNameTraces+=telId;fNameTraces+=rctag;fNameTraces+=rcnum;
  
  
  TFile fTraces(fNameTraces + ".root","recreate");
  int hcount = 0;
  TCanvas* cTraces = new TCanvas("cTraces","cTraces",800,600);
  TGaxis::SetMaxDigits(4); // force traces to use scientific notation
  TLegend* lc = new TLegend(.15,.15,.25,.85);
  //lc->SetNColumns(numColors);

  TString sTitle("");sTitle+= rctag;sTitle+=" ";sTitle+= rcnum;
  if (rctag == "Row") sTitle+=";time (s);photon count / "; sTitle+=round(timebinsize*1e7)*1e2; sTitle+=" ns";
  if (rctag == "Column") sTitle+=";Time Bins;photon count / 1 #mus";
  THStack* sc = new THStack("sc",sTitle);
  
  for (int iPixel = 0; iPixel < pixelCuts.size(); iPixel++) {
    TString hName("h"); hName+=fLoop;hName+="_m";hName+=telId;hName+="_r";hName+=detTel.GetPixel(iPixel+1).GetRow();hName+="_c";hName+=detTel.GetPixel(iPixel+1).GetColumn();
    // TString hTraceName("hTrace");
    // hTraceName+=eyeId;hTraceName+=telId;hTraceName+=iPixel;
    TH1F* hTrace;
    if (rctag == "Column") hTrace = new TH1F(hName,hName, 100, 0, 1000);//check if time is converted to bin properly
    if (rctag == "Row") hTrace  = new TH1F(hName,hName, maxtime/timebinsize, 0, maxtime);
    TString bounds = "nphotons*";
    bounds += "(";
    bounds += pixelCuts[iPixel]->GetName();//check if the cuts are defined wrt proper pixel. 
    bounds += "&&time<";
    bounds += maxtime;
    bounds +=")";
    
    
    if(rctag == "Column" && detTel.GetPixel(iPixel+1).GetColumn()==rcnum) fSimDataTree->Project(hName,"TMath::Floor(time*1e7)+280",bounds);//transform to time bin + shift by 280
    
    if(rctag == "Row" && detTel.GetPixel(iPixel+1).GetRow()==rcnum) fSimDataTree->Project(hName,"time",bounds);
    
    if(hTrace->GetEntries()==0){
      delete hTrace;
    }else{
      hTrace->SetLineColor(colorList[hcount]);
      hTrace->SetLineWidth(2);
      
      TString lcName("");
      if(rctag == "Column")lcName+=detTel.GetPixel(iPixel+1).GetRow();
      if(rctag == "Row") lcName+=detTel.GetPixel(iPixel+1).GetColumn();
      lc->AddEntry(hTrace,lcName,"l");
      sc->Add(hTrace);
      if(hcount==numColors)break;
      hcount++;
    }
  }
  
  
  sc->Draw("nostack");
  lc->Draw();
  cTraces->SaveAs(fNameTraces+".png");   	  

  //  delete DataTree;
  fTraces.Write();
  fTraces.Close();
  return eSuccess;
}

TTree* ELVESSimulator::makeTelDataTree(int telId, Vector telAxis, int eyeId,  CoordinateSystemPtr  eyeCS ){
  TString tName("t"); tName+=eyeId; tName+=telId; tName+=".root";
  TTree* tTMP = new TTree(tName,tName);
  double tTime, tElevation, tAzimuth, tNPhotons;
  TBranch *bTime = tTMP->Branch("time",&tTime);
  TBranch *bElevation = tTMP->Branch("elevation",&tElevation);
  TBranch *bAzimuth = tTMP->Branch("azimuth",&tAzimuth);
  TBranch *bNPhotons = tTMP->Branch("nphotons",&tNPhotons);
  
  for(int iCell=0; iCell<fNPhotonsToCreate;iCell+=fPhotonDiscr){
    tElevation = 90.-ELVESData[iCell].positions.GetTheta(eyeCS)/deg;
    tAzimuth = (ELVESData[iCell].positions.GetPhi(eyeCS)-telAxis.GetPhi(eyeCS))/deg;
    tTime = ELVESData[iCell].time;
    tNPhotons = ELVESData[iCell].nphotonsnormalized;
    tTMP->Fill();
  }//entries

  return tTMP;
}

VModule::ResultFlag
ELVESSimulator::ELVESSimDataCreator(){
  
  
  //initialize two different coordinate systems. One at the location where the ELVES was simualted in spherical coordinates. The other at the location where we want the event to be.
  const ReferenceEllipsoid ellipsoid(ReferenceEllipsoid::Get(ReferenceEllipsoid::eWGS84));

  UTMPoint elvesSimulatedLocationUTM(0.0, 0.0, SimulationParameters.AltSource, ellipsoid);//altitude does matter. This is just for reference
  UTMPoint elvesWantedLocationUTM(fELVESCenterLat, fELVESCenterLon, SimulationParameters.AltSource, ellipsoid);
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
  //    double radiusVar, phiVar, time, n, thetaVar;
    float radiusVar, phiVar, time, n, thetaVar;

  
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
  
  int previousProgress = -1;//for progress bar
  TGraph* hAtmo = new TGraph(fNumTreeEntries);
  TGraph* hGeom = new TGraph(fNumTreeEntries);

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
    UTMPoint locationCellUTM(ellipsoid.LatitudeLongitudeHeightToPoint((90.*deg)-thetaVarTMP, phiVarTMP, radiusVarTMP-6370*km),ellipsoid);
    
    Point posTMP =   Point(locationCellUTM.GetPoint(SimulatedLocationCS).GetX(SimulatedLocationCS),
			   locationCellUTM.GetPoint(SimulatedLocationCS).GetY(SimulatedLocationCS),
			   locationCellUTM.GetPoint(SimulatedLocationCS).GetZ(SimulatedLocationCS),
			   WantedLocationCS);
    
    ELVESSimData ESDtmp;
    ESDtmp.timeEye1 = time+posTMP.GetR(eye1CS)/(kSpeedOfLight*pow(10,9));
    ESDtmp.timeEye2 = time+posTMP.GetR(eye2CS)/(kSpeedOfLight*pow(10,9));
    ESDtmp.timeEye3 = time+posTMP.GetR(eye3CS)/(kSpeedOfLight*pow(10,9));
    ESDtmp.timeEye4 = time+posTMP.GetR(eye4CS)/(kSpeedOfLight*pow(10,9));

    //the correction involves the einstein coefficient 2E7, the \Delta t, r, theta, phi, and the arc length is calculated at the 90 km altitude. This is to convert the number density to number of photons in individual cells. Also the output of the simulation integrates over XTimeIntegratedSteps, so a division needs to be done here. 
    ESDtmp.nphotons = (n/SimulationParameters.TimeIntegratedSteps)*2E7*SimulationParameters.SizeStepTime *
      SimulationParameters.SizeStepRadiusHigh * SimulationParameters.SizeStepPhi * 6460E3 * SimulationParameters.SizeStepTheta * 6460E3 ;

    //this is applying the geometric correction and atmoaspheric correction independently for the individual eyes, wrt to what the elves looks like in their corrd. sys. ATMO TBD!!! 

    double atmocorrEye1, atmocorrEye2, atmocorrEye3, atmocorrEye4;
    double KYParam_a = 0.50572, KYParam_b = 6.07995*deg, KYParam_c = 1.6364; //Kasten and Young 1989.. careful all defined in degrees, while points are in rad
    double VODTot = 0.6;
    if(fdoatmocorr){    
      atmocorrEye1=TMath::Exp(-(VODTot)/(TMath::Sin((TMath::Pi()/2.)-posTMP.GetTheta(eye1CS))+KYParam_a*pow(((TMath::Pi()/2.)-posTMP.GetTheta(eye1CS) + KYParam_b)/deg,-KYParam_c)));
      hAtmo->SetPoint(i,90.-posTMP.GetTheta(eye1CS)/deg, atmocorrEye1);

      atmocorrEye2=TMath::Exp(-(VODTot)/(TMath::Sin((TMath::Pi()/2.)-posTMP.GetTheta(eye2CS))+KYParam_a*pow(((TMath::Pi()/2.)-posTMP.GetTheta(eye2CS) + KYParam_b)/deg,-KYParam_c)));
      atmocorrEye3=TMath::Exp(-(VODTot)/(TMath::Sin((TMath::Pi()/2.)-posTMP.GetTheta(eye3CS))+KYParam_a*pow(((TMath::Pi()/2.)-posTMP.GetTheta(eye3CS) + KYParam_b)/deg,-KYParam_c)));
      atmocorrEye4=TMath::Exp(-(VODTot)/(TMath::Sin((TMath::Pi()/2.)-posTMP.GetTheta(eye3CS))+KYParam_a*pow(((TMath::Pi()/2.)-posTMP.GetTheta(eye4CS) + KYParam_b)/deg,-KYParam_c)));
    }else{
      atmocorrEye1=1; atmocorrEye2=1; atmocorrEye3=1;atmocorrEye4=1;
    }

    double geomcorrEye1, geomcorrEye2, geomcorrEye3,geomcorrEye4;
    if(fdogeomcorr){
      geomcorrEye1 = ((1.8*1.8*kPi)/(4*kPi*posTMP.GetR(eye1CS)*posTMP.GetR(eye1CS)));
      hGeom->SetPoint(i,posTMP.GetR(eye1CS), geomcorrEye1);
      geomcorrEye2 = ((1.8*1.8*kPi)/(4*kPi*posTMP.GetR(eye2CS)*posTMP.GetR(eye2CS)));
      geomcorrEye3 = ((1.8*1.8*kPi)/(4*kPi*posTMP.GetR(eye3CS)*posTMP.GetR(eye3CS)));
      geomcorrEye4 = ((1.8*1.8*kPi)/(4*kPi*posTMP.GetR(eye4CS)*posTMP.GetR(eye4CS)));      
    }else{
      geomcorrEye1=1; geomcorrEye2=1; geomcorrEye3=1;geomcorrEye4=1;
    }
    ESDtmp.nphotonsEye1 = ESDtmp.nphotons*atmocorrEye1*geomcorrEye1;
    ESDtmp.nphotonsEye2 = ESDtmp.nphotons*atmocorrEye2*geomcorrEye2;
    ESDtmp.nphotonsEye3 = ESDtmp.nphotons*atmocorrEye3*geomcorrEye3;
    ESDtmp.nphotonsEye4 = ESDtmp.nphotons*atmocorrEye4*geomcorrEye4;
    ESDtmp.positions = posTMP;//get it for each eye at time of photon creation
    ELVESData.push_back(ESDtmp);
    displayProgress(i,fNumTreeEntries,previousProgress);
  }
  cout << endl;
  if(fdoatmocorr){
    TCanvas *cAtmo = new TCanvas("cAtmo","cAtmo",800,600);
    hAtmo->SetTitle("Atmospheric Attenuation for Individual Grid Cells");
    hAtmo->GetXaxis()->SetTitle("Elevation Angle (Degrees)");
    hAtmo->GetYaxis()->SetTitle("Correction Factor");
    hAtmo->Draw("A*");
    cAtmo->SetLogy();
    cAtmo->SaveAs("AtmoCorr.png");    
  }
  if(fdogeomcorr){
    TCanvas *cGeom = new TCanvas("cGeom","cGeom",800,600);
    hGeom->SetTitle("Geometric Correction for Individual Grid Cells");
    hGeom->GetXaxis()->SetTitle("Distance from Individual Grid Cells (m)");
    hGeom->GetYaxis()->SetTitle("Correction Factor");
    hGeom->Draw("A*");
    cGeom->SetLogy();
    cGeom->SaveAs("GeomCorr.png");    
  }

  fNPhotonsToCreate = ELVESData.size();
  cout << endl <<  "Number of Simulated Source Points: " << fNPhotonsToCreate << endl;
  
  //sorting the struct to get the ranges in time for each eye and the num density
  
  sort(ELVESData.begin(), ELVESData.end(), by_nphotons());
  fNPhotonsMIN = ELVESData[0].nphotons;
  fNPhotonsMAX = ELVESData[ELVESData.size()-1].nphotons;
  cout << "nPhotons Range at Source Point: " << ELVESData[0].nphotons  << " " << ELVESData[ELVESData.size()-1].nphotons << endl;
  sort(ELVESData.begin(), ELVESData.end(), by_nphotonsEye1());
  fNPhotonsMINEye1 = ELVESData[0].nphotonsEye1;
  //fNPhotonsMINEye1 = 1;
  fNPhotonsMAXEye1 = ELVESData[ELVESData.size()-1].nphotonsEye1;
  cout << "nPhotons Range E1: " << ELVESData[0].nphotonsEye1  << " " << ELVESData[ELVESData.size()-1].nphotonsEye1 << endl;
  sort(ELVESData.begin(), ELVESData.end(), by_nphotonsEye2());
  fNPhotonsMINEye2 = ELVESData[0].nphotonsEye2;
  //fNPhotonsMINEye2 = 1;
  fNPhotonsMAXEye2 = ELVESData[ELVESData.size()-1].nphotonsEye2;
  cout << "nPhotons Range E2: " << ELVESData[0].nphotonsEye2  << " " << ELVESData[ELVESData.size()-1].nphotonsEye2 << endl;
  sort(ELVESData.begin(), ELVESData.end(), by_nphotonsEye3());
  fNPhotonsMINEye3 = ELVESData[0].nphotonsEye3;
  //fNPhotonsMINEye3 = 1;
  fNPhotonsMAXEye3 = ELVESData[ELVESData.size()-1].nphotonsEye3;
  cout << "nPhotons Range E3: " << ELVESData[0].nphotonsEye3  << " " << ELVESData[ELVESData.size()-1].nphotonsEye3 << endl;
  sort(ELVESData.begin(), ELVESData.end(), by_nphotonsEye4());
  fNPhotonsMINEye4 = ELVESData[0].nphotonsEye4;
  //  fNPhotonsMINEye4 = 1;
  fNPhotonsMAXEye4 = ELVESData[ELVESData.size()-1].nphotonsEye4;
  cout << "nPhotons Range E4: " << ELVESData[0].nphotonsEye4  << " " << ELVESData[ELVESData.size()-1].nphotonsEye4 << endl;
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
ELVESSimulator::makePixels(int eyeId, int telId) {
  const fdet::Eye& eye = Detector::GetInstance().GetFDetector().GetEye(eyeId);
  const fdet::Telescope& tel = eye.GetTelescope(telId);
  const CoordinateSystemPtr& eyeCS = eye.GetEyeCoordinateSystem();
  Double_t telPhi = tel.GetAxis().GetPhi(eyeCS);
  
  pixelCuts.clear();
  TCutG* cut;
  for (unsigned int channelId = 1; channelId <= tel.GetLastPixelId(); ++channelId) {
    const fdet::Channel& detChannel = tel.GetChannel(channelId);
    unsigned int pixelId = detChannel.GetPixelId();
    const fdet::Pixel& pixel = tel.GetPixel(pixelId);
    Vector pixAxis = pixel.GetDirection();    
    Double_t centerX = pixAxis.GetPhi(eyeCS)/deg - telPhi/deg;
    Double_t centerY = 90.-pixAxis.GetTheta(eyeCS)/deg;
    Double_t r = 1.5/1.8;

    // create TCutG hexagon of radius r around (centerX, centerY)
    TString cutname("cutname");
    cutname += pixelId;
    cut  = new TCutG(cutname, 7);
    double startAngle = kPi/6;
    for (int point = 0; point < 7; point++) {
      cut->SetPoint(point, centerX+r*cos(startAngle+point*kPi/3), centerY+r*sin(startAngle+point*kPi/3));
    }
    cut->SetVarX("azimuth");
    cut->SetVarY("elevation");
    pixelCuts.push_back(cut);
  }
  
  return eSuccess;
}
