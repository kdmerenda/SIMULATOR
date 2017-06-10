//ELVES Simulation Converter from EMP Sim to Offline Raw File.

#include "ELVESPhotonGenerator.h"
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

#include <TLegend.h>
#include <TMap.h>
#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TFrame.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <stdio.h>
#include <TCutG.h>
#include <TColor.h>
#include <TStyle.h>
#include <TGaxis.h>

using namespace std;
using namespace utl;
using namespace fwk;
using namespace evt;
using namespace fevt;
using namespace fdet;
using namespace det;

VModule::ResultFlag
ELVESPhotonGenerator::Init()
{

  INFO("ELVESPhotonGenerator::Init()");

  old_file_name = "/home/mwaibel/workoffline/ELVESSim/PA030517.root";

  //  TFile *fIn = new TFile("/home/kmerenda/simelve/runs/forOffline/Auger3D_1e06/saved/PA110216.root");
  TFile *fIn = new TFile(old_file_name);
  fTree = (TTree*)fIn->Get("ELVES");
  //  fTreeParameters = (TTree*)fIn->Get("ELVESParameters");

  return eSuccess;
}

VModule::ResultFlag
ELVESPhotonGenerator::Run(evt::Event& event)
{

  /*
    =======================================
    BEGIN coordinate transformation code
    =======================================
  */

//  createNewTFile("/home/mwaibel/workoffline/ELVESSim/PA121016.root", 1, 1);

  const ReferenceEllipsoid ellipsoid(ReferenceEllipsoid::Get(ReferenceEllipsoid::eWGS84));
  const CoordinateSystemPtr referenceCS = Detector::GetInstance().GetReferenceCoordinateSystem();
  UTMPoint elvesSimulatedLocationUTM(0.0, 0.0, 5.0*km, ellipsoid);
  UTMPoint elvesWantedLocationUTM(sourceLat, sourceLon, 5.0*km, ellipsoid);
  const CoordinateSystemPtr WantedLocationCS = fwk::LocalCoordinateSystem::Create(elvesWantedLocationUTM.GetPoint());
  const CoordinateSystemPtr SimulatedLocationCS = fwk::LocalCoordinateSystem::Create(elvesSimulatedLocationUTM.GetPoint());

  TVector3* position = 0;
  Double_t time;
  Double_t n;
  fTree->SetBranchAddress("position", &position);
  fTree->SetBranchAddress("nN22P", &n);
  fTree->SetBranchAddress("time", &time);
  const int nEntries = fTree->GetEntries();
  vector<utl::Point> PositionPoints(nEntries);
  vector<Double_t> nvec(nEntries);
  vector<Double_t> timevec(nEntries);

//  Double_t maxTime = 0;
//  Double_t minTime = 10;
/*
  Int_t previousProgress = -1;
  cout << "Progress through first loop:" << endl;
  int loopTimes = nEntries; // Say how much of the data to go through
  for(int i=0; i<loopTimes ;i++){
    fTree->GetEntry(i);
    nvec[i]=n;
    timevec[i]=time;
    UTMPoint locationCellUTM(position->Theta()-(90.*deg)+sourceLat, position->Phi()+sourceLon, position->Mag()-6370*km,ellipsoid); // make a UTMPoint of the data point
    PositionPoints[i] = locationCellUTM.GetPoint(referenceCS);

    displayProgress(i, loopTimes, previousProgress);
  }
  cout << endl;
*/

  /*
    ==================================
    END: coordinate transformation code
    ==================================
  */


  /*
    =======================================
    BEGIN: event creation
    =======================================
  */
  Detector& detector = Detector::GetInstance();
  const FDetector& detFD = detector.GetFDetector();
  
  if (!event.HasFEvent())
    event.MakeFEvent();

  fevt::FEvent& fEvent = event.GetFEvent();
     
  const TimeStamp& simTime = detector.GetTime();
  
  // set the FD event times
  fEvent.GetHeader().SetTime(simTime);
  event.GetHeader().SetTime(simTime);
  
  for (FDetector::EyeIterator iEye = detFD.EyesBegin();
       iEye != detFD.EyesEnd() ; ++iEye) {
    const unsigned int eyeId = iEye->GetId();
    //   cout << eyeId << endl;
    
    fEvent.MakeEye(eyeId, fevt::ComponentSelector::eInDAQ);
    fevt::Eye& eyeEvent = fEvent.GetEye(eyeId, fevt::ComponentSelector::eInDAQ);
    for (fdet::Eye::TelescopeIterator iTel = iEye->TelescopesBegin();
      iTel != iEye->TelescopesEnd(); ++iTel) {
      const unsigned int telId = iTel->GetId();
      // cout << telId << " I'm in the loop "  << endl;
      
      if (eyeEvent.HasTelescope(telId, fevt::ComponentSelector::eInDAQ)) {
        ERROR("Telescope already there");
      }
      eyeEvent.MakeTelescope(telId, fevt::ComponentSelector::eInDAQ);
      fevt::Telescope& telEvent = eyeEvent.GetTelescope(telId, fevt::ComponentSelector::eInDAQ);
      telEvent.MakeSimData();
      //      cout << telId << endl; displays all telescopes numbers. 
      fevt::TelescopeSimData& telSim = telEvent.GetSimData();
      telSim.SetPhotonsStartTime(simTime);
    }
  }  

  /*
    =======================================
    END: event creation
    =======================================
  */



  /*
    =======================================
    BEGIN: photon generator
    =======================================
  */
  
  // utl::RandomEngine *fRandomEngine3;
  // RandomEngine::RandomEngineType* const randomEngine = &fRandomEngine3->GetEngine();
  for (fevt::FEvent::EyeIterator iEye = fEvent.EyesBegin(fevt::ComponentSelector::eInDAQ);
       iEye != fEvent.EyesEnd(fevt::ComponentSelector::eInDAQ); ++iEye) {
    const int eyeId = iEye->GetId();
    //    cout << eyeId << endl;
    const fdet::Eye& eye = Detector::GetInstance().GetFDetector().GetEye(eyeId);
    const CoordinateSystemPtr eyeCS = eye.GetEyeCoordinateSystem();

    for (fevt::Eye::TelescopeIterator iTel = iEye->TelescopesBegin(fevt::ComponentSelector::eInDAQ);
         iTel != iEye->TelescopesEnd(fevt::ComponentSelector::eInDAQ); ++iTel) {
      const int telId = iTel->GetId();
      //      const int index = Index(eyeId, telId);
      // cout << telId << endl;
      fevt::TelescopeSimData& telSim = iTel->GetSimData();
      telSim.ClearPhotons();

      const fdet::Telescope& detTel = detFD.GetTelescope(*iTel);
//      const double rDia = detTel.GetDiaphragmRadius();
//      const double binWidth = detTel.GetCamera().GetFADCBinSize();
//      const double normWavelength = detFD.GetReferenceLambda();
//      const double fRDiaMin = rDia/10. ;//minimum radius on diaphragm, lets say a tenth of the radius
//      const CoordinateSystemPtr& telCS = detTel.GetTelescopeCoordinateSystem();
      double detFOV = detTel.GetCamera().GetFieldOfView()/deg;//~23 deg
      detFOV += 10;//add 10 degrees to FOV so that more telescopes will be triggered. 
      double elvesTheta  = 90.-elvesWantedLocationUTM.GetPoint(referenceCS).GetTheta(referenceCS)/deg;
      double elvesPhi = elvesWantedLocationUTM.GetPoint(referenceCS).GetPhi(referenceCS)/deg;
      double elvesR = elvesWantedLocationUTM.GetPoint(referenceCS).GetR(referenceCS)/km;

      Vector telAxis = detTel.GetAxis();
      double telFOVLow =  telAxis.GetPhi(referenceCS)/deg - detFOV/2.;
      double telFOVHigh =  telAxis.GetPhi(referenceCS)/deg + detFOV/2.;
      
      //if we are outside detector field of view, go to next telescope
      ostringstream info;  
      if (elvesPhi > telFOVLow && elvesPhi < telFOVHigh) {
        cout << "Eye: " <<  eyeId << " Tel: " << telId << endl;
        cout << "Tel Axis wrt refCS (Az,El): " << telAxis.GetPhi(referenceCS)/deg << " "  <<   90.-telAxis.GetTheta(referenceCS)/deg <<  endl;
        cout << "FOV: " << detFOV << endl;
        cout << "ELVES Coord wrt refCS (r, Az, El): " << elvesR << "  " << elvesPhi  << " " << elvesTheta << endl;

        createNewTFile(iEye, iTel);

      } else {
        info << "Eye: " <<  eyeId << " Tel: " << telId << " - ELVES not in field of view";
        INFO(info);
        continue;
      }

    }  
  }

  /*
    =======================================
    END: photon generator
    =======================================
  */

/*
  for (int i = 0; i < telescope_filenames.size(); ++i) {
    createHistogramPlots(telescope_filenames.at(i), 0.00001);
  }
*/
  makePixels(4, 4);
  if (make_histograms) createHistogramPlots(telescope_filenames.at(2), 0.00001);
	if (make_row_traces) makeRowTrace(8, telescope_filenames.at(2));
	if (make_column_traces) {
		for (int i = 5; i <= 15; i+=5) {
			makeColumnTrace(i, telescope_filenames.at(2));
		}
	}

  return eSuccess;
}


VModule::ResultFlag
ELVESPhotonGenerator::Finish()
{
  
  INFO("ELVESPhotonGenerator::Finish()");
  
  return eSuccess;
}


// This function creates a new TFile with all the data of an elve in the coordinate system of the eye and telescope that the function
// is called with.
void ELVESPhotonGenerator::createNewTFile(fevt::FEvent::EyeIterator iEye, fevt::Eye::TelescopeIterator iTel) {
  // Create a new name for the new TFile
  TString new_file_name = old_file_name(0, old_file_name.Length() - 5);
  new_file_name += "eye";
  new_file_name += iEye->GetId();
  new_file_name += "tel";
  new_file_name += iTel->GetId();
  new_file_name += ".root";

  telescope_filenames.push_back(new_file_name);

  ostringstream info;
  // Check if the new file already exists
  if (exists(new_file_name.Data())) {
    info << new_file_name << " already exists, continuing.";
    INFO(info);
    return;
  }
  info << "Creating a new root file for eye " << iEye->GetId() << " tel " << iTel->GetId() << ".";
  INFO(info);

  // Get the eye and telescope CS to use in the coordinate transformations
  const int eyeId = iEye->GetId();
  const fdet::Eye& eye = Detector::GetInstance().GetFDetector().GetEye(eyeId);
  const CoordinateSystemPtr eyeCS = eye.GetEyeCoordinateSystem();

  Detector& detector = Detector::GetInstance();
  const FDetector& detFD = detector.GetFDetector();
//  const int telId = iTel->GetId();
  const fdet::Telescope& detTel = detFD.GetTelescope(*iTel);
//  const CoordinateSystemPtr& telCS = detTel.GetTelescopeCoordinateSystem();
  Double_t telPhi = detTel.GetAxis().GetPhi(eyeCS);

  // Open a new file and create all the branches we need in it
  TFile* fNewFile = new TFile(new_file_name, "RECREATE");
  TTree* fNewTree = new TTree("ELVES", "ELVES");
  Double_t  nphotons, time_det, time_atm, elevation, azimuth, ndensity, atm_corr, geo_corr, nN22P;
  fNewTree->Branch("elevation", &elevation);
  fNewTree->Branch("azimuth", &azimuth);
  fNewTree->Branch("nphotons", &nphotons);
  fNewTree->Branch("nN22P", &nN22P);
  fNewTree->Branch("time_det", &time_det);
  fNewTree->Branch("time_atm", &time_atm);
  fNewTree->Branch("atm_corr", &atm_corr);
  fNewTree->Branch("geo_corr", &geo_corr);

  // Access the data in the old file
  TVector3* position = 0;
  fTree->SetBranchAddress("position", &position);
  fTree->SetBranchAddress("nN22P", &ndensity);
  fTree->SetBranchAddress("time", &time_atm);
  const int nEntries = fTree->GetEntries();

//	vector<ELVESSimData>* elvePoints = new vector<ELVESSimData>(nEntries);
	vector<ELVESSimData> elvePoints(nEntries);

  const ReferenceEllipsoid ellipsoid(ReferenceEllipsoid::Get(ReferenceEllipsoid::eWGS84));
  const CoordinateSystemPtr referenceCS = Detector::GetInstance().GetReferenceCoordinateSystem();
  UTMPoint elvesSimulatedLocationUTM(0.0, 0.0, 5.0*km, ellipsoid);
  UTMPoint elvesWantedLocationUTM(sourceLat, sourceLon, 5.0*km, ellipsoid);
  const CoordinateSystemPtr WantedLocationCS = fwk::LocalCoordinateSystem::Create(elvesWantedLocationUTM.GetPoint());
  const CoordinateSystemPtr SimulatedLocationCS = fwk::LocalCoordinateSystem::Create(elvesSimulatedLocationUTM.GetPoint());

  Int_t previousProgress = -1;
  int loopTimes = nEntries; // Say how much of the data to go through
//  int loopTimes = 100000; // Say how much of the data to go through

  // These loops save the elve points in the new coordinate system up to 600 us and save it into a new TTree file.
  for (int i = 0; i < loopTimes; ++i) {
    fTree->GetEntry(i);

		// get the point in the elve with respect to the source location
    UTMPoint locationCellUTM(ellipsoid.LatitudeLongitudeHeightToPoint((90.*deg)-position->Theta(), position->Phi(), position->Mag()-6370*km), ellipsoid);
    Point posTMP = Point(locationCellUTM.GetPoint(SimulatedLocationCS).GetX(SimulatedLocationCS),
                         locationCellUTM.GetPoint(SimulatedLocationCS).GetY(SimulatedLocationCS),
                         locationCellUTM.GetPoint(SimulatedLocationCS).GetZ(SimulatedLocationCS),
                         WantedLocationCS);

		elvePoints[i].elevation = 90.-posTMP.GetTheta(eyeCS)/deg;
		elvePoints[i].azimuth   = posTMP.GetPhi(eyeCS)/deg-telPhi/deg;
		elvePoints[i].nN22P     = ndensity;
    // nphotons = einstein coefficient * nN22P * time_step
		elvePoints[i].nphotons  = 2*10000000*ndensity*sim_time_step;
		elvePoints[i].time_det  = time_atm + posTMP.GetR(eyeCS)/(kSpeedOfLight*pow(10,9));
		elvePoints[i].time_atm  = time_atm;

    displayProgress(i, 2*loopTimes, previousProgress, "Reading old root file");
  }

	// sort all the points by when they get to the detector
	previousProgress = -1;
  displayProgress(loopTimes, 2*loopTimes, previousProgress, "Sorting points");
	std::sort(elvePoints.begin(), elvePoints.end(), by_time);

  // add all of the points to the new root file
	Double_t min_time = elvePoints[0].time_det;
	for (int i = 0; i < loopTimes; ++i) {
		if (elvePoints.at(i).time_det - min_time < 0.0006) {
			elevation = elvePoints[i].elevation;
			azimuth   = elvePoints[i].azimuth;
			nphotons  = elvePoints[i].nphotons;
			time_det  = elvePoints[i].time_det - min_time;
			time_atm  = elvePoints[i].time_atm;
			nN22P     = elvePoints[i].nN22P;
			atm_corr  = attenuateAtmospheric(elvePoints[i]);
			geo_corr  = attenuateGeometric(elvePoints[i]);
			fNewTree->Fill(); // add the new event to the tree
		}
	  displayProgress(i+loopTimes, 2*loopTimes, previousProgress, "Creating new root file");
	}
	cout << endl;

  fNewTree->Print(); // add all filled events to the file

  fNewFile->Write(); // write the file
  fNewFile->Close(); // close the file

  return;
}

// return the atmospheric attenuation factor
Double_t ELVESPhotonGenerator::attenuateAtmospheric(ELVESSimData &point) {
	// atmospheric attenuation
	return  1/pow(kE,0.6*1/cos((90-point.elevation)*kPi/180));
}

// Return the geometric attenuation factor
Double_t ELVESPhotonGenerator::attenuateGeometric(ELVESSimData &point) {
  Double_t pixel_width = 0.0456; // distance from side to side
	Double_t pixel_area = 6/sqrt(3)*pixel_width/2; // hexagon area is 1/2*apothem*perimeter
	Double_t sphere_surface_area = 4*kPi*pow(kSpeedOfLight*pow(10, 9)*point.time_det, 2);

	// geometric attenuation
	return pixel_area/sphere_surface_area;
}



inline bool ELVESPhotonGenerator::exists (const std::string& name) {
  ifstream f(name.c_str());
  return f.good();
}


// This function prints out a progress bar based on how far through the data the program is.
void ELVESPhotonGenerator::displayProgress (Int_t currentLoop, Int_t totalLoops, Int_t &previousProgress, TString activity = "") {
  float progress = (float)currentLoop / (float)totalLoops;
  if (int(progress*100.0) > previousProgress) {
		ioctl(STDOUT_FILENO, TIOCGWINSZ, &size);
		int window_size = 0;
		if (size.ws_col > 100) window_size = 100;
		else window_size = size.ws_col;
		for (int i = 0; i < window_size-1; ++i) cout << " ";
		cout << "\r";

    previousProgress = int(progress*100.0);
    int barWidth = window_size - activity.Length() - 15;
    cout << "[";
    int pos = barWidth * progress;
    for (int j = 0; j < barWidth; ++j) {
      if (j < pos) cout << "=";
      else if (j == pos) cout << ">";
      else cout << " ";
    }
    cout << "] " << int(progress * 100.0) << " % | " <<  activity << "\r";
    cout.flush();
  }
  return;
}

void ELVESPhotonGenerator::createHistogramPlots(TString file_name, Double_t time_step) {
  TFile* elvesFile = new TFile(file_name);
  TTree* elvesTree = (TTree*)elvesFile->Get("ELVES");

  // We only need to create the canvas and axis once
  TCanvas* c = new TCanvas("c", "Elve In Telescope CS", 600, 600);

  // initialize the histogram
  // TH2D(name, title, nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);

  TH2D* elves_hist = new TH2D("elves_hist","", 50, -20, 20, 50, 0, 30);
	elves_hist->SetXTitle("Azimuth (degrees)");
	elves_hist->SetYTitle("Elevation (degrees)");
//	Double_t bin_maximum = 1.00475*pow(10,-5);
//  TH2D* elves_hist = new TH2D("elves_hist","", 50, min_azi, max_azi, 50, min_ele, max_ele);

  TString time_reference = "time_det";
  double max_time = elvesTree->GetMaximum(time_reference);
  Int_t graphNum = 0;
//  for (Double_t d = 0; d < time_step; d += time_step) {
  for (Double_t d = 0; d <= max_time; d += time_step) {
    elves_hist->Reset();

    TString hist_name = "Telescope CS, Time = ";
    hist_name += Int_t(d*1000000);
    hist_name += "us";
    elves_hist->SetTitle(hist_name);
//		elves_hist->SetMaximum(bin_maximum);

    // fills the histogram
    TString bounds = "nphotons*";
		if (do_geometric_correction) bounds += "geo_corr*";
		if (do_atmospheric_correction) bounds += "atm_corr*";
		bounds += "(" + time_reference + " > ";
    bounds += d;
    bounds += " && " + time_reference + " < ";
    bounds += d + time_step;
    bounds += ")";
    // "elevation:azimuth" goes "y-axis:x-axis"
    elvesTree->Project("elves_hist", "elevation:azimuth", bounds);

		gStyle->SetOptStat("");

    // draw the histogram
    c->Clear();
    elves_hist->Draw("col");
    drawPixels(*c);
    c->Update();

    // save the histogram
    TString plotname = "graphs/h_";
    plotname += file_name(file_name.Length() - 13, 8);
    plotname += "_";
    if (graphNum < 10) plotname += "0";
//    if (graphNum < 100) plotname += "0";
    plotname += graphNum++;
    plotname += ".png";
    c->SaveAs(plotname);
  }

	c->Close();
	return;
}


void ELVESPhotonGenerator::makeColumnTrace(int col, TString file_name) {
	TGaxis::SetMaxDigits(4); // force traces to use scientific notation

  TFile* elvesFile = new TFile(file_name);
  TTree* elvesTree = (TTree*)elvesFile->Get("ELVES");

  TCanvas* c = new TCanvas("c", "Column Trace", 600, 600);
	TLegend* legend = new TLegend(0.5, 0.7, 0.8, 0.85);
	legend->SetNColumns(2);

  TString time_reference = "time_det";
  double max_time = elvesTree->GetMaximum(time_reference);

	vector<TH1D*> cuts;
	TH1D* frame = new TH1D("frame", "", 50, 0, max_time);

	int color = 1;
	// TH1D(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup);
	TString title = "Column ";
	title += col;
	frame->SetTitle(title);
	frame->SetXTitle("time (s)");
	frame->SetYTitle("nphotons / 2us");
	frame->Draw("PLC");

	gStyle->SetOptStat("");

	Double_t trace_max = 0;
//  for (int i = 4; i <= 13; ++i) { // only loop over large enough histograms
	for (int i = 1; i <= 22; ++i) {
	  TString bounds = "nphotons*";
  	if (do_geometric_correction) bounds += "geo_corr*";
  	if (do_atmospheric_correction) bounds += "atm_corr*";
		bounds += "(";
 	  bounds += pixelCuts.at(22*(col-1)+i)->GetName();
 	  bounds += ")";
		cuts.push_back(new TH1D(pixelCuts.at(i)->GetName(), "", 300, 0, max_time));
		cuts.back()->SetAxisRange(0, max_time, "X");
 	  elvesTree->Project(pixelCuts.at(i)->GetName(), "time_det", bounds);

		if (trace_max < cuts.back()->GetMaximum()) trace_max = cuts.back()->GetMaximum();
		cout << cuts.back()->GetMaximum() << endl;
	}

	frame->SetAxisRange(0, trace_max*1.05, "Y");

	for (int i = 1; i <= 22; ++i) {
		if (cuts.at(i-1)->GetMaximum() > trace_max * 0.05) {
			TString legendLabel = "Row ";
			legendLabel += i;
			legend->AddEntry(cuts[i-1], legendLabel);

			if (color == 5) color++;
			cuts[i-1]->SetLineColor(color++);
			cuts[i-1]->SetLineWidth(3);
 		 	cuts[i-1]->Draw("C same");
 		 	c->Update();
		 	cout << "Row " << i << " drawn." << endl;
		}
	}

	legend->Draw();
	c->Update();

	TString save_string = "column";
	save_string += col;
	save_string += ".png";
	c->SaveAs(save_string);

	c->Close();
	return;
}


// make the traces for a single row. The traces are paired up by hand
void ELVESPhotonGenerator::makeRowTrace(int row, TString file_name) {
	if (row != 8) {
		cout << "functionality for rows other than 8 not implemented" << endl;
		return;
	}

	TGaxis::SetMaxDigits(4); // force traces to use scientific notation

  TFile* elvesFile = new TFile(file_name);
  TTree* elvesTree = (TTree*)elvesFile->Get("ELVES");

  TCanvas* c = new TCanvas("c", "Row Trace", 600, 600);
	TLegend* legend = new TLegend(0.5, 0.65, 0.8, 0.85);
	legend->SetNColumns(2);

  TString time_reference = "time_det";
  double max_time = elvesTree->GetMaximum(time_reference);

	vector<TH1D*> cuts;
	TH1D* frame = new TH1D("frame", "", 50, 0, max_time/2);

	TString frame_title = "Row ";
	frame_title += row;
	frame->SetTitle(frame_title);

	frame->SetXTitle("time (s)");
	frame->SetYTitle("nphotons / 2us");

	frame->Draw("PLC");

	Double_t trace_max = 0;
	for (int i = 1; i <= 440; ++i) {
		if (i%22 == row) {
			if (i/22==1||i/22==18||i/22==3||i/22==16||i/22==5||i/22==14||i/22==7||i/22==12||i/22==9||i/22==10) continue;
			TString bounds = "nphotons*";
	  	if (do_geometric_correction) bounds += "geo_corr*";
	  	if (do_atmospheric_correction) bounds += "atm_corr*";
			bounds += "(";
			bounds += pixelCuts.at(i)->GetName();
			bounds += ")";
			cuts.push_back(new TH1D(pixelCuts.at(i)->GetName(), "", 150, 0, max_time));
			cuts.back()->SetAxisRange(0, max_time, "X");
			elvesTree->Project(pixelCuts.at(i)->GetName(), "time_det", bounds);
			if (trace_max < cuts.back()->GetMaximum()) trace_max = cuts.back()->GetMaximum();
			cuts.back()->Smooth();
		
			gStyle->SetOptStat("");
			if (i/22>=10) cuts.back()->SetLineStyle(2);
			else cuts.back()->SetLineStyle(1);

			// make opposite sides different shades of the same color
			if (i/22==0 || i/22==19) cuts.back()->SetLineColor(kBlack);
			else if (i/22==2 || i/22==17) cuts.back()->SetLineColor(kBlue);
			else if (i/22==4 || i/22==15) cuts.back()->SetLineColor(kMagenta);
			else if (i/22==6 || i/22==13) cuts.back()->SetLineColor(kGreen);
			else if (i/22==8 || i/22==11) cuts.back()->SetLineColor(kRed);
			cuts.back()->SetLineWidth(2);
 	 	  cuts.back()->Draw("C same");
 	 	  c->Update();
			cout << cuts.back()->GetMaximum() << endl;
	 	  cout << "Column " << i/22 << " drawn." << endl;
		}
	}
	frame->SetAxisRange(0, trace_max*1.05, "Y"); // scale for row 8

	legend->AddEntry(cuts.at(0), "Col 1");
	legend->AddEntry(cuts.at(9), "Col 20");
	legend->AddEntry(cuts.at(1), "Col 3");
	legend->AddEntry(cuts.at(8), "Col 18");
	legend->AddEntry(cuts.at(2), "Col 5");
	legend->AddEntry(cuts.at(7), "Col 16");
	legend->AddEntry(cuts.at(3), "Col 7");
	legend->AddEntry(cuts.at(6), "Col 14");
	legend->AddEntry(cuts.at(4), "Col 9");
	legend->AddEntry(cuts.at(5), "Col 12");
	legend->Draw();

	c->Update();
	TString save_name = "row";
	save_name += row;
	save_name += ".png";
	c->SaveAs(save_name);

  c->Close();

  return;
}

void ELVESPhotonGenerator::drawPixels(TCanvas& c) {
	for (unsigned int i = 0; i < pixelCuts.size(); i++) {
		pixelCuts[i]->Draw("L");
	}
	c.Update();
	return;
}

void ELVESPhotonGenerator::makePixels(int eyeId, int telId) {
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

  return;
}
