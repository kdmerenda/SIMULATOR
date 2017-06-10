#ifndef _ELVESPhotonGeneratorNS_ELVESPhotonGenerator_h_
#define _ELVESPhotonGeneratorNS_ELVESPhotonGenerator_h_

#include <fwk/VModule.h>
#include <TTree.h>
#include <TCanvas.h>
#include <string>
#include <fevt/Eye.h>
#include <fevt/TelescopeSimData.h>
#include <fevt/FEvent.h>
#include <vector>
#include <TCutG.h>
#include <sys/ioctl.h>


namespace evt {
  class Event;
}

  /**
    \class ELVESPhotonGenerator

    \Module to convert EMP data to Offline Raw.

    \author Kevin-Druis Merenda
    \date 9 Jan 2017
    \version $Id: ELVESPhotonGenerator.h 28525 2016-02-17 18:03:31Z darko $
  */

  class ELVESPhotonGenerator : public fwk::VModule {

  public:

    fwk::VModule::ResultFlag Init();
    fwk::VModule::ResultFlag Run(evt::Event&);
    fwk::VModule::ResultFlag Finish();
    
    std::string GetSVNId() const
    { return std::string("$Id: ELVESPhotonGenerator.h 28525 2016-02-17 18:03:31Z darko $"); }

  private:

	  const bool do_geometric_correction = false;
	  const bool do_atmospheric_correction = false;
    const bool make_histograms = true;
	  const bool make_row_traces = false;
	  const bool make_column_traces = false;

		const Double_t n_grid_cells = 200; // number of horizontal and vertical grid cells in the simulation
		const Double_t radial_grid_size = 500*100; // size of delta_r of the grid cells (centimeters)
		const Double_t sim_radius = 300; // size of simulation radius (kilometers)
		const Double_t radius_earth = 6730; // radius of earth (kilometers)
		const Double_t sim_time_step = 500*pow(10,-9); // timestep of the simulation (nanoseconds)

		struct winsize size;

		struct ELVESSimData {
			Double_t elevation;
			Double_t azimuth;
			Double_t nphotons;
			Double_t time_det;
			Double_t time_atm;
			Double_t nN22P;
		};

    std::vector<TCutG*> pixelCuts;

		const double sourceLat = -36*utl::deg;
		const double sourceLon = -64*utl::deg;

		static bool by_time(ELVESSimData lhs, ELVESSimData rhs) {
			return lhs.time_det < rhs.time_det;
		}

    TString old_file_name;
    
    TTree* fTree;
    TTree* fTreeParameters;
    Double_t minTime;
    Double_t maxTime;

    std::vector<TString> telescope_filenames;

    void createNewTFile(fevt::FEvent::EyeIterator, fevt::Eye::TelescopeIterator);
    void createHistogramPlots(TString, Double_t);
    inline bool exists (const std::string&);
    void displayProgress(Int_t, Int_t, Int_t&, TString);
		Double_t attenuateAtmospheric(ELVESSimData&);
		Double_t attenuateGeometric(ELVESSimData&);

		void drawPixels(TCanvas&);
		void makePixels(int, int);

		void makeRowTrace(int, TString);
    void makeColumnTrace(int, TString);

    REGISTER_MODULE("ELVESPhotonGenerator", ELVESPhotonGenerator);

  };




#endif

// Configure (x)emacs for this file ...
// Local Variables:
// mode: c++
// compile-command: "make -k"
// End:
