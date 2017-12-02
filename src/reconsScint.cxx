#include "FitTools.h"
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
// Copied from recons.cxx
// but lighter program: plane recons only, and no correlation.

#define RAD2DEG 180./M_PI

int main(int argc, char* argv[])
{

    int nrun = atoi(argv[1]);
    std::cout << "Run " << nrun << " now analysed..." << std::endl;

    char* text_dir = getenv( "TREND_TEXT_PATH" );
    if ( text_dir == NULL )  {
        std::cout << "Environment variable TREND_TEXT_PATH is not defined. Using current directory instead." << std::endl;
        //setenv("TREND_TEXT_PATH",".",0); 
        //text_dir = getenv( "TREND_TEXT_PATH" );
        text_dir="./";
    }

    // Load antennas positions
    int antenna;
    int init_ant;
    double x;
    double y;
    double z;    
    std::vector<double> coord(300,0);
    std::ostringstream file_detpos;
    file_detpos << text_dir << "/coord_detectors.txt";
    FILE* fid_detpos  = fopen( file_detpos.str().c_str(),   "r" );
    std::cout << "Loading detector positions file " << file_detpos.str().c_str() << std::endl;
    int i = 0;
	while ( fscanf( fid_detpos, "%d %lf %lf %lf ", &antenna,&x,&y,&z )>0 )  {
	   if (i==0)  {  // first line
  	 		init_ant = antenna;
       }
       coord.at(i)=x;
       coord.at(i+1)=y;
       coord.at(i+2)=z;
       i=i+3;
    }	
      	
      // Load data files (I/Os)
      std::ostringstream file_in;
 	file_in << text_dir << "/R" << nrun << "_coinctable.txt";
 	FILE* fid_in  = fopen( file_in.str().c_str(),	"r" );
 	if (fid_in==0)  {
 	      std::cout << "Could not find file " << file_in.str().c_str() << std::endl;
 	      return 0;
 	}
 	else  {
 	      std::cout << "File " << file_in.str().c_str() << " successfully opened." <<  std::endl;
 	}
  	std::ostringstream file_out;
  	file_out << text_dir << "/R" << nrun << "_hybrid.txt";
	std::cout << "Writing plane wave recons results to file " << file_out.str() << std::endl;
  	FILE* fid_out = fopen( file_out.str().c_str(), "w" );
 
 	FitTools fBox;    // The fit toolbox.
	// Fit settings
 	fBox.fitModel()  = PLAN_WAVE;
 	fBox.fixedSpeed() = true;	 // true: wave speed is a fixed parameter in the fit / false: wave speed is an additional free parameter in the fit
 	fBox.cr() = 1;        // Speed value to use if fixed. By default it is set to 1.0 at initialisation.
	fBox.sigma_t()     = 3.0; // error on times, in m.

 
  	// Initialize variables
  	int j = 1;
	int k =1;
  	int iCoinc=0;
  	int iCoinc_prev=0;
  	int Id_antenna=0;
  	int index=0;
  	int evt_antenna=0;
  	int flag=0;
  	double trig_time=0;
  	double trig_time_prev=0;
  	double trig=0;
  	double correl=0;
  	double correlation=0;
  	double correlation_prev=1;
  	double amp=0;
	double gain=0;
	double calamp1=0;		
	double calamp2=0;		
	
  	while ( fscanf( fid_in, "%lf %d %d  %d %lf %lf %lf %d %lf %lf %lf %lf", &trig_time,&Id_antenna,&evt_antenna,&iCoinc,&trig,&correl,&correlation,&flag,&amp,&gain,&calamp2,&calamp1 )>0 )  {  // New coinctable.txt format
  	  if (floor(iCoinc/100)==double(iCoinc)/100 && trig==0)  {
  	     std::cout << "Processing reconstruction for coinc " << iCoinc << "..." << std::endl;
  	  }	  
  	  if (iCoinc_prev != iCoinc || i==0 )  {  // New coinc
  	      
	      // Reconstruction for previous coinc
  	      fBox.Na() = j;
	      if (fBox.Na()>3)  {   
  		      double chi2 = fBox.scan();  // Perform a 3D scan for the fBox.Na() first antennas
  		      //std::cout << "Source reconstructed at position (" << fBox.xs() << ", " << fBox.ys() << ", " << fBox.zs() << ")" << std::endl;
			  fprintf( fid_out, "%6.0d %12.2le %3.0d %12.2le %12.2le %12.2le %12.2le %12.2le %12.2le\n", iCoinc_prev, trig_time_prev, fBox.Na(), fBox.theta()*RAD2DEG, fBox.theta_error()*RAD2DEG, fBox.phi()*RAD2DEG, fBox.phi_error()*RAD2DEG, chi2,fBox.chi2_significance() ); // plan wave parameters
			  double thetatrue = 180-fBox.theta()*RAD2DEG;
		  }
  		  
  	      // Re-initialize variables
  	      j=1;
	      k=0;
  	      iCoinc_prev = iCoinc;
  	      trig_time_prev = trig_time;
  	  }
  	  else  {   // Same coinc
  		j++;
  	  }
	  
 	  index = Id_antenna-init_ant;  // Note:this means that ALL antennas have to be given in table starting from init_ant.
	  
  	  // Now read antenna positions
	  fBox.xa(j) = coord.at(index*3);
  	  fBox.ya(j) = coord.at(index*3+1);
  	  fBox.za(j) = coord.at(index*3+2);
	  fBox.ta(j) = trig*299792458*5e-9;
  	  
  	}  // while
	
	// Reconstruction for last coinc
  	fBox.Na() = j;
	if ( fBox.Na()>3  ) {  // Timing performed for intercorrelation
  	   double chi2 = fBox.scan();  // Perform the source reconstruction for the fBox.Na() first antennas
  	   fprintf( fid_out, "%6.0d %12.2le %3.0d %12.2le %12.2le %12.2le %12.2le %12.2le %12.2le\n", iCoinc_prev, trig_time_prev, fBox.Na(), fBox.theta()*RAD2DEG, fBox.theta_error()*RAD2DEG, fBox.phi()*RAD2DEG, fBox.phi_error()*RAD2DEG, chi2,fBox.chi2_significance() ); // plan wave parameters
	}
	
  	// Close I/Os
  	fclose( fid_out );
    fclose( fid_in  );
  //}  //loop on runs
  return 0;
}

