#include "FitTools.h"
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>

#define RAD2DEG 180./M_PI

int main(int argc, char* argv[])
{
    int nrun = atoi(argv[1]);
    bool bCorrel = atoi(argv[2]);    // Set to 0
    std::cout << "Run " << nrun << " now analysed..." << std::endl;

    char* text_dir = getenv( "TREND_TEXT_PATH" );
    if ( text_dir == NULL )  {
        std::cout << "Environment variable TREND_TEXT_PATH is not defined. Using current directory instead." << std::endl;
        text_dir="./";
    }

    // Load antennas positions
    int antenna;
    int init_ant;
    double x;
    double y;
    double z;
    std::vector<double> coord(900,0);
    std::ostringstream file_antennapos;
    file_antennapos << text_dir << "/coord_antennas.txt"; 
    FILE* fid_antennapos  = fopen( file_antennapos.str().c_str(),   "r" );
    std::cout << "Loading antennas positions file " << file_antennapos.str().c_str() << std::endl;
    int i = 0;
    while ( fscanf( fid_antennapos, "%d  %lf %lf %lf ", &antenna,&x,&y,&z )>0 )  {
       if (i==0)  {  // first line
  	 init_ant = antenna;
       }
       coord.at(i)=x;
       coord.at(i+1)=y;
       coord.at(i+2)=z;
       i=i+3;
       std::cout << "Antenna " << antenna << ", ( " << x << ", " << y << ", " << z << ")" << std::endl;
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
    FILE* fid_out ;
    FitTools fBox;    // The fit toolbox.
    //std::ostringstream file_out;
    
  	// Define variables
    int j;
    int k;
    int iCoinc;
    int iCoinc_prev;
    int Id_antenna;
    int index;
    int evt_antenna;
    int flag;
    double trig_time;
    double trig_time_prev;
    double trig;
    double correl;
    double correlation;
    double correlation_prev;
    double amp;
    double gain;
    double calamp;
    // Loop on analysis technic
    for (int typ=1; typ<3; typ++)  {
      	
      fid_in  = fopen( file_in.str().c_str(),	"r" );  // Need to be done inside loop otherwise stuck at the end of the coinctable file...
      std::ostringstream file_out;
	// Fit settings
 	if (typ==1)  {
 	  file_out << text_dir << "/R" << nrun << "_sphrecons.txt";
	  std::cout << "Now writing spherical wave recons results to file " << file_out.str() << std::endl;
  	  fBox.fitModel()  = POINT_SOURCE; // set either to PLAN_WAVE or POINT_SOURCE
 	}
 	else  {
 	  file_out << text_dir << "/R" << nrun << "_planerecons.txt";
	  std::cout << "Now writing plane wave recons results to file " << file_out.str() << std::endl;
  	  fBox.fitModel()  = PLAN_WAVE;
 	}
 	fBox.fixedSpeed() = true;	 // true: wave speed is a fixed parameter in the fit / false: wave speed is an additional free parameter in the fit
 	fBox.cr() = 1;        // Speed value to use if fixed. By default it is set to 1.0 at initialisation.
	fBox.sigma_t()     = 3.0; // error on times, in m.
        fid_out = fopen( file_out.str().c_str(), "w" );
 
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
	double calamp2=0;
	double calamp1=0;
	
  	while ( fscanf( fid_in, "%lf %d %d  %d %lf %lf %lf %d %lf %lf %lf %lf", &trig_time,&Id_antenna,&evt_antenna,&iCoinc,&trig,&correl,&correlation,&flag,&amp,&gain,&calamp2,&calamp1 )>0 )  {  // New coinctable.txt format
  	  if (floor(iCoinc/1000)==double(iCoinc)/1000 && trig==0)  {
  	     std::cout << "Processing reconstruction for coinc " << iCoinc << "..." << std::endl;
  	  }	  
  	  if (iCoinc_prev != iCoinc || i==0 )  {  // New coinc
  	      //std::cout << "Now writting recons from coinc " << iCoinc_prev << " to file." << std::endl;
	      // Reconstruction for previous coinc
  	      fBox.Na() = j;
	      if (fBox.Na()>3 && ( bCorrel==0 || bCorrel==1 && fBox.ta(2)+fBox.ta(3)>0 ) )  {   
  		  if ( bCorrel==0 || bCorrel==1 && fBox.ta(2)+fBox.ta(3)+fBox.ta(4)>0) {  // Timing performed for intercorrelation
 
  		      double chi2 = fBox.scan();  // Perform a 3D scan for the fBox.Na() first antennas
  		    
		      if (typ==2)  {	// Plan wave	        
  			fprintf( fid_out, "%6.0d %12.5le %3.0d %12.8le %12.5le %12.8le %12.5le %12.5le %12.5le\n", iCoinc_prev, trig_time_prev, fBox.Na(), fBox.theta()*RAD2DEG, fBox.theta_error()*RAD2DEG, fBox.phi()*RAD2DEG, fBox.phi_error()*RAD2DEG, chi2,fBox.chi2_significance() ); // plan wave parameters
			double thetatrue = 180-fBox.theta()*RAD2DEG;
		      }
  		      else  {  // Spherical wave
  			fprintf( fid_out, "%6.0d %12.5le %3.0d %18.15le %18.15le %18.15le %18.15le %12.5le %12.5le\n", iCoinc_prev, trig_time_prev, fBox.Na(), fBox.xs(), fBox.ys(), fBox.zs(), fBox.ts(), chi2,fBox.chi2_significance() ); // output best guess result
  		        //std::cout << "Coinc " << iCoinc_prev << ", typ  = " << typ << ", Nants = " << fBox.Na() << ", time 2nd antenna = " << fBox.ta(2) << ", source reconstructed at position (" << fBox.xs()  << ", " << fBox.ys() << ", " << fBox.zs() << ")" << std::endl;
  		      }
  		  }
  	      }
 
  	      // Re-initialize variables
  	      j=1;
	      k=0;
  	      iCoinc_prev = iCoinc;
  	      trig_time_prev = trig_time;
  	      if (bCorrel)  {
	      	correlation_prev = correlation;
	      }
  	  }
  	  else  {   // Same coinc
  		j++;
  	  }
	  
 	  index = Id_antenna-init_ant;  // Note:this means that ALL antennas have to be given in table starting from init_ant.
           //printf("Id_antenna=%d\n",Id_antenna);
           //printf("init_ant=%d\n",init_ant);
           //printf("index=%d\n",index);
	  
  	  // Now read antenna positions
	  fBox.xa(j) = coord.at(index*3);
  	  fBox.ya(j) = coord.at(index*3+1);
  	  fBox.za(j) = coord.at(index*3+2);
	  // Load time delay
	  if (bCorrel)  {
  	      fBox.ta(j) = correl*299792458*5e-9;  // in meters
  	  }
  	  else {
  	      fBox.ta(j) = trig*299792458*5e-9;
  	  }
	  
  	}  // while
	
	// Reconstruction for last coinc
  	fBox.Na() = j;
	if ( fBox.Na()>3  && ( bCorrel==0 || bCorrel==1 && fBox.ta(2)+fBox.ta(3)>0 ) ) {  // Timing performed for intercorrelation
  	   double chi2 = fBox.scan();  // Perform the source reconstruction for the fBox.Na() first antennas
  	   if (typ==2)  {   // Plan wave
  	   fprintf( fid_out, "%6.0d %12.5le %3.0d %12.8le %12.5le %12.8le %12.5le %12.5le %12.5le\n", iCoinc_prev, trig_time_prev, fBox.Na(), fBox.theta()*RAD2DEG, fBox.theta_error()*RAD2DEG, fBox.phi()*RAD2DEG, fBox.phi_error()*RAD2DEG, chi2,fBox.chi2_significance() ); // plan wave parameters
           //fprintf( fid_out, "%6.0d %12.5le %3.0d %18.15le %18.15le %18.15le %18.15le %12.5le %12.5le %12.5le\n", iCoinc_prev, trig_time_prev, fBox.Na(), fBox.theta()*RAD2DEG, fBox.theta_error()*RAD2DEG, fBox.phi()*RAD2DEG, fBox.phi_error()*RAD2DEG, chi2,fBox.chi2_significance() ); // plan wave parameters
			}
  	   else  {  // Spherical wave
  	   fprintf( fid_out, "%6.0d %12.5le %3.0d %18.15le %18.15le %18.15le %18.15le %12.5le %12.5le\n", iCoinc_prev, trig_time_prev, fBox.Na(), fBox.xs(), fBox.ys(), fBox.zs(), fBox.ts(), chi2,fBox.chi2_significance() ); // output best guess result
           //fprintf( fid_out, "%6.0d %12.5le %3.0d %15.12le %15.12le %15.12le %15.12le %12.5le %12.5le\n", iCoinc_prev, trig_time_prev, fBox.Na(), fBox.xs(), fBox.ys(), fBox.zs(), fBox.ts(), chi2,fBox.chi2_significance() ); // output best guess result
  		    }	   
  	}
	
  	// Close I/Os
  	fclose( fid_out );
      fclose( fid_in  );
    }  // loop on typ
  return 0;
}
