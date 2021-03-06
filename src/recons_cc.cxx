#include "FitTools.h"
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>

#define RAD2DEG 180./M_PI

int main(int argc, char* argv[])
{

  bool bCorrel = 1;  // Use inter-correlation results rather than standard technique
  int bCand = 1;  // Use coinctable.txt file fors candidate search or standards 

  FitTools fBox;    // The fit toolbox.
  /////********************************//////
  for (int i=1; i<argc;i++) {  // loop on runs

    int nrun = atoi(argv[i]);
    std::cout << "Run " << nrun << " now analysed..." << std::endl;

    // Load antennas positions
    int antenna;
    int init_ant;
    double x;
    double y;
    double z;
    std::vector<double> coord(300,0);
    std::ostringstream file_antennapos;
    file_antennapos << "../config/"; 
    if ( nrun > 0 && nrun <= 173 )  {  // North prototype
      file_antennapos << "coord_north.txt";
    }
    if ( nrun > 173 && nrun <= 241 )  {	// east prototype - 1rst positions
  	file_antennapos << "coord_east.txt.R240";
    }
    if ( nrun > 241 & nrun <= 550 )  {  // east prototype - 1rst positions (larger span)
  	file_antennapos << "coord_east.txt";
    }
    if ( nrun > 550 & nrun <= 624 )  {  // east prototype - 2nd positions (South)
  	file_antennapos << "coord_east2.txt";
    }
    if ( nrun > 624 & nrun<800 )  {  // east prototype - 3rd positions (South + 114 higher)
  	file_antennapos << "coord_east3.txt";
    }
    if ( nrun > 800 & nrun<1000)  {  // east prototype - back to 2nd positions (South + 114 at ground)
  	file_antennapos << "coord_east2.txt";
    }
    if ( nrun>1000 & nrun<1300)  {  // east prototype - 4th positions (November tests)
      	file_antennapos << "coord_crosspoint_test.txt";
    }
    if ( nrun>1500  & nrun<1727)  {  // CrossPoint
      	file_antennapos << "coord_crosspoint.txt";
    }
    if ( nrun>1726  & nrun<1740)  {  // CrossPoint witout 139
      	file_antennapos << "coord_crosspoint2.txt";
    }
    if ( nrun>1740 & nrun<2400)  {  // CrossPoint with new position for 149
      	file_antennapos << "coord_crosspoint3.txt";
    }
    if ( nrun>2400)  {  // CrossPoint with new position for 149
      	file_antennapos << "coord_east40.txt";
    }

    if ( nrun == 0 || nrun == 1000 )  { // simu
  	//file_antennapos << "coord_simus.txt";
	///file_antennapos << "coord_north.txt";
	//file_antennapos << "coord_east.txt";   // January
	file_antennapos << "coord_east2.txt";  // March (114 at ground)
	//file_antennapos << "coord_east3.txt";  // March (114 high)
	//file_antennapos << "coord_crosspoint.txt";  
	//file_antennapos << "coord_east20.txt";  // 17 antennas at East pos...
	
    }
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
    }
	
    /////********************************//////
    // Loop on analysis technic
    for (int typ=1; typ<3; typ++)  {
      	
    	// Load data files (I/Os)
    	std::ostringstream file_in;
 	file_in << "../coinctables/R" << nrun << "_coinctable0_" << bCand << ".txt";
 	FILE* fid_in  = fopen( file_in.str().c_str(),	"r" );
 	if (fid_in==0)  {
 	      std::cout << "Could not find file " << file_in.str().c_str() << std::endl;
 	      return 0;
 	}
 	else  {
 	      std::cout << "File " << file_in.str().c_str() << " successfully opened." <<  std::endl;
 	}
 
 	// Fit settings
 	//===
 	if (typ==1)  {
 	  fBox.fitModel()  = POINT_SOURCE; // set either to PLAN_WAVE or POINT_SOURCE
 	}
 	else  {
 	  fBox.fitModel()  = PLAN_WAVE;
 	}
 	fBox.fixedSpeed() = true;	 // true: wave speed is a fixed parameter in the fit / false: wave speed is an additional free parameter in the fit
 	fBox.cr() = 1;        // Speed value to use if fixed. By default it is set to 1.0 at initialisation.
	fBox.sigma_t()     = 3.0; // error on times, in m.

  	std::ostringstream file_out;
  	if (typ==1)  {
  	  file_out << "../results/R" << nrun << "_sphrecons" << bCand << "_" << bCorrel << ".txt";
	  std::cout << "Writing spherical wave recons results to file " << file_out.str() << std::endl;
  	}
  	else  {
  	  file_out << "../results/R" << nrun << "_planerecons" << bCand << "_" << bCorrel << ".txt";
	  std::cout << "Writing plan wave recons results to file " << file_out.str() << std::endl;
  	}
  	FILE* fid_out = fopen( file_out.str().c_str(), "w" );
 
  	// Loop on antenna data
  	//===
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
	
  	while ( fscanf( fid_in, "%lf %d %d  %d %lf %lf %lf %d %lf %lf", &trig_time,&Id_antenna,&evt_antenna,&iCoinc,&trig,&correl,&correlation,&flag,&amp,&gain )>0 )  {  // New coinctable.txt format
  	  if (floor(iCoinc/100)==double(iCoinc)/100 && trig==0)  {
  	     std::cout << "Processing reconstruction for coinc " << iCoinc << "..." << std::endl;
  	  }	  
  	  if (iCoinc_prev != iCoinc || i==0 )  {  // New coinc
  	      
	      // Reconstruction for previous coinc
  	      fBox.Na() = j;
	      if (fBox.Na()>3 && correlation_prev>0)  {   
  		  if ( bCorrel==0 || bCorrel==1 && fBox.ta(2)+fBox.ta(3)+fBox.ta(4)>0) {  // Timing performed for intercorrelation
 
  		      double chi2 = fBox.scan();  // Perform a 3D scan for the fBox.Na() first antennas
  		      //std::cout << "Source reconstructed at position (" << fBox.xs() << ", " << fBox.ys() << ", " << fBox.zs() << ")" << std::endl;
  		      
		      if (typ==2)  {	// Plan wave	        
  			fprintf( fid_out, "%6.0d %12.5le %3.0d %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n", iCoinc_prev, trig_time_prev, fBox.Na(), fBox.theta()*RAD2DEG, fBox.theta_error()*RAD2DEG, fBox.phi()*RAD2DEG, fBox.phi_error()*RAD2DEG, chi2,fBox.chi2_significance() ); // plan wave parameters
			double thetatrue = 180-fBox.theta()*RAD2DEG;
		      }
  		      else  {  // Spherical wave
  			fprintf( fid_out, "%6.0d %12.5le %3.0d %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n", iCoinc_prev, trig_time_prev, fBox.Na(), fBox.xs(), fBox.ys(), fBox.zs(), fBox.ts(), chi2,fBox.chi2_significance() ); // output best guess result
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
  	   fprintf( fid_out, "%6.0d %12.5le %3.0d %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n", iCoinc_prev, trig_time_prev, fBox.Na(), fBox.theta()*RAD2DEG, fBox.theta_error()*RAD2DEG, fBox.phi()*RAD2DEG, fBox.phi_error()*RAD2DEG, chi2,fBox.chi2_significance() ); // plan wave parameters
			}
  	   else  {  // Spherical wave
  	   fprintf( fid_out, "%6.0d %12.5le %3.0d %12.5le %12.5le %12.5le %12.5le %12.5le %12.5le\n", iCoinc_prev, trig_time_prev, fBox.Na(), fBox.xs(), fBox.ys(), fBox.zs(), fBox.ts(), chi2,fBox.chi2_significance() ); // output best guess result
  		    }	   
  	}
	
  	// Close I/Os
  	fclose( fid_out );
      fclose( fid_in  );
    }  // loop on typ
  }  //loop on runs
  return 0;
}
