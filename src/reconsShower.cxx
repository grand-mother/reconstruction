#include "FitTools.h"
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>

#define RAD2DEG 180./M_PI

int main(int argc, char* argv[])
{
      int nrun = atoi(argv[1]);
      bool bCorrel = atoi(argv[2]);
      bool bCal = atoi(argv[3]);

      std::cout << "Run " << nrun << " now analysed... bCal = " << bCal << std::endl;

      char* text_dir = getenv( "TREND_TEXT_PATH" );
      if ( text_dir == NULL )  {
        std::cout << "Environment variable TREND_TEXT_PATH is not defined. Using current directory instead." << std::endl;
        text_dir="./";
      }
	      	
    	
	// Antenna positions
      int antenna;
      int init_ant;
      double x;
      double y;
      double z;
      std::vector<double> coord(300,0);
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
         //std::cout << "ant Id = " << antenna << "( " << x << ", " << y << ", " << z << ")" << std::endl;
      }
	
   	  
      // Load data file (I/Os)
      std::ostringstream file_in;
      file_in << text_dir << "/R" << nrun << "_coinctable.txt"; 
      std::cout << "Reading coinctable file " << file_in.str().c_str() << std::endl;
      FILE* fid_in  = fopen( file_in.str().c_str(),   "r" );
      if (fid_in==0)  {
  	    std::cout << "Could not find file " << file_in.str().c_str() << std::endl;
  	    return 0;
      }
  	// Output file
	std::ostringstream file_out;
  	file_out << text_dir << "/R" << nrun << "_shower.txt";
	std::cout << "Writing shower recons results to file " << file_out.str() << std::endl;
  	FILE* fid_out = fopen( file_out.str().c_str(), "w" );


	// Fit settings
 	FitTools fBox;    // The fit toolbox.
  	fBox.fitModel()  = PLAN_WAVE; // set either to PLAN_WAVE or POINT_SOURCE
 	fBox.fixedSpeed() = true;	 // true: wave speed is a fixed parameter in the fit / false: wave speed is an additional free parameter in the fit
 	fBox.cr() = 1;        // Speed value to use if fixed. By default it is set to 1.0 at initialisation.
	fBox.sigma_t()     = 3.0; // error on times, in m.

  	// Init variables
	int k = 0;
  	int iCoinc=0;
  	int iCoinc_prev=0;
  	int Id_antenna=0;
  	int index=0;
  	int evt_antenna=0;
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
	int sat=0;
 		
	// Read coinctable file
	while ( fscanf( fid_in, "%lf %d %d  %d %lf %lf %lf %d %lf %lf %lf %lf", &trig_time,&Id_antenna,&evt_antenna,&iCoinc,&trig,&correl,&correlation,&sat,&amp,&gain,&calamp1,&calamp2 )>0 )  { 
          if (iCoinc_prev != iCoinc || i==0 )  {  // New coinc! Perform recons for previous one.  	      
	      
            // Reconstruction for previous coinc
  	      fBox.Na() = k;
	      if (fBox.Na()>3)  {   
  		  if ( bCorrel==0 || bCorrel==1 && fBox.ta(2)+fBox.ta(3)+fBox.ta(4)>0) {  // Timing performed for intercorrelation
			fBox.fitModel()  = PLAN_WAVE;
  		      double chi2 = fBox.scan();  // Perform a 3D scan for the fBox.Na() first antennas
  		      //std::cout << "Coinc " << iCoinc_prev << " reconstructed at theta = " << 180-fBox.theta()*180/3.14159 << " deg, phi = " << fBox.phi()*180/3.14159-180 << " deg" << std::endl;
  		      fBox.fitModel() = EXPONENTIAL_AMPLITUDE;
		      //std::cout << "Performing shower recons for coinc " << iCoinc_prev << ", Na = " << fBox.Na() << std::endl;
  		      double chi2exp = fBox.scan();  // Perform the shower geometry reconstruction for the xTDoA.Na() first antennas
 
  		      //std::cout << "Coinc " << iCoinc_prev << ": Core position = (" << fBox.x0() << ", " << fBox.y0()  << ", " <<  fBox.z0() << ")" << std::endl;
  		      //std::cout << "Exponential parameters: S0 = " << fBox.s0() << " dB, attenuation = " << 8.6859/fBox.a0()  << " m. " << std::endl;
	
		    	fprintf( fid_out, " %d %d  %15.5le %15.5le %15.5le %15.5le %15.5le %15.5le \n", iCoinc_prev,fBox.Na(),fBox.x0(), fBox.y0(), fBox.z0(),fBox.s0(), 8.6859/fBox.a0(), chi2exp );  //8.6859 = 20/log(10)
 

			//std::cout << "Source reconstructed at position (" << fBox.xs() << ", " << fBox.ys() << ", " << fBox.zs() << ")" << std::endl;
  		      
  		  }  //if bCorrel
  	      } // if Na
 
  	      // Re-initialize variables
	      // k=0;
  	      iCoinc_prev = iCoinc;
  	      trig_time_prev = trig_time;
  	      k = 0;
              if (bCorrel)  {
	      	correlation_prev = correlation;
	      }
          }  // if iCoinc_prev

	    
            if (calamp2!= 0)  {  // this is a shower candidate
  	         
		   
  		   // Load antenna positions
  		   k++;
  		   index = Id_antenna-init_ant;  // Note:this means that ALL antennas have to be given in table starting from init_ant.
  		   fBox.xa(k) = coord.at(index*3);
  		   fBox.ya(k) = coord.at(index*3+1);
  		   fBox.za(k) = coord.at(index*3+2);
  		   // Load time delay
  		   if (bCorrel)  {
  			fBox.ta(k) = correl*299792458*5e-9;  // in meters
  		   }
  		   else {
  			fBox.ta(k) = trig*299792458*5e-9;
  		   }
  		   // Load amplitude
  		   if (bCal == 1) {
                     if (calamp2>0)  { // use PSD calibration
		 	fBox.sa(k) = 20*log(calamp2)/log(10);
                     }
  		   }
		   else  {  // use rms noise calib..
		 	fBox.sa(k) = 20*log(calamp1)/log(10);
		   }
		   //std::cout << "Coinc " << iCoinc << std::endl;
                   //std::cout << "k=" << k << ", Antenna " << Id_antenna << ": xpos = " << coord.at(index*3) << ", ypos = " << coord.at(index*4) << ", t = " << correl << ", amplitude raw " << amp << ", amplitude = " << calamp1  << std::endl;

		   fBox.sat(k) = sat;

		   }  //if calamp
               else  {
                 k = 0;
               }

      }  // while sur coinctable file
      
      fclose( fid_in  );
      
      // Close I/Os
      fclose( fid_out );
      return 0;
} // main

