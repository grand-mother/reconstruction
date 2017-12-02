#ifndef FITTOOLS_H
#define FITTOOLS_H 1

#include "../port_i/src/port_i.h"
#include <math.h>
#include <stdio.h>
#include <setjmp.h>


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// DEFINE
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define MAX_ANTENNA 200  
                         // Note that the memory space for input antennas
                         // data is 'static'. We allocate more space than
                         // we require. This is convenient since the 
			 // minimiser is written in FORTRAN.
			 
#define N_ANTENNA_DATA 6 
                         // Number of antenna related data: position x,y,z,
                         // signal arrival time t, signal amplitude S and
			 // saturation flag.

#define ADD_USER_SPACE 3 
                         // Additional user space for fit data. e.i. to
                         // store the cascade direction in the amplitude fit.
			 
#define MAX_PARAMETERS 5 
                         // Maximum number of parameters for any model.


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			 
// ENUMERATIONS
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
enum fitType
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
{ 
  POINT_SOURCE		=	1, 
  PLAN_WAVE		=	2, 
  EXPONENTIAL_AMPLITUDE	=	4, 
  SPHERICAL_AMPLITUDE	=	8 
};


extern "C" {
  double gradsol_( double*, const int&, const double&, const int& );
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			 
// Toolbox class for fitting (FitTools)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class FitTools
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
public:
  			FitTools();
			~FitTools();

  // References to Point Source parameters (Spherical loss model also)
  //===
  double&		xs() { return _xs; }  // source x position
  double&		ys() { return _ys; }  // source y position
  double&		zs() { return _zs; }  // source z position
  double&		ts() { return _ts; }  // source emitting time
  
  // References to Plane Wave parameters
  //===
  double&		theta() { return _theta; }  // wave vector theta angle 
  double&		phi()   { return _phi;   }  // wave vector phi angle
  
  // Reference to wave speed
  //===
  double&		cr()    { return _cr;    }  // wave normalised speed
  
  // References to Exponential Amplitude Model
  //===
  double&		x0() { return _x0; }  // Impact x coordinate 
  double&		y0() { return _y0; }  // Impact y coordinate
  double&		z0() { return _z0; }  // Impact z coordinate: not in the fit; defined as barycenter of za values.
  double&		a0() { return _a0; }  // Signal characteristic loss (in dB/m: multiply by ln(10)/20 to get exponential loss factor in 1/m)
  
  // Reference to source amplitude (Exponential & Spherical loss models)
  //===
  double&		s0() { return _s0; }  // Signal source amplitude (in dB)
  
  // References to antennas
  //===
  int&			Na() { return _Na; }  // Number of antennas for this minimisation
  double&		xa( const int& n );   // nth antenna x position 
  double&		ya( const int& n );   // nth antenna y position
  double&		za( const int& n );   // nth antenna z position
  double&		ta( const int& n );   // nth antenna signal arrival time
  double&		sa( const int& n );   // nth antenna signal amplitude ( in dB ) 
  double&		sat( const int& n );  // nth antenna saturation flag
  
  // References to error terms (one must call the propagate() method 1st to compute these terms)
  //===
  double&		sigma_t()         { return _sigma_t;      }  // Error on time measurements (in m).
  double		error( const int& i, const int& j ); // Error propagation coefficient on parameter i from the time measurement on antenna j.                                                       
  double		error( const int& i );               // Error coefficient on parameter i for identical and independently distributed  
                                                             // errors on inputs with unit sigma. Multiply by the common sigma value of input.
							     // measurements to get the true parameters error.
							     
  double		theta_error()       { return ( _errorII[ 0 ]*_sigma_t ); }  // Typical error on theta.
  double		phi_error()         { return ( _errorII[ 1 ]*_sigma_t ); }  // Typical error on phi.
  double		chi2_significance() { return _chi2_significance;         }  // Chi2 significance.
  
  int&			gchi2_algorithm() { return _gchi2_algorithm; } // algorithm used for the generalised chi2 computation.
  
  // I/O to fit options
  //===
  fitType&		fitModel()        { return _fitModel;      } // Model to fit
  bool&			fixedSpeed()      { return _fixedSpeed;    } // Is speed fixed or is it a free parameter in the fit?
  bool&			fixedSource()     { return _fixedSource;   } // Is source position fixed or is it a free parameter? (spherical fit only)
  bool&			fixedZs()         { return _fixedZs;       } // Is source z coordinate fixed or is it a free parameter? (spherical fit only)
  
  bool&			computeErrors()   { return _computeErrors; } // To compute errors or not to ...
  
  // Fit methods
  //===
  double		fit( int na=-1 );                     // Single fit
  double		scan( int na=-1 );                    // Multiple fit with a scan of parameter space for initial conditions
  void			errorPropagation( double chi2=-1.0 ); // Compute errors propagation and, if provided, the chi2 significance (only PLAN WAVE fit for now).
  
  double		randn();                                                       // Gaussian normal centred random variable.
  double		erfinv( double );                                              // Inverse error function.
  void 			eigenv( double** A, double* V, const int& );                   // Compute the eigen values V of matrix A. Overwrites eigen vectors in A. 
  double		gchi2cdf( double*, double*, int*, const int&, const double& ); // Generalised chi2 cdf.

private:
  PORT_i		_xPORT; // C++ interface to the FORTRAN minimiser 
  
  fitType		_fitModel;
  bool			_fixedSpeed;
  bool			_fixedSource;
  bool			_fixedZs;
  bool			_computeErrors;

  double		_Xa[ N_ANTENNA_DATA*MAX_ANTENNA + ADD_USER_SPACE ];
  bool 			_sat[ MAX_ANTENNA ];
  int			_Na;
  double		_errorM[ MAX_PARAMETERS ][ MAX_ANTENNA ];
  double		_errorII[ MAX_PARAMETERS ];
  double		_sigma_t;
  double		_chi2_significance;
  
  
  double		_dummyReal;
  bool			_dummyBool;
  
  void			_initFit();
  void			_copyFitParameters();
  int			_nParameters;
  
  double		_scan_TDoA();
  double		_scan_EAM();
  double		_scan_SAM();
  
  void			_tred2( double**, const int&, double*, double* );
  void 			_tqli( double*, double*, const int&, double** );
  double		_pythag ( const double&, const double& );
  
  double 		_imhof( double*, int*, double*, const int&, const double&, double&, 
  			  const double&, const double&, const double&, int& );
  double		_imhofbd( double*, int*, double*, const int&, const double&, 
  			  const double&, int& );
  double 		_imhofint( double*, int*, double*, const int&, const double&, const double&, 
  			  const double&, int& );
  double 		_imhoffn( double*, int*, double*, const int&, const double&, const double& );
  double 		_ruben( double*, int*, double*, const int&, const double&, const double&, int& );
  double 		_qf( double*, double*, int*, const int&, const double&, const double&, 
  			  const int&, const double&, double*, int* ); 
  void 			_qf_counter();
  double 		_qf_exp1( double ); 
  double 		_qf_log1( const double&, const bool& );
  void 			_qf_order();
  double 		_qf_errbd( double, double* );
  double 		_qf_ctff( double, double* );
  double 		_qf_truncation( double, double );
  void 			_qf_findu( double*, double );
  void 			_qf_integrate( int, double, double, bool );
  double 		_qf_cfe( double );
  
  double 		_qf_sigsq, _qf_lmax, _qf_lmin, _qf_mean, _qf_c;
  double 		_qf_intl, _qf_ersm;
  int 			_qf_count, _qf_r, _qf_lim;  
  bool 			_qf_ndtsrt, _qf_fail;
  int 			*_qf_n,*_qf_th; 
  double 		*_qf_lb,*_qf_nc;
  jmp_buf 		_qf_env;
  
  bool			_randn_parity;
  double  		_randn_d1;
  double  		_randn_d2;
  double  		_randn_sq;
  
  int			_gchi2_algorithm;
  int			_mult[  MAX_ANTENNA ];
  double		_lambda[ MAX_ANTENNA ];
  double		_delta[ MAX_ANTENNA ];
  double		_lambda_153[ MAX_ANTENNA+1 ];
  double*		_p_lambda_153;
  double		_trace_155[ 7 ];
  
  double		_xs;
  double		_ys;
  double		_zs;
  double		_ts;
  double		_cs;
  
  double		_theta;
  double		_phi;
  double		_cr;
  
  double		_x0;
  double		_y0;
  double		_z0;
  double		_s0;
  double		_a0;
};

// Objective functions to minimise and their gradients
//===
void PSF_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM ); // Point source
void PSF_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM );
void PWF_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM ); // Plane wave
void PWF_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM );
void EAM_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM ); // Exponential amplitude model
void EAM_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM );
void SAM_function_F( int& N, double* X, int& NF, double* F, int* Na, double* Xa, void* UFPARM ); // Spherical wave amplitude model
void SAM_function_G( int& N, double* X, int& NF, double* G, int* Na, double* Xa, void* UFPARM );

#endif
