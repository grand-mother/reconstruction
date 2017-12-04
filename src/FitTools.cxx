#include "FitTools.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FitTools::FitTools():
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Class constructor. Main task is tuning the PORT_i interfaces handling
//  the minimisation.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_xPORT( MAX_PARAMETERS ) // There are up to MAX_PARAMETERS parameters to minimise.
{


  // Initialisation of antennas times and positions
  //===
  memset( _Xa, 0x0, (N_ANTENNA_DATA*MAX_ANTENNA+ADD_USER_SPACE)
    *sizeof(double) );
  _Na = 0;

  // Initialisation of error matrix.
  //===
  memset( _errorM,  0x0, MAX_PARAMETERS*MAX_ANTENNA*sizeof( double ) );
  memset( _errorII, 0x0, MAX_PARAMETERS*sizeof( double ) );
  _sigma_t = 0.0;

  // Initialisation of gchi2 vars
  //===
  _gchi2_algorithm = 204;
  _lambda_153[ 0 ] = 0.0;
  _p_lambda_153 = _lambda_153+1;

  // Minimiser tuning
  //===
  _xPORT.max_func_evals() = 1000;
  _xPORT.max_iterations() = 1000;
  _xPORT.printing_off(); // Mute the minimiser, comment this for debug

  // Default fit settings
  //===
  _fitModel      = POINT_SOURCE;
  _fixedSource   = false;
  _fixedSpeed    = true;
  _fixedZs       = false;
  _cr            = 1.0;
  _computeErrors = true;

  // Default parameter initialisation
  //===
  _xs = _ys = _zs = _ts = 0.0;
  _theta = _phi = 0.0;
  _x0 = _y0 = _z0 = _s0 = _a0 = 0.0;

  // Initialise random engine
  //===
  srand48( time( NULL ) );
  _randn_parity = false;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FitTools::~FitTools()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void FitTools::_initFit()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define MAX_R 1e10
{
  // Reset constrains
  //===
  for ( int i = 0; i < _xPORT.N(); i++ )
  {
    _xPORT.B( 1, i+1 ) = -MAX_R;
    _xPORT.B( 2, i+1 ) =  MAX_R;
  }

  if ( _fitModel == POINT_SOURCE )
  {
    _nParameters = 4;

    // Connect the objective function to minimise and its gradient
    //===
    _xPORT.f = &PSF_function_F;
    _xPORT.g = &PSF_function_G;

    // Set bounding: Ts <= 0
    //===
    _xPORT.B( 2, 4 ) = 0.0; // Ts parameter (4th) max value (2) is 0.0.
                            // This is to ensure a causal solution.

    // Initialise fit parameters
    //===
    _xPORT.X( 1 ) = _xs;
    _xPORT.X( 2 ) = _ys;
    _xPORT.X( 3 ) = _zs;
    _xPORT.X( 4 ) = _ts;
    _xPORT.X( 5 ) = _cr;
  }
  else if ( _fitModel == PLAN_WAVE )
  {
    _nParameters = 2;

    // Connect the objective function to minimise and its gradient
    //===
    _xPORT.f = &PWF_function_F;
    _xPORT.g = &PWF_function_G;

    // Set bounding: 0 <= theta <= pi
    //===
    _xPORT.B( 1, 1 ) = 0.0;    // theta parameter (1st) min value (1) is 0.0
    _xPORT.B( 2, 1 ) = M_PI;   // theta parameter (1st) max value (2) is pi
    _xPORT.B( 1, 2 ) = 0.0;    // phi parameter (2nd) min value (1) is 0.0
    _xPORT.B( 2, 2 ) = 2*M_PI; // phi parameter (2nd) max value (2) is 2*pi

    // Initialise fit parameters
    //===
    _xPORT.X( 1 ) = _theta;
    _xPORT.X( 2 ) = _phi;
    _xPORT.X( 3 ) = _cr;
  }
  else if ( _fitModel == EXPONENTIAL_AMPLITUDE )
  {
    _nParameters = 4;

    // Connect the objective function to minimise and its gradient
    //===
    _xPORT.f = &EAM_function_F;
    _xPORT.g = &EAM_function_G;

    // Set bounding: a0 >= 0
    //===
    _xPORT.B( 1, 4 ) = 0.0;    // a0 parameter (4th) min value (1) is 0.0

    // Initialise fit parameters
    //===
    _xPORT.X( 1 ) = _x0;
    _xPORT.X( 2 ) = _y0;
    _xPORT.X( 3 ) = _s0;
    _xPORT.X( 4 ) = _a0;

    // Initialise cascade direction
    //===
    _Xa[ N_ANTENNA_DATA*MAX_ANTENNA     ] = -sin( _phi )*sin( _theta );
    _Xa[ N_ANTENNA_DATA*MAX_ANTENNA + 1 ] = cos( _phi )*sin( _theta );
    _Xa[ N_ANTENNA_DATA*MAX_ANTENNA + 2 ] = cos( _theta );
  }
  else if ( _fitModel == SPHERICAL_AMPLITUDE )
  {
    _nParameters = 4;

    // Connect the objective function to minimise and its gradient
    //===
    _xPORT.f = &SAM_function_F;
    _xPORT.g = &SAM_function_G;

    // Set bounding: a0 >= 0
    //===
    _xPORT.B( 1, 1 ) = 0.0; // a0 parameter (1st) min value (1) is 0.0

    // Initialise fit parameters
    //===
    _xPORT.X( 1 ) = _s0;
    _xPORT.X( 2 ) = _xs;
    _xPORT.X( 3 ) = _ys;
    _xPORT.X( 4 ) = _zs;

    if ( _fixedSource ) // Fit only source amplitude
      _nParameters = 1;
    else if ( _fixedZs )
      _nParameters = 3; // Fit source amplitude and x,y coordinates (not z).
  }

  if ( !_fixedSpeed && ( _fitModel == POINT_SOURCE || _fitModel == PLAN_WAVE ) )
  {
    //  Add wave speed as last parameter
    //===
    _nParameters++;

    // Set bounding: cr >= 0
    //===
    _xPORT.B( 1, _nParameters ) = 0.0; // Speed parameter min value is 0.0
  }
}
#undef MAX_R


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void FitTools::_copyFitParameters()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( _fitModel == POINT_SOURCE )
  {
    // Get fit parameters
    //===
    _xs = _xPORT.X( 1 );
    _ys = _xPORT.X( 2 );
    _zs = _xPORT.X( 3 );
    _ts = _xPORT.X( 4 );

    if ( _nParameters == 5 )
      _cr = _xPORT.X( 5 );
  }
  else if ( _fitModel == PLAN_WAVE )
  {
    // Get fit parameters
    //===
    _theta = _xPORT.X( 1 );
    _phi   = _xPORT.X( 2 );
    _cr    = _xPORT.X( 3 );
  }
  else if ( _fitModel == EXPONENTIAL_AMPLITUDE )
  {
    // Get fit parameters
    //===
    _x0 = _xPORT.X( 1 );
    _y0 = _xPORT.X( 2 );
    _s0 = _xPORT.X( 3 );
    _a0 = _xPORT.X( 4 );
  }
  else if ( _fitModel == SPHERICAL_AMPLITUDE )
  {
    // Get fit parameters
    //===
    _s0 = _xPORT.X( 1 );
    _xs = _xPORT.X( 2 );
    _ys = _xPORT.X( 3 );
    _zs = _xPORT.X( 4 );
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& FitTools::xa( const int& n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the nth antenna x position.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( ( n <= 0 ) || ( n > MAX_ANTENNA ) )
  {
    _dummyReal = 0.0;
    return _dummyReal;
  }
  else
    return _Xa[ N_ANTENNA_DATA*( n - 1 ) ];
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& FitTools::ya( const int& n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the nth antenna y position.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( ( n <= 0 ) || ( n > MAX_ANTENNA ) )
  {
    _dummyReal = 0.0;
    return _dummyReal;
  }
  else
    return _Xa[ N_ANTENNA_DATA*( n - 1 ) + 1 ];
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& FitTools::za( const int& n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the nth antenna z position.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( ( n <= 0 ) || ( n > MAX_ANTENNA ) )
  {
    _dummyReal = 0.0;
    return _dummyReal;
  }
  else
    return _Xa[ N_ANTENNA_DATA*( n - 1 ) + 2 ];
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& FitTools::ta( const int& n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the nth antenna signal arrival time.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( ( n <= 0 ) || ( n > MAX_ANTENNA ) )
  {
    _dummyReal = 0.0;
    return _dummyReal;
  }
  else
    return _Xa[ N_ANTENNA_DATA*( n - 1 ) + 3 ];
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& FitTools::sa( const int& n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the nth antenna signal amplitude (in dB).
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( ( n <= 0 ) || ( n > MAX_ANTENNA ) )
  {
    _dummyReal = 0.0;
    return _dummyReal;
  }
  else
    return _Xa[ N_ANTENNA_DATA*( n - 1 ) + 4 ];
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double& FitTools::sat( const int& n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the nth antenna saturation flag
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( ( n <= 0 ) || ( n > MAX_ANTENNA ) )
  {
    _dummyReal = 0.0;
    return _dummyReal;
  }
  else
    return _Xa[ N_ANTENNA_DATA*( n - 1 ) + 5 ];
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::error( const int& i, const int& j )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the error propagation coefficient on parameter i resulting
//  from measurement j.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( ( j <= 0 ) || ( j > MAX_ANTENNA ) ||
       ( i <= 0 ) || ( i > MAX_PARAMETERS ) )
    return ( 0.0 );
  else
    return _errorM[ i-1 ][ j-1 ];
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::error( const int& i )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Reference to the error coefficient on parameter i for identical and
//  independently distributed errors on inputs with unit sigma. Multiply
//  by the common sigma value of input measurements to get the true
//  parameters error.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( ( i <= 0 ) || ( i > MAX_PARAMETERS ) )
    return( 0.0 );
  else
    return _errorII[ i-1 ];
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::fit( int na )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Perform a single fit with user suplied initial conditions. Returns the
//  fit chi2. In case of faillure -1 is returned. The fit best guess is
//  stored in fit parameters vector.
//
//  The number na of antenna can be provided as input argument. The current
//  value of _Na is assumed otherwise.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  // Update number of antennas if required
  //===
  if ( na > 0 )
    _Na = na;

  double t0, z0;
  if ( _fitModel == POINT_SOURCE || _fitModel == PLAN_WAVE )
  {
    // Regularise the time origin
    //===
    t0 = this->ta( 1 );
    for ( int i = 2; i <= _Na; i++ )
      if ( this->ta( i ) < t0 )
        t0 = this->ta( i );
    for ( int i = 1; i <= _Na; i++ )
      this->ta( i )-= t0;
  }
  else if ( _fitModel == EXPONENTIAL_AMPLITUDE )
  {
    // Regularise z coordinate
    //===
    z0 = 0.0;
    for ( int i = 1; i <= _Na; i++ )
      z0+= this->za( i )/_Na;
    for ( int i = 1; i <= _Na; i++ )
      this->za( i )-= z0;
    this->z0() = z0;
  }

  // Initialise according to fit model
  //===
  _initFit();

  // Do the fit
  //===
  double chi2 = -1;
  if ( _xPORT.MNGB( &_Na, _Xa, NULL, _nParameters ) )
  {
    chi2 = _xPORT.value_f();
    this->_copyFitParameters();
  }

  if ( _fitModel == POINT_SOURCE || _fitModel == PLAN_WAVE )
  {
    // Restore the time origin
    //===
    for ( int i = 1; i <= _Na; i++ )
      this->ta( i )+= t0;
    if ( _fitModel == POINT_SOURCE )
    {
      this->ts()+= t0;
      _xPORT.X( 4 )+= t0;
    }
  }
  else if ( _fitModel == EXPONENTIAL_AMPLITUDE )
  {
    // Restore the z origin
    //===
    for ( int i = 1; i <= _Na; i++ )
      this->za( i )+= z0;
  }

  // Normalise chi2
  //===
  if ( ( _fitModel == POINT_SOURCE || _fitModel == PLAN_WAVE ) &&
        ( _sigma_t > 0.0 ) )
    chi2 /= _sigma_t*_sigma_t;
  //std::cout << "Chi2 normalized = " << chi2 << std::endl;
  return chi2;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::scan( int na )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Generic routine for scan like algorithms. The scan tries several fit with
//  various initial conditions. The minimum chi2 over all fits is returned.
//  In case of faillure -1 is returned. The fit best guess is stored as
//  parameter vector.
//
//  The number na of antenna can be provided as input argument. The current
//  value of _Na is assumed otherwise.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  // Update number of antennas if required
  //===
  if ( na > 0 )
    _Na = na;

  // Initialise according to fit model
  //===
  _initFit();

  // Route on specific scan algorithm
  //===
  double chi2;
  if ( _fitModel == POINT_SOURCE || _fitModel == PLAN_WAVE )
  {
    chi2 = this->_scan_TDoA();
    //std::cout << "Best Chi2 = " << chi2  << ", theta = " << this->theta()*180/3.1415927 << " deg , phi = " << this->phi()*180/3.1415927 << " deg. " << std::endl;
    //std::cout << "Best Chi2 = " << chi2  << ", source pos = ( = " << this->xs() << ", " << this->ys() << ", " << this->zs() << "). Time =  " << this->ts()  << std::endl;
  }
  else if ( _fitModel == EXPONENTIAL_AMPLITUDE )
    chi2 = this->_scan_EAM();
  else if ( _fitModel == SPHERICAL_AMPLITUDE )
    chi2 = this->_scan_SAM();

  if ( _computeErrors )
  {
    if ( ( _fitModel == POINT_SOURCE || _fitModel == PLAN_WAVE ) &&
      ( _sigma_t <= 0.0 ) )
      this->errorPropagation();
    else
      this->errorPropagation( chi2 );
  }

  return( chi2 );
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_scan_TDoA()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Scan Algorithm for FitTools algorithms. Tries several fit with various
//  initial conditions. For example for point source localisation initial
//  conditions are taken as a scan of the space in spherical coordinates:
//  r, theta, phi. Steps by 15 degrees in theta, half-quadrant (45 degrees)
//  in phi and logarithmically in r, from 1 m to 10 km.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define N_SPEED  5
#define N_RANGE  4
#define N_THETA  5
#define N_PHI    8
{
  // Scan values
  //===
  double speed_v[] = {  0.20, 0.40, 0.60, 0.80, 1.0  };
  double range_v[] = {  100.,300., 1000.,3000.};
  double theta_v[] = {  15., 45., 60., 75., 90. };
  double phi_v[]   = {  0.,  45., 90., 135., 180., 225., 270., 315.};
  //double phi_v[]   = {  0.,  15., 30., 45., 60., 75., 90., 105., 120., 135., 150., 165., 180., 195., 210., 225., 240., 255., 270., 285., 300., 315., 330., 345.};

  double deg2rad = M_PI/180.;

  // Tune scan parameters according to wave type fit
  //===
  double x0, y0, z0;
  int ir_max, ic_max;
  if ( _fitModel == POINT_SOURCE )
  {
    x0 = y0 = z0 = 0.0;
    for ( int i = 1; i <= _Na; i++ )
    {
      x0+= this->xa( i )/_Na;
      y0+= this->ya( i )/_Na;
      z0+= this->za( i )/_Na;
    }
    ir_max = N_RANGE;
  }
  else if ( _fitModel == PLAN_WAVE )
    ir_max = 1;

  if ( _fixedSpeed )
    ic_max = 1;
  else
    ic_max = N_SPEED;

  // Do the scan
  //==
  double chi2 = -1;
  double disx = 0;
  double disy = 0;
  double disz = 0;
  double dist = 1e60;
  double dista = 0;
  int m;

  std::vector<double> S( _nParameters, 0 );
  for ( int ic = 0; ic < ic_max; ic++ )
    for ( int ir = 0; ir < ir_max; ir++ )
      for ( int it = 0; it < N_THETA; it++ )
        for ( int ip = 0; ip < N_PHI; ip++ )
  {
    if ( !_fixedSpeed )
      this->cr() = speed_v[ ic ];

    if ( _fitModel == POINT_SOURCE )
    {
      this->xs() = -range_v[ ir ]*sin( phi_v[ ip ]*deg2rad )
        *sin( theta_v[ it ]*deg2rad ) + x0;
      this->ys() = range_v[ ir ]*cos( phi_v[ ip ]*deg2rad )
        *sin( theta_v[ it ]*deg2rad ) + y0;
      this->zs() = range_v[ ir ]*cos( theta_v[ it ]*deg2rad ) + z0;

      // Distance to the closest antenna from the hypothetical
      // source position
      //===
      for ( int i = 0; i < this->Na(); i++ )
      {
        disx  = this->xs()-this->xa( i*N_ANTENNA_DATA ) ;
	disy  = this->ys()-this->ya( i*N_ANTENNA_DATA + 1 );
        disz  = this->zs()-this->za( i*N_ANTENNA_DATA + 2 );
        dista = sqrt( disx*disx + disy*disy + disz*disz );
        if ( dista < dist )
          dist = dista;
      }
      this->ts() = -dist/this->cr();
    }
    else if ( _fitModel == PLAN_WAVE )
    {
      this->theta() = theta_v[ it ]*deg2rad;
      this->phi()   = phi_v[ ip ]*deg2rad;
    }

    double d = this->fit();

    if ( ( d >= 0 && d < chi2 ) || ( chi2 < 0 ) )
    {
      //std::cout << "Chi2 improved!" << std::endl;
      //std::cin >> m;  //pause
      chi2   = d;
      for ( int ipar = 0; ipar < _nParameters; ipar++ )
        S[ ipar ] = _xPORT.X( ipar+1 );
    }
  }

  // Restore best parameter guess
  //===
  for ( int ipar = 0; ipar < _nParameters; ipar++ )
    _xPORT.X( ipar+1 ) = S[ ipar ];
  this->_copyFitParameters();

  return chi2;
}
#undef N_PHI
#undef N_THETA
#undef N_RANGE
#undef N_SPEED


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_scan_SAM()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Scan algorithm for amplitude reconstruction model according to
//  1/distance_to_source loss (spherical waves). Initial conditions are
//  taken as a scan of the source position around the reconstructed one
//  in the spherical fit.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define N_RANGE  11
{
  double Dx, Dy, Dz, x0, y0, z0, s0, smax, xmax, ymax, zmax, disx, disy,
         disz, dist;

  int nx, ny, nz;

  Dx = 100.0; // m
  Dy = 100.0; // m
  Dz = 100.0; // m

  if ( _fixedSource )
  {
    nx = ny = nz = 1;
  }
  else if ( _fixedZs )
  {
    nx = ny = N_RANGE;
    nz = 1;
  }
  else
  {
    nx = ny = nz = N_RANGE;
  }

  x0 = this->xs();
  y0 = this->ys();
  z0 = this->zs();

  // Get antenna with max amplitude
  //===
  smax = this->sa( 0 );
  xmax = this->xa( 0 );
  ymax = this->ya( 0 );
  zmax = this->za( 0 );
  for ( int i = 1; i < _Na; i++ )
  {
    if ( this->sa( i ) > smax )
    {
      smax = this->sa( i );
      xmax = this->xa( i );
      ymax = this->ya( i );
      zmax = this->za( i );
    }
  }

  // Scan values
  //===
  double range_v[]  = { -4.0, -2.0, -1.0, -0.50, -0.25, 0.0, 0.25, 0.50,
                         1.0, 2.0, 4.0 };

  // Do the scan
  //===
  double chi2 = -1;

  std::vector<double> S( _nParameters, 0 );
  for ( int ix = 0; ix < nx; ix++ )
    for ( int iy = 0; iy < ny; iy++ )
      for ( int iz = 0; iz < nz; iz++ )
  {
    if ( !_fixedSource )
    {
      this->xs() = x0 + range_v[ ix ]*Dx; // New source x-coord
      this->ys() = y0 + range_v[ iy ]*Dy; // New source y-coord
    }
    if ( !_fixedZs && !_fixedSource )
      this->zs() = z0 + range_v[ iz ]*Dz; // New source z-coord

    disx = this->xs() - xmax;
    disy = this->ys() - ymax;
    disz = this->zs() - zmax;
    dist = sqrt( disx*disx + disy*disy + disz*disz ); // Distance from the source position to the antenna with the highest amplitude

    this->s0() = smax + 20*log(dist)/log(10);         // New source amplitude (assuming max amplitude is correct) in dB

    double d = this->fit();
    if ( ( d >= 0 && d < chi2 ) || ( chi2 < 0 ) )
    {
      chi2   = d;
      for ( int ipar = 0; ipar < _nParameters; ipar++ )
        S[ ipar ] = _xPORT.X( ipar+1 );
    }
  }

  // Restore best parameter guess
  //===
  for ( int ipar = 0; ipar < _nParameters; ipar++ )
    _xPORT.X( ipar+1 ) = S[ ipar ];
  this->_copyFitParameters();

  return chi2;
}
#undef N_RANGE


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_scan_EAM()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Scan algorithm for amplitude reconstruction model according to exponential
//  loss with the lateral distance to the shower axis. Initial conditions are
//  taken as a scan of the impact point over the field of antennas. Loss is
//  scanned for a characteristic length in the range 50-500m.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define N_RANGE  11
#define N_LAMBDA 7
{
  static const double cte = 20.0/log( 10.0 );

  // Get array half widths along x and y directions and get the coordinates
  // of the closest antenna to shower impact point
  //===
  double dmin, dmax, Dx, Dy, x0, y0, s0;
  int i0;

  // Half width along x
  //===
  dmin = dmax =  this->xa( 1 );
  for ( int i = 1; i <= _Na; i++ )
  {
    if ( this->xa( i ) < dmin ) dmin = this->xa( i );
    else if ( this->xa( i ) > dmax ) dmax = this->xa( i );
  }
  Dx = 0.5*( dmax - dmin );

  // Half width along y
  //===
  dmin = dmax =  this->ya( 1 );
  for ( int i = 2; i <= _Na; i++ ) {
    if ( this->ya( i ) < dmin ) dmin = this->ya( i );
    else if ( this->ya( i ) > dmax ) dmax = this->ya( i );
  }
  Dy = 0.5*( dmax - dmin );

  // Coordinates of the closest antenna to the shower impact point
  //===
  i0 = 1;
  dmax =  this->sa( 1 );
  for ( int i = 2; i <= _Na; i++ )
  {
    if ( this->sa( i ) > dmax )
    {
      dmax = this->sa( i );
      i0 = i;
    }
  }
  x0 = this->xa( i0 );
  y0 = this->ya( i0 );
  s0 = this->sa( i0 );

  // Scan values
  //===
  double lambda_v[] = {  50., 100., 150., 200., 300., 400., 500. };
  double range_v[]  = {  -4.0, -2.0, -1.0, -0.50, -0.25, 0.0,
                          0.25, 0.50, 1.0, 2.0, 4.0 };
  // Do the scan
  //===
  double chi2 = -1.0;
  std::vector<double> S( _nParameters, 0 );
  for ( int il = 0; il < N_LAMBDA; il++ )
    for ( int ix = 0; ix < N_RANGE; ix++ )
      for ( int iy = 0; iy < N_RANGE; iy++ )
  {
    double dx = range_v[ ix ]*Dx;
    double dy = range_v[ iy ]*Dy;
    this->x0() = x0 + dx;  // New core x-coord
    this->y0() = y0 + dy;  // New core y-coord
    this->s0() = s0*exp( sqrt( dx*dx + dy*dy )/lambda_v[ il ] );  // New core amplitude (assuming max amplitude is correct)
    this->a0() = cte/lambda_v[ il ];

    double d = this->fit();
    if ( ( d >= 0 && d < chi2 ) || ( chi2 < 0 ) )
    {
      chi2   = d;
      for ( int ipar = 0; ipar < _nParameters; ipar++ )
        S[ ipar ] = _xPORT.X( ipar+1 );
    }
  }

  // Restore best parameter guess
  //===
  for ( int ipar = 0; ipar < _nParameters; ipar++ )
    _xPORT.X( ipar+1 ) = S[ ipar ];
  this->_copyFitParameters();
  return chi2;
}
#undef N_RANGE
#undef N_LAMBDA


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void FitTools::errorPropagation( double chi2obs )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Compute LO error propagation. Build up the error matrix relating input
//  errors on times, DT to output errors on parameters, DT, as:
//
//  DP_i = error( i, j )*DT_j.
//
//  For identical and independently distributed time errors of standard
//  deviation sigma_t(), the standard deviation of the error on parameter
//  P_i goes as:
//
//  sigma_i/sigma_t = error( i ) = sqrt( sum( error( i, j )^2 ), j=1...Na ).
//
//  If a chi2obs value is provided, the routine additionaly computes the
//  significance of the observed chi2 value, assuming identical independently
//  Gausian distributed time errors of standard deviation sigma_t(). The
//  computation is done according to the LO propagation of the error terms.
//
//  NB: Note that the routine is yet only implemented for the plan wave
//  reconstruction algorithm.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( _fitModel != PLAN_WAVE )
    return;

  double A[ 2 ][ 2 ], B[ 2], dK[ 2 ][ 3 ], d1, d2, dx, dy, dz;
  int i, j;

  // Clear error terms
  //===
  memset( _errorM,  0x0, MAX_PARAMETERS*MAX_ANTENNA*sizeof( double ) );
  memset( _errorII, 0x0, MAX_PARAMETERS*sizeof( double ) );

  double ct = cos( _theta );
  double st = sin( _theta );
  double cp = cos(   _phi );
  double sp = sin(   _phi );

  dK[ 0 ][ 0 ] = -ct*sp;
  dK[ 0 ][ 1 ] =  ct*cp;
  dK[ 0 ][ 2 ] =    -st;

  dK[ 1 ][ 0 ] = -st*cp;
  dK[ 1 ][ 1 ] = -st*sp;
  dK[ 1 ][ 2 ] =    0.0;

  // Compute A matrix ( AX = B, with X = parameters error vector
  // and B = measurements error vector )
  //===
  memset( A, 0x0, 4*sizeof( double ) );
  for ( j = 0; j < _Na-1; j++ )
    for ( i = j+1; i < _Na; i++ )
  {
    dx = _Xa[ j*N_ANTENNA_DATA     ] - _Xa[ i*N_ANTENNA_DATA     ];
    dy = _Xa[ j*N_ANTENNA_DATA + 1 ] - _Xa[ i*N_ANTENNA_DATA + 1 ];
    dz = _Xa[ j*N_ANTENNA_DATA + 2 ] - _Xa[ i*N_ANTENNA_DATA + 2 ];

    d1 = dx*dK[ 0 ][ 0 ] + dy*dK[ 0 ][ 1 ] + dz*dK[ 0 ][ 2 ];
    d2 = dx*dK[ 1 ][ 0 ] + dy*dK[ 1 ][ 1 ];

    A[ 0 ][ 0 ] += d1*d1;
    A[ 1 ][ 1 ] += d2*d2;

    A[ 1 ][ 0 ] += d1*d2;
  }
  A[ 0 ][ 1 ] = A[ 1 ][ 0 ];

  // Compute error matrix terms
  //===
  d1 = A[ 0 ][ 0 ]*A[ 1 ][ 1 ] - A[ 0 ][ 1 ]*A[ 1 ][ 0 ];
  if ( d1 == 0.0 )
    return;
  d1 = 1.0/d1;

  for ( j = 0; j < _Na; j++ )
  {
    memset( B, 0x0, 2*sizeof( double ) );
    for ( i = 0; i < _Na; i++ )
    {
      dx = _Xa[ j*N_ANTENNA_DATA     ] - _Xa[ i*N_ANTENNA_DATA     ];
      dy = _Xa[ j*N_ANTENNA_DATA + 1 ] - _Xa[ i*N_ANTENNA_DATA + 1 ];
      dz = _Xa[ j*N_ANTENNA_DATA + 2 ] - _Xa[ i*N_ANTENNA_DATA + 2 ];
      B[ 0 ] += dx*dK[ 0 ][ 0 ] + dy*dK[ 0 ][ 1 ] + dz*dK[ 0 ][ 2 ];
      B[ 1 ] += dx*dK[ 1 ][ 0 ] + dy*dK[ 1 ][ 1 ];
    }

    _errorM[ 0 ][ j ] = ( B[ 0 ]*A[ 1 ][ 1 ] - B[ 1 ]*A[ 0 ][ 1 ] )*d1;
    _errorM[ 1 ][ j ] = ( B[ 1 ]*A[ 0 ][ 0 ] - B[ 0 ]*A[ 1 ][ 0 ] )*d1;

    _errorII[ 0 ] += _errorM[ 0 ][ j ]*_errorM[ 0 ][ j ];
    _errorII[ 1 ] += _errorM[ 1 ][ j ]*_errorM[ 1 ][ j ];
  }
  _errorII[ 0 ] = sqrt( _errorII[ 0 ] );
  _errorII[ 1 ] = sqrt( _errorII[ 1 ] );


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Chi2 significance
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  //  Following a LO expansion of the error terms the chi2 writes:
  //
  //  chi2 = DT^t*C^Dt
  //
  //  where C is a symmetric positive define matrix of correlation terms.
  //  We first diagonalise C over an orthonormal basis such that:
  //
  //  chi2 = sum( lambda_i*u_i^2, i=1...Na )
  //
  //  where the u_i's are centred normal distributed. Then we use a
  //  generalised chi2 algorithm to compute the cdf corresponding to the
  //  observation.
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ( chi2obs == -1.0 )
  {
    _chi2_significance = -1.0;
    return;
  }

  // Build the generalised chi2 correlation matrix
  //===
  double xs[ 3 ], S[ 2 ][ MAX_ANTENNA ];
  double **C = new double*[ _Na ];
  for ( i = 0; i < _Na; i++ )
    C[ i ] = new double[ _Na ];

  memset( xs, 0x0, 3*sizeof( double ) );
  for ( i  = 0; i < _Na; i++ )
  {
    xs[ 0 ]+= _Xa[ i*N_ANTENNA_DATA     ];
    xs[ 1 ]+= _Xa[ i*N_ANTENNA_DATA + 1 ];
    xs[ 2 ]+= _Xa[ i*N_ANTENNA_DATA + 2 ];
  }

  for ( i  = 0; i < _Na; i++ )
  {
    S[ 0 ][ i ] =
      ( xs[ 0 ] - _Na*_Xa[ i*N_ANTENNA_DATA     ] )*dK[ 0 ][ 0 ] +
      ( xs[ 1 ] - _Na*_Xa[ i*N_ANTENNA_DATA + 1 ] )*dK[ 0 ][ 1 ] +
      ( xs[ 2 ] - _Na*_Xa[ i*N_ANTENNA_DATA + 2 ] )*dK[ 0 ][ 2 ];
    S[ 1 ][ i ] =
      ( xs[ 0 ] - _Na*_Xa[ i*N_ANTENNA_DATA     ] )*dK[ 1 ][ 0 ] +
      ( xs[ 1 ] - _Na*_Xa[ i*N_ANTENNA_DATA + 1 ] )*dK[ 1 ][ 1 ];
  }

  for ( j  = 0; j < _Na-1; j++ )
    for ( i  = j+1; i < _Na; i++ )
  {
    C[ i ][ j ] = C[ j ][ i ] =
      _errorM[ 0 ][ i ]*_errorM[ 0 ][ j ]*A[ 0 ][ 0 ] +
      _errorM[ 1 ][ i ]*_errorM[ 1 ][ j ]*A[ 1 ][ 1 ] +
      ( _errorM[ 0 ][ i ]*_errorM[ 1 ][ j ] +
        _errorM[ 1 ][ i ]*_errorM[ 0 ][ j ] )*A[ 0 ][ 1 ] +
      _errorM[ 0 ][ i ]*S[ 0 ][ j ] + _errorM[ 0 ][ j ]*S[ 0 ][ i ] +
      _errorM[ 1 ][ i ]*S[ 1 ][ j ] + _errorM[ 1 ][ j ]*S[ 1 ][ i ] +
      -1.0;
  }

  for ( i  = 0; i < _Na; i++ )
  {
    C[ i ][ i ] =
      _errorM[ 0 ][ i ]*_errorM[ 0 ][ i ]*A[ 0 ][ 0 ] +
      _errorM[ 1 ][ i ]*_errorM[ 1 ][ i ]*A[ 1 ][ 1 ] +
      2.0*_errorM[ 0 ][ i ]*_errorM[ 1 ][ i ]*A[ 0 ][ 1 ] +
      2.0*_errorM[ 0 ][ i ]*S[ 0 ][ i ] +
      2.0*_errorM[ 1 ][ i ]*S[ 1 ][ i ] +
      _Na - 1.0;
  }

  // Compute the generalised chi2 eigen values
  //===
  this->eigenv( C, _lambda, _Na );

  // Protect eigen values and chi2 against roundof errors
  //===
  double lmax, sum;

  lmax = _lambda[ 0 ];
  sum  = 0.0;
  for ( i  = 0; i < _Na; i++ )
  {
    if ( _lambda[ i ] > lmax )
      lmax = _lambda[ i ];
    sum+= _lambda[ i ];
  }

  for ( i  = 0; i < _Na; i++ )
    if ( _lambda[ i ]/lmax < 1e-3 )
      _lambda[ i ] = 0.0;
  if ( chi2obs/sum < 1e-3 )
    chi2obs = 0.0;

  // Build the generalised chi2 coefficients and ndofs
  //===
  int n = 0;
  for ( i = 0; i <_Na; i++ )
  {
    if ( _lambda[ i ] == 0.0 )
     continue;

    _lambda[ n ] = _lambda[ i ];
    _mult[ n ]   = 1;

    for ( j = i+1; j < _Na; j++ )
      if ( fabs( _lambda[ n ] - _lambda[ j ] )/fabs( _lambda[ n ] + _lambda[ j ] ) < 1e-3 )
      {
        _lambda[ j ] = 0.0;
	_mult[ n ]++;
      }
    n++;
  }
  memset( _delta, 0x0, _Na*sizeof( double ) );

  // Compute the chi2 significance
  //===
  _chi2_significance = erfinv( gchi2cdf( _lambda, _delta, _mult, n, chi2obs ) )*sqrt( 2 );

  for ( i = 0; i < _Na; i++ )
    delete[] C[ i ];
  delete[] C;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::randn()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Gaussian normal centred random variable.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( _randn_parity )
  {
    _randn_parity = false;
    return ( _randn_d2 );
  }
  else
    _randn_parity = true;

  do
  {
    _randn_d1 = 2*drand48() - 1;
    _randn_d2 = 2*drand48() - 1;
    _randn_sq = _randn_d1*_randn_d1 + _randn_d2*_randn_d2;
  }
  while ( _randn_sq > 1. || _randn_sq == 0 );

  _randn_sq = sqrt( -2*log( _randn_sq )/_randn_sq );

  _randn_d1*= _randn_sq;
  _randn_d2*= _randn_sq;

  return ( _randn_d1 );
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::erfinv( double p )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Inverse of the error function. Adaptation from an algorithm giving the
//  lower tail quantile standard normal distribution function.
//
//  The original algorithm uses a minimax approximation by rational functions
//  and the result has a relative error whose absolute value is less than
//  1.15e-9.
//
//  Ref: Peter J. Acklam,  http://www.math.uio.no/~jacklam
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  // Coefficients in rational approximations.
  static const double a[] =
  {
    -3.969683028665376e+01,
     2.209460984245205e+02,
    -2.759285104469687e+02,
     1.383577518672690e+02,
    -3.066479806614716e+01,
     2.506628277459239e+00
  };

  static const double b[] =
  {
    -5.447609879822406e+01,
     1.615858368580409e+02,
    -1.556989798598866e+02,
     6.680131188771972e+01,
    -1.328068155288572e+01
  };

  static const double c[] =
  {
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e+00,
    -2.549732539343734e+00,
     4.374664141464968e+00,
     2.938163982698783e+00
  };

  static const double d[] =
  {
     7.784695709041462e-03,
     3.224671290700398e-01,
     2.445134137142996e+00,
     3.754408661907416e+00
  };

  double q, r;

  if ( p <= 0 )
    return ( 0.0 );
  else if ( p >= 1)
    return ( HUGE_VAL );

  p = 0.5*( 1 + p );

  if ( p > 0.97575 ) // Rational approximation for upper region
  {
    q  = sqrt(-2*log(1-p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
	    ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1) /
	    sqrt(2);
  }
  else // Rational approximation for central region
  {
    q = p - 0.5;
    r = q*q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
	   (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1) /
	   sqrt(2);
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void FitTools::eigenv( double** A, double* V, const int& n )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Compute the eigen values V of the n by n square matrix A. The eigen
//  vectors are overwriten to matrix A.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double work[ MAX_ANTENNA ];

  _tred2( A, n, V, work );
  _tqli( V, work, n, A );
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void FitTools::_tred2( double **a, const int& n, double* d, double* e )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Householder reduction of a double, symmetric matrix a[0..n-1][0..n-1]. On
//  output, a is replaced by the orthogonal matrix Q effecting the
//  transformation. d[0..n-1] returns the diagonal elements of the tridiagonal
//  matrix, and e[0..n-1] the off-diagonal elements, with e[0]=0. Several
//  statements, as noted in comments, can be omitted if only eigenvalues are
//  to be found, in which case a contains no useful information on output.
//  Otherwise they are to be included.
//
//  Ref: Numerical Recipes in C, 11.2.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  int l, k, j, i;
  double scale, hh, h, g, f;

  for ( i = n-1; i >= 1; i-- )
  {
    l = i-1;
    h = scale = 0.0;
    if ( l > 0 )
    {
      for ( k = 0; k <= l; k++ )
        scale += fabs( a[ i ][ k]  );
      if ( scale == 0.0 ) // Skip transformation.
        e[ i ] = a [ i ][ l ];
      else
      {
        for ( k = 0; k <= l; k++ )
	{
          a[ i ][ k ] /= scale; // Use scaled a's for transformation.
          h += a[ i ][ k ]*a[ i ][ k ]; // Form sigma in h.
        }
        f = a[ i ][ l ];
        g = ( f >= 0.0 ? -sqrt( h ) : sqrt( h ) );
        e[ i ] = scale*g;
        h -= f*g; // Now h is equation (11.2.4).
        a[ i ][ l ] = f - g; //Store u in the ith row of a.
        f = 0.0;
        for ( j = 0; j <= l; j++ )
	{
          /* Next statement can be omitted if eigenvectors not wanted */
          a[ j ][ i ] = a[ i ][ j ]/h; // Store u/H in ith column of a.
          g = 0.0; // Form an element of A.u in g.
          for ( k = 0; k <= j; k++ )
            g += a[ j ][ k ]*a[ i ][ k ];
          for ( k = j+1; k <= l; k++ )
            g += a[ k ][ j ]*a[ i ][ k ];
          e[ j ] = g/h; // Form element of p in temporarily unused element of e.
	  f += e[ j ]*a[ i ][ j ];
        }
        hh = f/( h + h ); // Form K, equation (11.2.11).
        for ( j = 0; j <= l; j++ ) // Form q and store in e overwriting p.
	{
          f = a[ i ][ j ];
          e[ j ] = g = e[ j ] - hh*f;
          for ( k = 0; k <= j; k++ ) // Reduce a, equation (11.2.13).
          a[ j ][ k ] -= ( f*e[ k ] + g*a[ i ][ k ] );
        }
      }
    }
    else
      e[ i ] = a[ i ][ l ];
    d[ i ] = h;
  }
  /* Next statement can be omitted if eigenvectors not wanted */
  d[ 0 ] = 0.0;
  e[ 0 ] = 0.0;
  /* Contents of this loop can be omitted if eigenvectors not
     wanted except for statement d[ i ] = a [ i ][ i ]; */
  for ( i = 0; i < n; i++ ) // Begin accumulation of transformation matrices.
  {
    l = i - 1;
    if ( d[ i ] ) // This block skipped when i=0
    {
      for ( j = 0; j <= l; j++ )
      {
        g = 0.0;
        for ( k = 0; k <= l; k++ ) // Use u and u/H stored in a to form P.Q.
          g += a[ i ][ k ]*a[ k ][ j ];
        for ( k = 0; k <= l; k++ )
          a[ k ][ j ] -= g*a[ k ][ i ];
      }
    }
    d[ i ] = a[ i ][ i ]; // This statement remains.
    a[ i ][ i ] = 1.0;    // Reset row and column of a to identity
    for ( j = 0; j <= l; j++ ) // matrix for next iteration.
      a[ j ][ i ] = a[ i ][ j ] = 0.0;
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void FitTools::_tqli( double* d, double* e, const int& n, double **z )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  QL algorithm with implicit shifts, to determine the eigenvalues and
//  eigenvectors of a double, symmetric, tridiagonal matrix, or of a double,
//  symmetric matrix previously reduced by tred2 (11.2). On input, d[0..n-1]
//  contains the diagonal elements of the tridiagonal matrix. On output, it
//  returns the eigenvalues. The vector e[0..n-1] inputs the subdiagonal
//  elements of the tridiagonal matrix, with e[0] arbitrary. On output e
//  is destroyed. When finding only the eigenvalues, several lines may be
//  omitted, as noted in the comments. If the eigenvectors of a tridiagonal
//  matrix are desired, the matrix z[0..n-1][0..n-1] is input as the identity
//  matrix. If the eigenvectors of a matrix that has been reduced by tred2
//  are required, then z is input as the matrix output by tred2. In either
//  case, the kth column of z returns the normalized eigenvector
//  corresponding to d[k].
//
//  Ref: Numerical Recipes in C, 11.3.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define SIGN( a, b ) ( (b) >= 0.0 ? fabs(a) : -fabs(a) )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;

  for ( i = 1; i < n; i++ ) // Convenient to renumber the elements of e.
    e[ i-1 ] = e[ i ];
  e[ n-1 ] = 0.0;

  for ( l = 0; l < n; l++ )
  {
    iter=0;
    do
    {
      for ( m = l; m < n-1; m++ ) // Look for a single small subdiagonal element to split the matrix.
      {
        dd = fabs( d[ m ] ) + fabs( d[ m+1 ] );
        if ( (double)( fabs( e[ m ] ) + dd ) == dd )
	  break;
      }
      if ( m != l )
      {
        if ( iter++ == 30 )
	{
	  printf( "Too many iterations in tqli\n" );
	  return;
	}
        g = ( d[ l+1 ] - d[ l ] )/( 2.0*e[ l ] ); // Form shift.
        r = _pythag( g, 1.0 );
        g = d[ m ] - d[ l ] + e[ l ]/( g + SIGN( r, g ) ); // This is dm - ks.
        s = c = 1.0;
        p = 0.0;
        for ( i = m-1; i >= l; i-- ) // A plane rotation as in the original QL, followed by Givens
	{                            // rotations to restore tridiagonal form.
          f = s*e[ i ];
          b = c*e[ i ];
          e[ i+1 ] = ( r = _pythag( f, g ) );
          if ( r == 0.0 ) // Recover from underflow.
	  {
	    d[ i+1 ] -= p;
            e[ m ] = 0.0;
            break;
          }
          s = f/r;
          c = g/r;
          g = d[ i+1 ] - p;
          r = ( d[ i ] - g )*s + 2.0*c*b;
          d[ i+1 ] = g + ( p = s*r );
          g = c*r - b;
          /* Next loop can be omitted if eigenvectors not wanted*/
          for ( k = 0; k < n; k++ ) // Form eigenvectors.
	  {
            f = z[ k ][ i+1 ];
            z[ k ][ i+1 ] = s*z[ k ][ i ] + c*f;
            z[ k ][ i ]   = c*z[ k ][ i ] - s*f;
          }
        }
        if ( r == 0.0 && i >= l )
	  continue;
        d[ l ] -= p;
        e[ l ]  = g;
        e[ m ] = 0.0;
      }
    }
    while ( m != l );
  }
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#undef SIGN
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_pythag ( const double& a, const double& b )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  sqrt(a*a+b*b) with protection against under/overflow
//
//  Ref: Numerical Recipes in C, 2.6.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
   double wa, wb, w;

   wa = fabs ( a );
   wb = fabs ( b );

   if ( wa > wb )
   {
      w = wa;
      wa = wb;
      wb = w;
   }

   if ( wb == 0.0 )
      return 0.0;
   else
   {
      w = wa / wb;
      return wb * sqrt ( 1.0 + w * w );
   }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::gchi2cdf( double* lambda, double* delta, int* ndof,
  const int& n, const double& x )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Computes the cdf of the generalised chi2 distribution:
//
//  gchi2 = sum[ lambda_i*nc_chi2( delta_i, ndof_i ), i=0,n-1 ]
//
//  where nc_chi2( delta, ndof ) is a non centred chi2 distributed variable
//  with ndof degrees of freedom and non-centrality parameter delta.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static const int	order_153	= 12;
  static const double	eps1_256	= 1e-4;
  static const double	eps2_256	= 1e-4;
  static const double	eps3_256	= 1e-4;

  static const double	eps_204		= 1e-4;
  static const double	eps_155		= 1e-4;
  static const int      lim_155         = 10000;

  double 		bound_256  	= 0.0;
  int    		ifault   	= 0;

  if ( ( n == 1 ) && ( _mult[ 0 ] == 1 ) )
    return( erf( x/sqrt( 2 )/lambda[ 0 ] ) );

  if ( _gchi2_algorithm == 256 )
    return ( _imhof( lambda, ndof, delta, n, x, bound_256,
      eps1_256, eps2_256, eps3_256, ifault ) );
  else if ( _gchi2_algorithm == 204 )
    return( _ruben( lambda, ndof, delta, n, x, eps_204, ifault ) );
  else if ( ( _gchi2_algorithm == 155 ) )
    return( _qf( lambda, delta, ndof, n, 0.0, x, lim_155, eps_155, _trace_155, &ifault ) );
  else if ( ( _gchi2_algorithm == 153 ) )
  {
    memcpy( _p_lambda_153, lambda, n*sizeof( double ) ); // lambda[ 0 ] = 0.0 is used to pass
                                                         // an additional parameter not relevant here. See AS153.f for details.
							 // In addition lambda values should be ordered and differents for this algorithm
							 // to work.
    return( gradsol_( _lambda_153, n, x, order_153 ) );
  }
  else
    return ( -1.0 );
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_imhof( double* lambda, int* mult, double* delta,
  const int& noterms, const double& arg, double& bound,
  const double& eps1, const double& eps2, const double& eps3, int& ifault )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  C++ translation of Koerts and Abrahamse's implementation of
//  Imhof's procedure for evaluating the probability that a diagonal
//  form in normal variables is less than a given value, arg.
//
//  Ref: Algorithm AS 256.3  Appl. Statist. (1990) Vol.39, No.2.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define SECURED_IMHOF 0
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  #if ( SECURED_IMHOF )
    // Check parameters
    //===
    if ( ( noterms < 1 ) || ( eps1 <= 0.0 ) ||
         ( eps2 <= 0.0 ) || ( eps3 <= 0.0 ) )
    {
      ifault = 2;
      return( -2.0 );
    }

    ifault = 0;
    for ( int i = 0; i < noterms; i++ )
    {
      if ( ( mult[ i ] < 1 ) || ( delta[ i ] < 0.0 ) )
      {
        ifault = 3;
        return( -(double)i );
      }
    }
  #endif

  // Main body
  //===
  if ( bound <= 0.0 )
    bound = _imhofbd( lambda, mult, delta, noterms, eps1, eps2, ifault );

  double p = _imhofint( lambda, mult, delta, noterms, arg, bound, eps3, ifault );
  if ( ( p < 1e-4 ) || ( ifault != 0 )  )
    p = 0.0;

  return( p );
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_imhofbd( double* lambda, int* mult, double* delta,
  const int& noterms, const double& eps1, const double& eps2, int& ifault )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  C++ translation of Koerts and Abrahamse's procedure
//  for evaluating Imhof's upper bound.
//
//  Ref: Algorithm AS 256.1  Appl. Statist. (1990) Vol.39, No.2.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  int i;
  double count, hold, range, sum1, sum2;

  #if ( SECURED_IMHOF )
    if ( ( noterms < 1 ) || ( eps1 <= 0.0 ) || ( eps2 <= 0.0 ) )
    {
      ifault = 2;
      return( -2.0 );
    }
  #endif

  count = sum1 = sum2 = 0.0;
  for( i = 0; i < noterms; i++ )
  {
    hold = fabs( lambda[ i ] );
    if ( hold > eps2 )
    {
      count += (double)mult[ i ];
      sum1  += mult[ i ]*log( hold );
    }
    sum2 += delta[ i ];
  }

  if ( count < 0.9 )
  {
    ifault = 4;
    return( -4.0 );
  }

  count *= 0.5;
  sum1   = 0.5*sum1 + log( M_PI*count );
  range  = exp( -( sum1 + 0.5*sum2 + log( eps1 ) )/count );

  if ( sum2 == 0.0 )
    range += 5.0/count;
  else
  {
    do
    {
      sum2 = 0.0;
      for ( i = 0; i < noterms; i++ )
      {
        hold  = range*lambda[ i ];
	hold *= hold;
        sum2 += delta[ i ]*hold/( 1.0 + hold );
      }
      hold = exp( sum1 + count*log( range ) + 0.5*sum2 );
      if ( hold*eps1 <= 1.0 )
        range += 5.0/count;
      else
        count = 0.0;
    }
    while ( count != 0.0 );
  }

  return( range );
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_imhofint( double* lambda, int* mult, double* delta,
  const int& noterms, const double& arg, const double& bound,
  const double& eps3, int& ifault )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  C++ translation of Koerts and Abrahamse's procedure
//  for evaluating Imhof's integral by Simpson's Rule.
//
//  Ref: Algorithm AS 256.2  Appl. Statist. (1990) Vol.39, No.2.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static const int maxit = 14;

  int i, j ,n;
  double eps4, int1, int2, step, sum1, sum2, sum4;
  bool cgd;

  #if( SECURED_IMHOF )
    if ( ( noterms < 1 ) || ( bound <= 0.0 ) || ( eps3 <= 0.0 ) )
    {
      ifault = 2;
      return( -2.0 );
    }
  #endif

  ifault = 5;
  n      = 2;
  step   = 0.5*bound;
  eps4   = 3.0*M_PI*eps3;

  sum1   = -arg;
  for ( i = 0; i < noterms; i++ )
    sum1 +=lambda[ i ]*( mult[ i ] + delta[ i ] );
  sum1 = 0.5*sum1 + _imhoffn( lambda, mult, delta, noterms, arg, bound );

  sum4 = _imhoffn( lambda, mult, delta, noterms, arg, step );
  int2 = ( sum1 + 4.0*sum4 )*step;

  sum2 = 0.0;
  for ( i = 1; i <= maxit; i++ )
  {
    n+= n;
    step  = 0.5*step;

    sum2 += sum4;
    sum4  = 0.0;
    for ( j = 1; j <= n; j+= 2 )
      sum4 += _imhoffn( lambda, mult, delta, noterms, arg, j*step );

    int1 = int2;
    int2 = ( sum1 + 2.0*sum2 + 4.0*sum4 )*step;

    if ( i > 3 )
      if ( fabs( int1 - int2 ) < eps4 )
      {
        if ( fabs( int2 ) > 1.5*M_PI )
	  ifault = 6;
	else
	  ifault = 0;
        break;
      }
  }

  return( 0.5 - int2/( 3.0*M_PI ) );
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_imhoffn( double* lambda, int* mult, double* delta,
  const int& noterms, const double& arg, const double& u )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// imhoffn evaluates Imhof's integrand.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  int i;
  double hold, hold2, hold3, rho, sum, theta;

  rho   = 0.0;
  sum   = 0.0;
  theta = -u*arg;

  for ( i  = 0; i < noterms; i++ )
  {
    hold   = u*lambda[ i ];
    hold2  = hold*hold;
    hold3  = 1.0 + hold2;
    theta += mult[ i ]*atan( hold ) + delta[ i ]*hold/hold3;
    sum   += delta[ i]*hold2/hold3;
    rho   += mult[ i ]*log( hold3 );
  }

  return( sin( 0.5*theta )/( u*exp( 0.5*sum + 0.25*rho ) ) );
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#undef SECURED_IMHOF
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_ruben( double* lambda, int* mult, double* delta,
  const int& n, const double& c, const double& eps, int& ifault )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Use of Rubin's (1962) method to evaluate the expression
//  Pr[d k(i) K}(m(i),k(i)}) < c]   where k(i) and c are positive constants,
//  and where K}(m(i),k(i)}) represents an independent chi-squared random
//  variable with m(j) degrees of freedom and non-centrality parameter k(i)}.
//
//   n =  number of chi-squared terms
//   c =  Critical chi-squared value
//   maxit = maximum number of iterations = 500 in Vol 33 No. 3 1984
//   eps = degree of accuracy
//
//  The program returns the cumulative probability value at the point c in
//  the distribution.
//
//  Ref: Program description: Journal of the Royal Statistical Society
//  (Series C) Vol 33 No 3  1984, R.W. Farebrother,  Algorithm AS204,
//  pp 332 - 339.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static const int    maxit  = 500;                  // Maximum number of iterations in eq. 3 of algorithm.
  static const double mode   = 0.0;                  // Set mode = 0.90625 for AS 204A. But 0.0 should be
                                                     // closer to optimal beta choice (see Rubin's).
  static const double lnspi2 = 0.5*log( 0.5*M_PI );

  double tol, beta, ao, aoinv, dnsty, z, rz, eps2, hold, hold2, sum, sum1, dans,
    lans, pans, prbty, temp;
  int itemp, i, j, k, m;
  double gamma[ MAX_ANTENNA ], theta[ MAX_ANTENNA ], a[ maxit ],
    b[ maxit ];
  bool ext, L;

  if ( ( n < 1 ) || ( c <= 0.0 ) || ( maxit < 1 ) || ( eps <= 0.0 ) )
  {
    ifault = 2;
    return( 0.0 );
  }

  tol = -200.0;
  // preliminaries
  beta = sum = lambda[ 0 ];
  for ( i = 0; i < n; i++ )
  {
     if ( ( lambda[ i  ] <= 0.0 ) || ( mult[ i ] < 1 ) || ( delta[ i ] < 0.0 ) )
     {
        ifault = -i;
        return( 0.0 );
     }
     if ( beta > lambda[ i  ] )
       beta = lambda[ i  ];
     if ( sum < lambda[ i ] )
       sum = lambda[ i ];
  }

  if ( mode > 0.0 )
    beta *= mode;
  else
    beta = 2.0/( 1.0/beta + 1.0/sum );

  k    = 0;
  sum  = 1.0;
  sum1 = 0.0;
  for ( i = 0; i < n; i++ )
  {
    hold = beta/lambda[ i ];
    gamma[ i ] = 1.0 - hold;
    for ( j = 0; j < mult[ i ]; j++ )
      sum *= hold;
    sum1 += delta[i];
    k    += mult[ i ];
    theta[ i ] = 1.0;
  }

  ao = exp( 0.5*( log( sum ) - sum1 ) );
  if ( ao <= 0.0 )
  {
    ifault = 1;
    return( 0.0 );
  }

  z = c/beta;
  /* Evaluate probability and density of chi-squared on
     k degrees of freedom. */

  itemp = ( k / 2 ) * 2;
  if (  k == itemp )
  {
    i    = 2;
    lans = -0.5*z;
    dans = exp( lans );
    pans = 1.0 - dans;
  }
  else
  {
    i    = 1;
    lans = -0.5*( z + log( z ) ) - lnspi2;
    dans = exp( lans );
    pans = erf( sqrt( 0.5*z ) );
  }

  k-= 2;

  while( i <= k )
  {
    if ( lans < tol )
    {
      lans += log( z/i );
      dans  = exp( lans );
     }
     else
     {
       temp = dans;
       dans = temp * z/i;
     }

     temp = pans;
     pans = temp - dans;
     i+= 2;
  }

  // Evaluate successive terms of expansion
  prbty = pans;
  dnsty = dans;
  eps2  = eps/ao;
  aoinv = 1.0/ao;
  sum   = aoinv - 1.0;

  for ( m = 0; m < maxit; m++ )
  {
    sum1 = 0.0;
    for ( i = 0; i < n; i++ )
    {
      hold = theta[ i ];
      theta[ i ] *= gamma[ i ];
      sum1 += theta[ i ]*mult[ i ] + m*delta[ i ]*( hold - theta[ i ] );
    }
    sum1 = b[ m ] = 0.5*sum1;
    for ( i = m-1; i >= 0; i-- )
        sum1 += b[ i ]*a[ m - i - 1 ];
    a[ m ] = sum1/( m + 1 );
    sum1 = a[ m ];

    k+= 2;

    if ( lans < tol )
    {
      lans += log( z/k );
      dans  = exp( lans );
    }
    else
    {
      temp = dans;
      dans = temp*z/k;
    }
    pans  -= dans;
    sum   -= sum1;
    dnsty += dans*sum1;
    sum1  *= pans;
    prbty += sum1;

    if ( prbty < -aoinv )
    {
      ifault = 3;
      return ( 0.0 );
    }

    if ( fabs( pans*sum ) < eps2 )
    {
      if ( fabs( sum1 ) < eps2 )
      {
        ifault = 0;
	break;
      }
    }
  }
  if ( m == maxit )
    ifault = 4;

  dnsty *= ao/(beta + beta);
  prbty *= ao;
  if ( ( prbty < 0.0 ) || ( prbty > 1.0 ) )
    ifault += 5;
  else if ( dnsty < 0.0 )
    ifault += 6;

  return( prbty );
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_qf( double* lb, double* nc, int* n, const int& r,
  const double& sigma, const double& c, const int& lim, const double& acc,
  double* trace, int* ifault )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Distribution function of a linear combination of non-central
//  chi-squared random variables :
//
//  input:
//  lb[j]            coefficient of j-th chi-squared variable
//  nc[j]            non-centrality parameter
//  n[j]             degrees of freedom
//  j = 0, 2 ... r-1
//  sigma            coefficient of standard normal variable
//  c                point at which df is to be evaluated
//  lim              maximum number of terms in integration
//  acc              maximum error
//
//  output:
//  ifault = 1       required accuracy NOT achieved
//           2       round-off error possibly significant
//           3       invalid parameters
//           4       unable to locate integration parameters
//           5       out of memory
//
// trace[0]         absolute sum
// trace[1]         total number of integration terms
// trace[2]         number of integrations
// trace[3]         integration interval in final integration
// trace[4]         truncation point in initial integration
// trace[5]         s.d. of initial convergence factor
// trace[6]         cycles to locate integration parameters
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define LOG28 .0866 // log( 2.0 ) / 8.0
//-------------------------------------------------------------------------
#define square( X ) \
  (X)*(X)
//-------------------------------------------------------------------------
#define cube( X ) \
  (X)*(X)*(X)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  int j, nj, nt, ntm;
  double acc1, almx, xlim, xnt, xntm;
  double utx, tausq, sd, intv, intv1, x, up, un, d1, d2, lj, ncj;
  double qfval;
  static const int rats[] = { 1, 2, 4, 8 };

  if ( setjmp( _qf_env ) != 0 ) // On call to longjmp, jump to endofproc
  {
    *ifault = 4;
    goto endofproc;
  }

  _qf_r   = r;
  _qf_lim = lim;
  _qf_c   = c;
  _qf_n   = n;
  _qf_lb  = lb;
  _qf_nc  = nc;

  memset( trace, 0x0, 7*sizeof( double ) );

  *ifault    = 0;
  _qf_count  = 0;
  _qf_intl   = 0.0;
  _qf_ersm   = 0.0;
  qfval      = -1.0;
  acc1       = acc;
  _qf_ndtsrt = true;
  _qf_fail   = false;
  xlim       = (double)_qf_lim;

  _qf_th = (int*)malloc( _qf_r*sizeof( int ) );
  if ( !_qf_th )
  {
    *ifault = 5;
    goto endofproc;
  }

  // Find mean, sd, max and min of lb,
  // check that parameter values are valid.
  _qf_sigsq = square( sigma );
  sd        = _qf_sigsq;
  _qf_lmax  = 0.0;
  _qf_lmin  = 0.0;
  _qf_mean  = 0.0;
  for ( j = 0; j < _qf_r; j++ )
  {
    nj  = _qf_n[  j ];
    lj  = _qf_lb[ j ];
    ncj = _qf_nc[ j ];
    if ( nj < 0  ||  ncj < 0.0 )
    {
      *ifault = 3;
      goto  endofproc;
    }
    sd  += square( lj )*( 2.0*nj + 4.0*ncj );
    _qf_mean += lj*( nj + ncj );
    if ( _qf_lmax < lj )
      _qf_lmax = lj ;
    else if ( _qf_lmin > lj )
      _qf_lmin = lj;
  }

  if ( sd == 0.0  )
  {
    qfval = ( _qf_c > 0.0 ) ? 1.0 : 0.0;
    goto  endofproc;
  }
  if ( ( _qf_lmin == 0.0 ) && ( _qf_lmax == 0.0 ) && ( sigma == 0.0 ) )
  {
    *ifault = 3;
    goto  endofproc;
  }
  sd   = sqrt( sd );
  almx = ( _qf_lmax < -_qf_lmin ) ? -_qf_lmin : _qf_lmax;

  // Starting values for findu, ctff
  utx = 16.0/sd;
  up  = 4.5/sd;
  un  = -up;

  // Truncation point with no convergence factor.
  _qf_findu( &utx, 0.5*acc1 );

  // Does convergence factor help?
  if ( _qf_c != 0.0  && ( almx > 0.07 * sd ) )
  {
    tausq = .25 * acc1 / _qf_cfe( _qf_c );
    if ( _qf_fail )
      _qf_fail = false ;
    else if ( _qf_truncation( utx, tausq ) < .2 * acc1 )
    {
      _qf_sigsq += tausq;
      _qf_findu( &utx, .25 * acc1 );
      trace[ 5 ] = sqrt( tausq );
    }
  }
  trace[ 4 ] = utx;
  acc1       = 0.5 * acc1;

  // Find RANGE of distribution, quit if outside this.
  l1:
    d1 = _qf_ctff( acc1, &up ) - _qf_c;
    if ( d1 < 0.0 )
    {
      qfval = 1.0;
      goto endofproc;
    }
    d2 = _qf_c - _qf_ctff( acc1, &un );
    if ( d2 < 0.0 )
    {
      qfval = 0.0;
      goto endofproc;
    }

    // Find integration interval.
    intv = 2.0 * M_PI / ( ( d1 > d2 ) ? d1 : d2 );

    // Calculate number of terms required for main and
    // auxillary integrations.
    xnt  = utx / intv;
    xntm = 3.0 / sqrt( acc1 );
    if ( xnt > xntm * 1.5 )
    {
      // Parameters for auxillary integration
      if ( xntm > xlim )
      {
        *ifault = 1;
	goto endofproc;
      }
      ntm   = (int)floor( xntm + 0.5 );
      intv1 = utx / ntm;
      x     = 2.0 * M_PI / intv1;
      if ( x <= fabs( _qf_c ) )
        goto l2;

      // Calculate convergence factor.
      tausq = 0.33 * acc1 / ( 1.1 * ( _qf_cfe( _qf_c - x ) + _qf_cfe( _qf_c + x ) ) );
      if ( _qf_fail )
        goto l2;
      acc1 = 0.67 * acc1;

      // Auxillary integration.
      _qf_integrate( ntm, intv1, tausq, false );
      xlim  -= xntm;
      _qf_sigsq += tausq;
      trace[ 2 ] += 1.0;
      trace[ 1 ] += ntm + 1.0;

      // Find truncation point with new convergence factor.
      _qf_findu( &utx, .25 * acc1 );
      acc1 = 0.75 * acc1;
      goto l1;
    }

    // Main integration.
    l2:
      trace[ 3 ] = intv;
      if ( xnt > xlim )
      {
        *ifault = 1;
	goto endofproc;
      }
      nt = (int)floor( xnt + 0.5 );
      _qf_integrate( nt, intv, 0.0, true );
      trace[ 2 ] += 1.0;
      trace[ 1 ] += nt + 1.0;
      qfval       = 0.5 - _qf_intl;
      trace[ 0]   = _qf_ersm;

      // Test whether round-off error could be significant
      // allow for radix 8 or 16 machines.
      up = _qf_ersm;
      x  = up + acc / 10.0;
      for ( j = 0; j < 4; j++ )
        if ( rats[ j ] * x == rats[ j ] * up )
	  *ifault = 2;

    endofproc:
      free( (char*)_qf_th );
      trace[ 6 ] = (double)_qf_count;
      return qfval;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void FitTools::_qf_counter()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Count number of calls to errbd, truncation, cfe.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  _qf_count++;

  if ( _qf_count > _qf_lim )
    longjmp( _qf_env, 1 );
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_qf_exp1( double x )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// To avoid underflows.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  return ( x < -50.0 ? 0.0 : exp( x ) );
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_qf_log1( const double& x, const bool& first )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  If ( first ) log( 1 + x ) ; else  log( 1 + x ) - x.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  if ( fabs( x ) > 0.1 )
    return ( first ? log( 1.0 + x ) : ( log( 1.0 + x ) - x ) );
  else
  {
    double s, s1, term, y, k;

    y     = x / ( 2.0 + x );
    term  = 2.0 * cube( y );
    k     = 3.0;
    s     = ( first ? 2.0 : - x )*y;
    y     = square( y );
    for ( s1 = s + term / k; s1 != s; s1 = s + term / k )
    {
      k    += 2.0;
      term *= y;
      s     = s1;
    }
    return s;
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void FitTools::_qf_order()
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Find order of absolute values of lb
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  int j, k;
  double lj;

  for ( j = 0; j < _qf_r; j++ )
  {
    lj = fabs( _qf_lb[ j ] );
    for ( k = j-1; k >= 0; k-- )
    {
      if ( lj > fabs( _qf_lb[ _qf_th[ k ] ] ) )
        _qf_th[ k + 1] = _qf_th[ k ];
      else goto l1;
    }
    k = -1;
    l1 :
      _qf_th[ k + 1 ] = j;
  }
  _qf_ndtsrt = false;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_qf_errbd( double u, double* cx )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Find bound on tail probability using mgf, cutoff
//  point returned to *cx.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  double sum1, lj, ncj, x, y, xconst;
  int j, nj;

  _qf_counter();

  xconst = u * _qf_sigsq;
  sum1 = u * xconst;
  u *= 2.0;

  for ( j = _qf_r-1; j >= 0; j-- )
  {
    nj  = _qf_n[  j ];
    lj  = _qf_lb[ j ];
    ncj = _qf_nc[ j ];
    x   = u * lj;
    y = 1.0 - x;
    xconst += lj * ( ncj / y + nj ) / y;
    sum1   += ncj * square( x / y )
           + nj * ( square( x ) / y + _qf_log1( -x, false ) );
  }
  *cx = xconst;

  return _qf_exp1( -0.5 * sum1 );
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_qf_ctff( double accx, double* upn )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Find ctff so that p(qf > ctff) < accx,  if (upn > 0),
//  p(qf < ctff) < accx, otherwise.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  double u1, u2, u, rb, xconst, c1, c2;

  u2 = *upn;
  u1 = 0.0;
  c1 = _qf_mean;
  rb = 2.0 * ( ( u2 > 0.0 ) ? _qf_lmax : _qf_lmin );
  for ( u = u2 / ( 1.0 + u2 * rb ); _qf_errbd( u, &c2 ) > accx; u = u2 / ( 1.0 + u2 * rb ) )
  {
    u1 = u2;
    c1 = c2;
    u2 = 2.0 * u2;
  }
  for ( u = ( c1 - _qf_mean ) / ( c2 - _qf_mean ); u < 0.9; u = ( c1 - _qf_mean ) / ( c2 - _qf_mean ) )
  {
    u = ( u1 + u2 ) / 2.0;
    if ( _qf_errbd( u / (1.0 + u * rb ), &xconst ) > accx )
    {  u1 = u;
       c1 = xconst;
    }
    else
    {
      u2 = u;
      c2 = xconst;
    }
  }
  *upn = u2;

  return c2;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_qf_truncation( double u, double tausq )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Bound integration error due to truncation at u.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  double sum1, sum2, prod1, prod2, prod3, lj, ncj,
    x, y, err1, err2;
  int j, nj, s;

  _qf_counter();

  sum1   = 0.0;
  prod2  = 0.0;
  prod3  = 0.0;
  s      = 0;
  sum2   = ( _qf_sigsq + tausq )*u*u;
  prod1  = 2.0*sum2;
  u     *= 2.0;
  for ( j = 0; j < _qf_r; j++ )
  {
    lj    = _qf_lb[ j ];
    ncj   = _qf_nc[ j ];
    nj    = _qf_n[  j ];
    x     = u*u*lj*lj;
    sum1 += ncj*x/( 1.0 + x );

    if ( x > 1.0 )
    {
      prod2 += nj*log( x );
      prod3 += nj*_qf_log1( x, true );
      s     += nj;
    }
    else
      prod1 += nj*_qf_log1( x, true );
  }
  sum1  *= 0.5;
  prod2 += prod1;
  prod3 += prod1;
  x = _qf_exp1( -sum1 - 0.25*prod2 ) / M_PI;
  y = _qf_exp1( -sum1 - 0.25*prod3 ) / M_PI;
  err1 = ( s  ==  0 )  ? 1.0 : x*2.0/s;
  err2 = ( prod3 > 1.0 ) ? 2.5*y : 1.0;
  if ( err2 < err1 )
    err1 = err2;
  x = 0.5*sum2;
  err2 = ( x  <=  y ) ? 1.0 : y/x;

  return ( err1 < err2 ) ? err1 : err2;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void FitTools::_qf_findu( double* utx, double accx )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Find u such that truncation( u ) < accx and truncation( u / 1.2 ) > accx.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  double u, ut;
  int i;
  double divis[] = { 2.0, 1.4, 1.2, 1.1 };

  ut = *utx;
  u  = 0.25*ut;
  if ( _qf_truncation( u, 0.0 ) > accx )
  {
    for ( u = ut; _qf_truncation( u, 0.0 ) > accx; u = ut )
      ut *= 4.0;
  }
  else
  {
    ut = u;
    for ( u = 0.25*u; _qf_truncation( u, 0.0 ) <=  accx; u = 0.25*u )
      ut = u;
  }
  for ( i = 0; i < 4; i++ )
  {
    u = ut/divis[ i ];
    if ( _qf_truncation( u, 0.0 )  <=  accx )
      ut = u;
  }
  *utx = ut;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void FitTools::_qf_integrate( int nterm, double interv, double tausq, bool mainx )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Carry out integration with nterm terms, at stepsize
//  interv.  if (! mainx ) multiply integrand by
//  1.0-exp(-0.5*tausq*u^2).
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  double inpi, u, sum1, sum2, sum3, x, y, z;
  int k, j, nj;

  inpi = interv / M_PI;
  for ( k = nterm; k >= 0; k-- )
  {
    u = ( k + 0.5 )*interv;
    sum1 = - 2.0*u*_qf_c;
    sum2 = fabs( sum1 );
    sum3 = -0.5*_qf_sigsq*square( u );
    for ( j = _qf_r-1; j >= 0; j-- )
    {
      nj    = _qf_n[ j ];
      x     = 2.0*_qf_lb[ j ]*u;
      y     = square( x );
      sum3 -= 0.25*nj*_qf_log1( y, true );
      y     = _qf_nc[ j ]* x/( 1.0 + y );
      z     = nj*atan( x ) + y;
      sum1 += z;
      sum2 += fabs( z );
      sum3 -= 0.5*x*y;
    }

    x = inpi*_qf_exp1( sum3 )/u;
    if ( !mainx )
      x *= 1.0 - _qf_exp1( -0.5*tausq*square( u ) );
    sum1      = sin( 0.5*sum1 )*x;
    sum2     *= 0.5*x;
    _qf_intl += sum1;
    _qf_ersm += sum2;
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double FitTools::_qf_cfe( double x )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Coef of tausq in error when convergence factor of
// exp1(-0.5*tausq*u^2) is used when df is evaluated at x.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  double axl, axl1, axl2, sxl, sum1, lj;
  int j, k, t;

  _qf_counter();

  if ( _qf_ndtsrt )
    _qf_order();
  axl  = fabs( x );
  sxl  = ( x > 0.0 ) ? 1.0 : -1.0;
  sum1 = 0.0;

  for ( j = _qf_r-1; j >= 0; j-- )
  {
    t = _qf_th[ j ];
    if ( _qf_lb[ t ] * sxl > 0.0 )
    {
      lj   = fabs( _qf_lb[ t ] );
      axl1 = axl - lj*( _qf_n[t] + _qf_nc[ t ] );
      axl2 = lj / LOG28;
      if ( axl1 > axl2 )
        axl = axl1;
      else
      {
        if ( axl > axl2 )
	  axl = axl2;
        sum1 = ( axl - axl1 ) / lj;
        for ( k = j-1; k >= 0; k-- )
          sum1 += ( _qf_n[ _qf_th[ k ] ] + _qf_nc[ _qf_th[ k ] ] );
        goto l;
      }
    }
  }

  l:
    if ( sum1 > 100.0 )
    {
      _qf_fail = true;
      return 1.0;
    }
    else
      return pow( 2.0, ( sum1 / 4.0 ) ) / ( M_PI*square( axl ) );
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#undef LOG28
#undef square
#undef cube
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PSF_function_F( int& N, double* X, int& NF, double* F, int* Na,
  double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The objective function to minimise for the point source fit. It is built
//  as a chi2 over all antennas, requiring:
//
//  | Xa - Xs | - cs*( Ta - Ts ) = 0
//
//  Time is assumed to be expressed in m and speed cs is normalised to C0.
//  Note that although the evaluation of |Xa-Xs| is longer than the one of
//  of (Xa-Xs)^2 building the chi2 from |Xa-Xs| is more stable numericaly.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double d1, d2, d3, fi;

  std::cout.precision(15);
  F[ 0 ] = 0.0;
  for ( int i = 0; i < Na[ 0 ]; i++ )
  {
    d1      = ( X[ 0 ] - Xa[ i*N_ANTENNA_DATA     ] );
    d2      = ( X[ 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ] );
    d3      = ( X[ 2 ] - Xa[ i*N_ANTENNA_DATA + 2 ] );
    fi      = X[ 4 ]*( X[ 3 ] - Xa[ i*N_ANTENNA_DATA + 3 ] ) +
                sqrt( d1*d1 + d2*d2 + d3*d3 );
    //std::cout << "** Antenna " << i << ": (" << Xa[ i*N_ANTENNA_DATA     ] << "," << Xa[ i*N_ANTENNA_DATA + 1 ] << "," << Xa[ i*N_ANTENNA_DATA + 2 ] << "): t = " << Xa[ i*N_ANTENNA_DATA + 3  ] << std::endl;
    //std::cout << "dPos = " <<  sqrt( d1*d1 + d2*d2 + d3*d3 ) << ", dTime = " <<  X[ 4 ]*( X[ 3 ] - Xa[ i*N_ANTENNA_DATA + 3 ] ) << ", sum = " << fi << std::endl;
    F[ 0 ] += fi*fi;
  }
  //std::cout << "Source: (" << X[ 0 ] <<", " << X[ 1 ] << ", " << X[ 2 ]  << "), time " << X[ 3 ]  << std::endl;
  //std::cout << "Chi2 = " << F[ 0 ] << std::endl;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PSF_function_G( int& N, double* X, int& NF, double* G, int* Na,
  double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The gradient of the objective function.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double d1, d2, d3, d, fi;

  G[ 0 ] = G[ 1 ] = G[ 2 ] = G[ 3 ] = G[ 4 ] = 0.0;
  for ( int i = 0; i < Na[ 0 ]; i++ )
  {
    d1  = ( X[ 0 ] - Xa[ i*N_ANTENNA_DATA     ] );
    d2  = ( X[ 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ] );
    d3  = ( X[ 2 ] - Xa[ i*N_ANTENNA_DATA + 2 ] );
    d   = sqrt( d1*d1 + d2*d2 + d3*d3 );
    fi  = X[ 4 ]*( X[ 3 ] - Xa[ i*N_ANTENNA_DATA + 3 ] ) + d;

    G[ 0 ] += ( X[ 0 ] - Xa[ i*N_ANTENNA_DATA	  ] )*fi/d;
    G[ 1 ] += ( X[ 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ] )*fi/d;
    G[ 2 ] += ( X[ 2 ] - Xa[ i*N_ANTENNA_DATA + 2 ] )*fi/d;
    G[ 3 ] += X[ 4 ]*fi;
  }
  G[ 0 ] *= 2.0;
  G[ 1 ] *= 2.0;
  G[ 2 ] *= 2.0;
  G[ 3 ] *= 2.0;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PWF_function_F( int& N, double* X, int& NF, double* F, int* Na,
  double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The objective function to minimise for plane wave fit. It is built
//  as a chi2 over all antennas, requiring:
//
//  ( Xa_j - Xa_i ).k - cr*( Ta_j - Ta_i ) = 0
//
//  over all pairs.
//
//  Time is assumed to be expressed in m and speed cr is normalised to C0
//
//  OMH 28/09/09: warning: due to modified conventions ( x=WE, y=SN ),
//  convertions from cartesian to spherical coordinates has been modified.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double st, ct, sp, cp, fi;

  ct = cos( X[ 0 ] );
  st = sin( X[ 0 ] );
  cp = cos( X[ 1 ] );
  sp = sin( X[ 1 ] );


  F[ 0 ] = 0.0;
  for ( int j = 0; j < Na[ 0 ]-1; j++ )
  {
    //std::cout << "** Antenna " << j << ": (" << Xa[ j*N_ANTENNA_DATA     ] << "," << Xa[ j*N_ANTENNA_DATA + 1 ] << "," << Xa[ j*N_ANTENNA_DATA + 2 ] << "): t = " << Xa[ j*N_ANTENNA_DATA + 3  ] << std::endl;
    for ( int i = j+1; i < Na[ 0 ]; i++ )
    {
    //std::cout << "Antenna " << i << ": (" << Xa[ i*N_ANTENNA_DATA     ] << "," << Xa[ i*N_ANTENNA_DATA + 1 ] << "," << Xa[ i*N_ANTENNA_DATA + 2 ] << "): t = " << Xa[ i*N_ANTENNA_DATA + 3  ] << std::endl;
    fi  = -( Xa[ j*N_ANTENNA_DATA     ] - Xa[ i*N_ANTENNA_DATA     ] )*st*sp;
    fi +=  ( Xa[ j*N_ANTENNA_DATA + 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ] )*st*cp;
    fi +=  ( Xa[ j*N_ANTENNA_DATA + 2 ] - Xa[ i*N_ANTENNA_DATA + 2 ] )*ct;
    fi -= X[ 2 ]*( Xa[ j*N_ANTENNA_DATA + 3 ] - Xa[ i*N_ANTENNA_DATA + 3 ] );
    //std::cout << "Chi2 = " << F[ 0 ] << "+ " <<  fi*fi << std::endl;
    F[ 0 ]+= fi*fi;
    }
  }
  //std::cout << "*** Nants = " << Na[ 0 ] << ", theta = " << X[ 0 ]*180/3.14159 << " deg, " << "phi = " << X[ 1 ]*180/3.14159 << " deg ->  Chi2 = " << F[0] << std::endl;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PWF_function_G( int& N, double* X, int& NF, double* G, int* Na,
  double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The gradient of the objective function.
//  OMH 28/09/09: warning: due to modified conventions (x=WE, y=SN),
//  conversions from cartesian to spherical corrdinates has been modified.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double st, ct, sp, cp, xij, yij, zij, tij, fi;

  ct = cos( X[ 0 ] );
  st = sin( X[ 0 ] );
  cp = cos( X[ 1 ] );
  sp = sin( X[ 1 ] );

  G[ 0 ] = G[ 1 ] = G[ 2 ] = 0.0;
  for ( int j = 0; j < Na[ 0 ]-1; j++ )
    for ( int i = j+1; i < Na[ 0 ]; i++ )
  {
    xij = Xa[ j*N_ANTENNA_DATA     ] - Xa[ i*N_ANTENNA_DATA     ];
    yij = Xa[ j*N_ANTENNA_DATA + 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ];
    zij = Xa[ j*N_ANTENNA_DATA + 2 ] - Xa[ i*N_ANTENNA_DATA + 2 ];
    tij = Xa[ j*N_ANTENNA_DATA + 3 ] - Xa[ i*N_ANTENNA_DATA + 3 ];

    fi  = ( -xij*sp + yij*cp )*st + zij*ct - X[2]*tij;

    G[ 0 ]+= ( ( -xij*sp + yij*cp )*ct - zij*st )*fi;
    G[ 1 ]+= ( ( -xij*cp - yij*sp )*st + zij*ct )*fi;
    G[ 2 ]-= tij*fi;
  }
  G[ 0 ]*= 2;
  G[ 1 ]*= 2;
  G[ 2 ]*= 2;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SAM_function_F( int& N, double* X, int& NF, double* F, int* Na,
  double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The objective function to minimise for spherical source amplitude fit.
//  It is built as a chi2 over all antennas, requiring:
//
//  si = s0 /|R-Ri|
//
//  where si is the signal on antenna i, Ri the position of antenna i and
//  R the source position.
//
//  Amplitude is assumed to be expressed in dB.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double xi, yi, zi, fi;

  F[ 0 ] = 0.0;
  for ( int i = 0; i < Na[ 0 ]; i++ )
  {
    xi = X[ 1 ] - Xa[ i*N_ANTENNA_DATA     ];
    yi = X[ 2 ] - Xa[ i*N_ANTENNA_DATA + 1 ];
    zi = X[ 3 ] - Xa[ i*N_ANTENNA_DATA + 2 ];
    fi = Xa[ i*N_ANTENNA_DATA + 4 ] - X[ 0 ] + 10.0*log10( xi*xi + yi*yi +
      zi*zi );

    if ( ( Xa[ i*N_ANTENNA_DATA + 5 ] >= 1 ) && ( fi < 0.0 ) ) // Saturated antenna and fit
      fi = 0.0;

    F[ 0 ]+= fi*fi;
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SAM_function_G( int& N, double* X, int& NF, double* G, int* Na,
  double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The gradient of the objective function.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double fi, di, xi, yi, zi;

  const static double cte = 40.0/log(10.0);

  G[ 0 ] = G[ 1 ] = G[ 2 ] = G[ 3 ] = 0.0;
  for ( int i = 0; i < Na[ 0 ]; i++ )
  {
    xi = X[ 1 ] - Xa[ i*N_ANTENNA_DATA     ];
    yi = X[ 2 ] - Xa[ i*N_ANTENNA_DATA + 1 ];
    zi = X[ 3 ] - Xa[ i*N_ANTENNA_DATA + 2 ];
    di = xi*xi + yi*yi + zi*zi;
    fi = Xa[ i*N_ANTENNA_DATA + 4 ] - X[ 0 ] + 10.0*log10( xi*xi + yi*yi +
      zi*zi );

    if ( ( Xa[ i*N_ANTENNA_DATA + 5 ] >= 1 ) && ( fi < 0.0 ) ) // Saturated antenna and fit
      fi = 0.0;

    G[ 0 ]+= fi;
    G[ 1 ]+= fi*xi/di;
    G[ 2 ]+= fi*yi/di;
    G[ 3 ]+= fi*zi/di;

  }
  G[ 0 ]*= -2.0;
  G[ 1 ]*= cte;
  G[ 2 ]*= cte;
  G[ 3 ]*= cte;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void EAM_function_F( int& N, double* X, int& NF, double* F, int* Na,
  double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The objective function to minimise for exponential amplitude fit. It is
//  built as a chi2 over all antennas, requiring:
//
//  si = s0 - a0*di
//
//  where si is the signal on antenna i and di the lateral distance to the
//  axis given as:
//
//  di^2 = | Xi - X0 |^2 - ( ( Xi - X0 ).u )^2
//
//  with X0 the point on the cascade axis that intercepts the plane z0 = 0,
//  and u the cascade direction
//
//  Amplitude si and loss a0 are assumed to be expressed in dB and dB/m.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double fi, di, xi, yi, zi, ux, uy, uz;

  ux = Xa[ N_ANTENNA_DATA*MAX_ANTENNA     ];
  uy = Xa[ N_ANTENNA_DATA*MAX_ANTENNA + 1 ];
  uz = Xa[ N_ANTENNA_DATA*MAX_ANTENNA + 2 ];

  F[ 0 ] = 0.0;
  for ( int i = 0; i < Na[ 0 ]; i++ )
  {
    xi = X[ 0 ] - Xa[ i*N_ANTENNA_DATA     ];
    yi = X[ 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ];
    zi =        - Xa[ i*N_ANTENNA_DATA + 2 ];
    di = xi*ux + yi*uy + zi*uz;
    di = xi*xi + yi*yi + zi*zi - di*di;
    if ( di < 0.0 ) // prevent rounding error
      di = 0.0;
    else
      di = sqrt( di );
    fi = X[ 3 ]*di + Xa[ i*N_ANTENNA_DATA + 4 ] - X[ 2 ];

    if ( ( Xa[ i*N_ANTENNA_DATA + 5 ] >= 1 ) && ( fi < 0.0 ) ) // Saturated antenna and fit
      fi = 0.0;

    F[ 0 ]+= fi*fi;
  }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void EAM_function_G( int& N, double* X, int& NF, double* G, int* Na,
  double* Xa, void* UFPARM )
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  The gradient of the objective function.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  static double fi, si, di, xi, yi, zi, ux, uy, uz;

  ux = Xa[ N_ANTENNA_DATA*MAX_ANTENNA     ];
  uy = Xa[ N_ANTENNA_DATA*MAX_ANTENNA + 1 ];
  uz = Xa[ N_ANTENNA_DATA*MAX_ANTENNA + 2 ];

  G[ 0 ] = G[ 1 ] = G[ 2 ] = G[ 3 ] = 0.0;
  for ( int i = 0; i < Na[ 0 ]; i++ ) {
    xi = X[ 0 ] - Xa[ i*N_ANTENNA_DATA     ];
    yi = X[ 1 ] - Xa[ i*N_ANTENNA_DATA + 1 ];
    zi =        - Xa[ i*N_ANTENNA_DATA + 2 ];
    si = xi*ux + yi*uy + zi*uz;
    di = xi*xi + yi*yi + zi*zi - si*si;
    if ( di < 0.0 ) // prevent rounding error
      di = 0.0;
    else
      di = sqrt( di );
    fi = X[ 3 ]*di + Xa[ i*N_ANTENNA_DATA + 4 ] - X[ 2 ];

    if ( ( Xa[ i*N_ANTENNA_DATA + 5 ] >= 1 ) && ( fi < 0.0 ) ) // Saturated antenna and fit
      fi = 0.0;

    G[ 0 ]+= fi/di*( xi - si*ux );
    G[ 1 ]+= fi/di*( yi - si*uy );
    G[ 2 ]+= fi;
    G[ 3 ]+= fi*di;
  }
  G[ 0 ]*= 2*X[ 3 ];
  G[ 1 ]*= 2*X[ 3 ];
  G[ 2 ]*= -2;
  G[ 3 ]*= 2;
}
