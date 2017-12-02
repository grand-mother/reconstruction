#include "FitTools.h"

int main()
{
  FitTools* fBox = new FitTools; // Fit toolbox
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // I) ANTENNAE COORDINATES
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  fBox->Na()    = 6;
  fBox->xa( 1 ) = -100.0; fBox->ya( 1 ) =  100.0; fBox->za( 1 ) = -10.0;
  fBox->xa( 2 ) = -100.0; fBox->ya( 2 ) = -100.0; fBox->za( 2 ) = -10.0;
  fBox->xa( 3 ) =  100.0; fBox->ya( 3 ) =  100.0; fBox->za( 3 ) =  10.0;
  fBox->xa( 4 ) =  100.0; fBox->ya( 4 ) = -100.0; fBox->za( 4 ) =  10.0;
  fBox->xa( 5 ) =    0.0; fBox->ya( 5 ) =  -25.0; fBox->za( 5 ) =   0.0;
  fBox->xa( 6 ) =   25.0; fBox->ya( 5 ) =    0.0; fBox->za( 5 ) =   0.0;


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // II) POINT SOURCE AND SPHERICAL AMPLITUDE MODEL
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  double x0, y0, z0, t0, s0, dx, dy, dz, di;
  
  x0 =  15.0;
  y0 =  50.0;
  z0 =   1.0;
  t0 = -20.0;
  s0 =  20.0;
  
  for ( int i = 1; i <=  fBox->Na(); i++ )
  {
    dx = x0 - fBox->xa( i );
    dy = y0 - fBox->ya( i );
    dz = z0 - fBox->za( i );
    di = sqrt( dx*dx + dy*dy +dz*dz );
    fBox->ta( i ) = t0 + di;
    fBox->sa( i ) = s0 - 20.0*log10( di );
  }
  
  double chi2 = -1.0;
  
  fBox->fitModel()    = POINT_SOURCE;
  fBox->fixedSpeed()  = true;
  fBox->cr()          = 1.0;
  chi2                = fBox->scan();
  printf( "1)  POINT SOURCE FIT:\n" );
  printf( "%10.3e %10.3e %10.3e %10.3e [ %10.3e ]\n",  
    fBox->xs(), fBox->ys(), fBox->zs(), fBox->ts(), chi2 );
    
  fBox->fitModel()    = SPHERICAL_AMPLITUDE;
  chi2                = fBox->scan();
  printf( "2a) SPHERICAL AMPLITUDE FIT:\n" );
  printf( "%10.3e %10.3e %10.3e %10.3e [ %10.3e ]\n",  
    fBox->xs(), fBox->ys(), fBox->zs(), fBox->s0(), chi2 );
  
  fBox->fixedZs()     = true;      // Zs is a cte parameter for the fit
  fBox->zs()          = z0 + 10.0; // put some 10 m systematic here
  chi2                = fBox->scan();
  printf( "2b) SPHERICAL AMPLITUDE FIT (Z FIXED):\n" );
  printf( "%10.3e %10.3e %10.3e %10.3e [ %10.3e ]\n", 
    fBox->xs(), fBox->ys(), fBox->zs(), fBox->s0(), chi2 );
    
  fBox->fixedSource() = true;       // Source coorddinates are cte for the fit
  fBox->xs()          = x0 + 10.0;  // put some 10 m systematic on x, y z
  fBox->ys()          = y0 + 10.0;
  fBox->zs()          = z0 + 10.0;
  chi2                = fBox->scan();
  printf( "2c) SPHERICAL AMPLITUDE FIT (X,Y,Z FIXED):\n" );
  printf( "%10.3e %10.3e %10.3e %10.3e [ %10.3e ]\n",  
    fBox->xs(), fBox->ys(), fBox->zs(), fBox->s0(), chi2 );
    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // III) PLAN WAVE
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  double theta, phi, a0, ct, st, cp, sp;
  
  theta = 45.0/180.0*M_PI;
  phi   = 30.0/180.0*M_PI;
  
  ct = cos( theta );
  st = sin( theta );
  cp = cos( phi );
  sp = sin( phi );
  
  fBox->sigma_t()     = 3.0; // error on times, in m.
  fBox->fitModel()    = PLAN_WAVE;
  fBox->fixedSpeed()  = true;
  fBox->cr()          = 1.0;
  
  // Test p-value by toy MC
  //===
  FILE* fid = fopen( "gchi2-test.txt", "w+" );
  for ( int j = 0; j < 1000; j++ )
  {
    for ( int i = 1; i <=  fBox->Na(); i++ )
      fBox->ta( i ) = st*( fBox->ya( i )*cp - fBox->xa( i )*sp ) + ct*fBox->za( i ) +
        fBox->randn()*fBox->sigma_t(); // Arrival times with gaussian errors 
	
    chi2 = fBox->scan();
    
    fprintf( fid, "%13.7e %13.7e\n", chi2, fBox->chi2_significance() );
  }
  fclose( fid );
  
  // Reconstruction with error propagation computation
  //===
  for ( int i = 1; i <=  fBox->Na(); i++ )
    fBox->ta( i ) = st*( fBox->ya( i )*cp - fBox->xa( i )*sp ) + ct*fBox->za( i ); // No errors
    
  chi2 = fBox->scan();
  
  printf( "3)  PLAN WAVE FIT:\n" );
  printf( "%10.3e +- %10.3e %10.3e +- %10.3e [ %10.3e ]\n",  
    fBox->theta()*180.0/M_PI, fBox->theta_error()*180.0/M_PI, 
    fBox->phi()*180.0/M_PI,  fBox->phi_error()*180.0/M_PI,
    fBox->chi2_significance() );
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // IV) Generalised chi2 cdf
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  double chi2obs, sig, lambda[ 5 ], delta[ 5 ];
  int ndof[ 5 ];
  int n;
  
  // Benchmark for algorithms comparison
  //===
  lambda[ 0 ] = 1.0; delta[ 0 ] = 0.0; ndof[ 0 ] = 1;
  lambda[ 1 ] = 3.0; delta[ 1 ] = 0.0; ndof[ 1 ] = 1;
  lambda[ 2 ] = 5.0; delta[ 2 ] = 0.0; ndof[ 2 ] = 1;
  lambda[ 3 ] = 7.0; delta[ 3 ] = 0.0; ndof[ 3 ] = 1;
  lambda[ 4 ] = 9.0; delta[ 4 ] = 0.0; ndof[ 4 ] = 1;
  n = 5;
  chi2obs = 20.0;
  
  int algo[] = { 153, 155, 204, 256 };
  printf( "4) GENERALISED CHI2 ALGORITHMS:\n" );
  printf( " %3s %10s %10s\n", "ALG", "G-CHI2", "CDF" );
  for ( int i = 0; i < 4; i++ )
  {
    fBox->gchi2_algorithm() = algo[ i ];
    sig = fBox->gchi2cdf( lambda, delta, ndof, n, chi2obs );
    printf( " %3d %10.3e %10.3e\n", algo[ i ], chi2obs, sig );
  }
  
  // For timing
  //===
  fBox->gchi2_algorithm() = 204;
  for ( int i = 0; i < 1000000; i++ )
    sig = fBox->gchi2cdf( lambda, delta, ndof, n, chi2obs );
  
  //  Time in micro-s for the differents algorithms:
  //
  //  AS 153     7.5 **
  //  AS 155   110.0 *
  //  AS 204     1.5 ****
  //  AS 256   405.0 *
  //
  //  ==> default algorithm is AS 204 (based on Ruben's expansion  
  //  in a series of chi2's). Note that algorithm AS 153 requires   
  //  all lambda_i's to be distinct and ndof's to be all equal to 1. 
  //  Hence it has a limited usage case.
    
  delete fBox;

  return( 0 );
}
