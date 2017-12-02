#include "FitTools.h"

int main()
{
  FitTools* fBox = new FitTools; // Fit toolbox
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // I) ANTENNAE COORDINATES
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  fBox->Na()    = 4;
  fBox->xa( 1 ) = 35.0; fBox->ya( 1 ) =  462.0; fBox->za( 1 ) =   2642.0;  //110
  fBox->xa( 2 ) = -23.0; fBox->ya( 2 ) =  399.0; fBox->za( 2 ) =   2642.0;  //113
  fBox->xa( 3 ) = 92.0; fBox->ya( 3 ) =  354.0; fBox->za( 3 ) =   2638.0;  //114
  fBox->xa( 4 ) = -25.0; fBox->ya( 4 ) =  269.0; fBox->za( 4 ) =   2634.0;  //116
  
  //fBox->Na()    = 6;
  //fBox->xa( 1 ) = 2178.0; fBox->ya( 1 ) =  -25.0; fBox->za( 1 ) =   0.0;  //110
  //fBox->xa( 2 ) = 2225.0; fBox->ya( 2 ) =  -52.0; fBox->za( 2 ) =   1.0;  //113
  //fBox->xa( 3 ) = 2303.0; fBox->ya( 3 ) =  -23.0; fBox->za( 3 ) =   3.0;  //114
  //fBox->xa( 4 ) = 2380.0; fBox->ya( 4 ) =  -47.0; fBox->za( 4 ) =   7.0;  //116
  //fBox->xa( 5 ) = 2420.0; fBox->ya( 5 ) =  -21.0; fBox->za( 5 ) =   7.0;  //117
  //fBox->xa( 6 ) = 2332.0; fBox->ya( 6 ) = -115.0; fBox->za( 6 ) =   6.0;  //115

  
 
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // II) PLAN WAVE
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  double theta, phi, dt, dp, ct, st, cp, sp, chi2;
  FILE* fid;
  
  dt = 1.0/180*M_PI;  
  dp = 1.0/180*M_PI;
  
  fBox->sigma_t()     = 3.0; // error on times, in m.
  fBox->fitModel()    = PLAN_WAVE;
  fBox->fixedSpeed()  = true;
  fBox->cr()          = 1.0;
  
  // Error propagation as theta
  //===
  phi = 195.4/180.0*M_PI;

  cp  = cos( phi );
  sp  = sin( phi );
  fid = fopen( "proto6a-theta.txt", "w+" );
  for ( int it = 51; it <= 52; it++ )
  {
    theta = 0.0 + dt*it;
    ct    = cos( theta );
    st    = sin( theta );
    for ( int i = 1; i <=  fBox->Na(); i++ )
      fBox->ta( i ) = st*( fBox->ya( i )*cp - fBox->xa( i )*sp ) + ct*fBox->za( i );
    
    chi2 = fBox->scan();
    
    fprintf( fid, "%10.3e %10.3e %10.3e %10.3e\n",  
      fBox->theta()*180.0/M_PI, fBox->theta_error()*180.0/M_PI, 
      fBox->phi()*180.0/M_PI,  fBox->phi_error()*180.0/M_PI );
  }
  fclose( fid );
  
  // Error propagation as phi
  //===
  theta = 45.0/180.0*M_PI;
  ct    = cos( theta );
  st    = sin( theta );
  fid   = fopen( "proto6a-phi.txt", "w+" );
  for ( int ip = 0; ip <= 360; ip++ )
  {
    phi = 0.0 + dp*ip;
    cp  = cos( phi );
    sp  = sin( phi );
    for ( int i = 1; i <=  fBox->Na(); i++ )
      fBox->ta( i ) = st*( fBox->ya( i )*cp - fBox->xa( i )*sp ) + ct*fBox->za( i );
    
    chi2 = fBox->scan();
    
    fprintf( fid, "%10.3e %10.3e %10.3e %10.3e\n",  
      fBox->theta()*180.0/M_PI, fBox->theta_error()*180.0/M_PI, 
      fBox->phi()*180.0/M_PI,  fBox->phi_error()*180.0/M_PI );
  }
  fclose( fid );
  
  // Error propagation 2D map
  //===
  fid   = fopen( "proto6a-2D.txt", "w+" );
  for ( int it = 0; it <= 90; it++ )
  {
    theta = 0.0 + dt*it;
    ct    = cos( theta );
    st    = sin( theta );
    for ( int ip = 0; ip <= 360; ip++ )
    {
      phi = 0.0 + dp*ip;
      cp  = cos( phi );
      sp  = sin( phi );
      for ( int i = 1; i <=  fBox->Na(); i++ )
        fBox->ta( i ) = st*( fBox->ya( i )*cp - fBox->xa( i )*sp ) + ct*fBox->za( i );
    
      chi2 = fBox->scan();
    
      fprintf( fid, " %10.3e %10.3e",  
        fBox->theta_error()*180.0/M_PI, 
        fBox->phi_error()*180.0/M_PI );
    }
    fprintf( fid, "\n" );
  }
  fclose( fid );
    
  delete fBox;

  return( 0 );
}
