#include "conversion.h"

#include <iostream>
#include <cmath>



vector<double> age2depth(const vector<double> & age, double Xmax)
{
  unsigned int size = age.size();
  vector<double> X(size);
  for(unsigned int i = 0; i < size; i++) X[i] = 2.*Xmax*1./(3./age[i]-1.);

  return X; 
}



double age2depth(double age, double Xmax)
{
  double X = 2.*Xmax*1./(3./age-1.);
  
  return X; 
}



vector<double> depth2age(const vector<double> & X, double Xmax)
{
  unsigned int size = X.size();
  vector<double> age(size);
  for(unsigned int i = 0; i < size; i++) age[i] = 3./(1.+2.*Xmax/X[i]);

  return age;
}



double depth2age(double X, double Xmax)
{
  double age = 3./(1.+2.*Xmax/X);

  return age;
}



double altitude2depth(double altitude)
{
  // Following the functional form used in CORSIKA, the depth profile is divided into four layers
  double depth = 0.;

  if( altitude < -5.801 ) {cout << "ERROR: altitude lower than -5.801 km. EXITING." << endl; exit(0);}

  if( altitude > 112.8 ) return depth;

  // atmospheric depth decreases linearly with altitude
  if( altitude >= 100.0 && altitude <= 112.8 ) {depth = 1.128292e-2-altitude*1./1.e4; return depth;} 

  double a[4] = {-1.865562e2 , -9.49199e1 , 6.1289e-1 , 0.0};
  double b[4] = {1.2226562e3 , 1.1449069e3 , 1.3055948e3 , 5.401778e2};
  double c[4] = {9.9418638 , 8.7815355 , 6.3614304 , 7.7217016};

  unsigned int par = 0;
  if( altitude >= 40.0 && altitude < 100.0 ) par = 3;
  if( altitude >= 10.0 && altitude < 40.0 ) par = 2;
  if( altitude >= 4.0 && altitude < 10.0 ) par = 1;
  if( altitude >= -5.801 && altitude < 4.0 ) par = 0;

  depth = a[par]+b[par]*exp(-altitude/c[par]);

  return depth;
}



double depth2altitude(double depth)
{
  double altitude = 0.;
  if( depth < 0.0 ) { cout << "ERROR: atmospheric depth is negative. EXITING." << endl; exit(0);}

  // atmospheric depth decreases linearly with altitude
  if( depth <= 0.0012829199 ) {altitude = 1.e4*(1.128292e-2-depth); return altitude;} 

  double a[4] = {-1.865562e2 , -9.49199e1 , 6.1289e-1 , 0.0};
  double b[4] = {1.2226562e3 , 1.1449069e3 , 1.3055948e3 , 5.401778e2};
  double c[4] = {9.9418638 , 8.7815355 , 6.3614304 , 7.7217016};

  unsigned int par = 0;
  if( depth <= 3.03950 && depth > 0.0012829199 ) par = 3;
  if( depth <= 271.700 && depth > 3.03950 ) par = 2;
  if( depth <= 631.100 && depth > 271.700 ) par = 1;
  if( depth <= 2004.79 && depth > 631.100 ) par = 0;

  altitude = -c[par]*log( (depth-a[par]) / b[par]);

  return altitude;
}





