#ifndef _CONVERSION_H
#define _CONVERSION_H

#include <vector>

using namespace std;

//! Shower age to slant depth [\f$ g . cm^{-2} \f$] 
vector<double> age2depth(const vector<double> & age, double Xmax);

//! Shower age to slant depth [\f$ g . cm^{-2} \f$] 
double age2depth(double age, double Xmax);

//! Slant depth [\f$ g . cm^{-2} \f$] to shower age
vector<double> depth2age(const vector<double> & X, double Xmax);

//! Slant depth [\f$ g . cm^{-2} \f$] to shower age 
double depth2age(double X, double Xmax);

//! Altitude [km] to depth [\f$ g . cm^{-2} \f$]
double altitude2depth(double altitude);

//! Depth [\f$ g . cm^{-2} \f$] to altitude [km]
double depth2altitude(double depth);


#endif
