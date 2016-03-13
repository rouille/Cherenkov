#ifndef _COMMON_H
#define _COMMON_H

#include <string>
#include <vector>

#include <TCanvas.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TGraph.h>

using namespace std;

/*! This is stricly identical to :
  \code
  vector<Type> NewObjectName(OriginalObject.size());
  for (unsigned int i = 0; i < OriginalObject.size(); i++) NewObjectName[i] = OriginalObject[i].FieldName;
  \encode
 */
#define DECLARE_VECTOR(Type,NewObjectName,OriginalObject,FieldName) std::vector<Type> NewObjectName;{size_t i=0;size_t size=OriginalObject.size();NewObjectName.reserve(size); while( i<size ){NewObjectName.push_back(OriginalObject[i].FieldName); i++;}}

/*! This is stricly identical to :
  \code
  Type * NewObjectName = new Type[OriginalObject.size()];
  for (unsigned int i = 0; i < OriginalObject.size(); i++) NewObjectName[i] = OriginalObject[i].FieldName;
  \endcode
 */
#define DECLARE_POINTER(Type,NewObjectName,OriginalObject,FieldName) Type * NewObjectName=new Type[OriginalObject.size()]; {size_t i=0;size_t size=OriginalObject.size(); while( i<size ){NewObjectName[i]=OriginalObject[i].FieldName;i++;}}

//! Defines some mathematical constants
namespace kMathConstants
{
  //! Radian to degree conversion
  const double RTOD =  57.2957795130823;

  //! Degree to radian conversion
  const double DTOR =  0.017453292519943295;

  //! \f$ 2\pi \f$
  const double TwoPi = 6.28318530717958623200;

  //! \f$ \frac{\pi}{2} \f$
  const double PiOver2 = 1.57079632679489661923;
}

//! Defines some physical constants
namespace kPhysicalConstants
{
  //! Speed of light in \f$ m . s^{-1} \f$
  const double CSPEED =  299792458;

  //! Mass of the electron in MeV
  const double Me = 0.511;

  //! Fine structure constant
  const double alpha = 0.007297352569;

  //! Crtical energy in air in eV 
  const double Ec = 86.e6;

  //! Radiation length in air in \f$ g . cm^{-2} \f$
 const double X0 = 37.1;

  //! Mean interaction of photons in unit of radiation length
  const double Tint = 9./7.; 
}

//! Tells you whether the file is available or not.
bool CheckFile(string fileName);

//! Plots a 2D graph
void SetPlotAttributes(TGraphErrors *, double xMin, double xMax, string name = "", string Xaxis = "X", string Yaxis = "Y");

//! Plots an histogram
void SetPlotAttributes(TH1F * Histo, string name, string Xaxis);

//! Binning
vector<double> Bins(unsigned int size, double min, double max, bool logarithmic = false);

/*! 
  This is a simple 5 points Newton-Cotes formula that allows numerical integration of the function given by x 
  and y. This (simple) algorithm requires the x values to be equally spaced and in increasing order. It is very 
  accurate as long as the function is well sampled (see Numerical Recipes in C for details).
 */
double Integrate_nc5(const vector<double> & x, const vector<double> & y); 

/* 
   Given the vectors x and y, wich tabulate a function (with the x's in order), this routine returns a linear
   interpolation at point u
 */
vector<double> Interpol(const vector<double> & x, const vector<double> & y, const vector<double> & u);

/* 
   Given the vectors x and y, wich tabulate a function (with the x's in order), this routine returns a linear
   interpolation at point u
 */
double Interpol(const vector<double> & x, const vector<double> & y, double u);


#endif
