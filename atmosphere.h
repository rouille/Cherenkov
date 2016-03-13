#ifndef _ATMOSPHERE_H_
#define _ATMOSPHERE_H_

#include <vector>
#include <cstring>

using namespace std;

//! Atmospheric model
class TAtmosphere
{
  public :
    //! Constructor
    TAtmosphere() {}

    //! Altitude in km
    double fAltitude;

    //! Density in \f$ g . cm^{-3} \f$
    double fDensity;

    //! Atmospheric depth in \f$ g . cm^{-2} \f$ 
    double fDepth;

    //! Refractive index (independent of the wavelength)
    double fRefractiveIndex;

    //! Refractive index - 1
    double fDelta;
};

vector<TAtmosphere> GetAtmosphere(string fileName);

#endif
