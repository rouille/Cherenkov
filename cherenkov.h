#ifndef _CHERENKOV_H
#define _CHERENKOV_H

#include "atmosphere.h"
#include "shower.h"

#include <vector>

using namespace std;



class TCherenkov
{
  public :
    //! Constructor
    TCherenkov(const vector<TAtmosphere> & atmosphere, TShower * shower, double waveMin, double waveMax);

    //! Destructor
    ~TCherenkov();

    //! Total number of Cherenkov photons produced
    void ComputeTotalNumberPhotons(vector<double> & T, vector<double> & Nc);
    
    //! Normalized angular distribution with respect to shower axis
    void ComputeAngularDistribution(vector<double> & T, vector<double> & angle, vector<vector<double> > & distribution);

  private :
    //! Atmosphere
    vector<TAtmosphere> fAtmosphere;

    //! Shower
    TShower * fShower;

    //! Minimum wavelength of Cherenkov photons produced
    double fWaveMin;
    
    //! Maximum wavelength of Cherenkov photons produced
    double fWaveMax;
    
    //! Energy threshold condition for Cherenkov in air (in MeV)
    double EnergyThreshold(double delta);

    //! Number of Cherenkov photons produced by a electron/positron 
    double Yield(double energy, double delta, double density);

    //! Normalized angular distribution of produced Cherenkov photons
    vector<double> AngularDistribution(vector<double> & angle, double age, double delta);
};

//! Energy threshold condition for Cherenkov radiation in air (in MeV)
double * CherenkovEnergyThreshold(unsigned int size, const double * delta);

#endif
