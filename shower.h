#ifndef _SHOWER_H_
#define _SHOWER_H_

#include <vector>
#include <cstring>

#include "TRandom3.h"
#include "TGraph.h"

using namespace std;


//! Shower generator
class TShower
{
  public :
    //! Constructor
    TShower(double energy, double * coord, unsigned int step = 800);

    //! Destructor
    ~TShower();

    //! Generates shower
    void GenerateShower();

    //! Returns #fEnergy
    double GetEnergy() const {return fEnergy;}

    //! Returns the zenith and azimuth angle of the incoming cosmic ray
    void GetIncomingDirection(double & theta, double & phi) const {theta = fTheta, phi = fPhi;}

    //! Get #fT1
    double GetT1() const {return fT1;}
    
    //! Get #fTmax
    double GetTmax() const {return fTmax;}
    
    //! Get #fTmax
    unsigned int GetStep() const {return fStep;}
    
    //! Get #fTmax
    bool GetStatus() const {return fStatus;}
    
    //! Get the longitudinal profile of the EAS
    void GetLongitudinalProfile(vector<double> & T, vector<double> & Ne); 

    //! Get the longitudinal profile of the EAS
    TGraph * GetLongitudinalProfile();
  
  private :
    //! Initializes #fRandom, #fX, #fX1, #fT and #fT1
    void Init();

    //! Random generator
    TRandom3 * fRandom;

    //! Energy in eV
    double fEnergy;

    //! Zenith angle
    double fTheta;

    //! Azimuth angle
    double fPhi;

    //! Step in radiation length
    unsigned int fStep;

    //! Depth of the first interaction in unit of radiation length
    double fT1;

    //! Depth at shower maximum in unit of radiation length
    double fTmax;
    
    //! Number of radiation length
    vector<double> fT;

    //! Number of electrons/positrons
    vector<double> fNe;

    //! Tells you if the shower has been generated or not
    bool fStatus;
};


//! Mean longitudinal development of the electron/positron component of photon initiated electromagnetic EAS
//! Greisen (1956)
vector<double> Greisen(vector<double> & T, double energy);

//! Mean longitudinal development of the electron/positron component of photon initiated electromagnetic EAS
//! Greisen (1956)
double Greisen(double T, double energy);

//! Electron energy spectrum between 1 MeV and 10 GeV in MeV
//! Nerling et al. (2006)
vector<double> ElectronEnergySpectrum(vector<double> & energy, double age);



#endif
