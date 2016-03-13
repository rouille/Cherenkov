#include <cmath>
#include <iostream>

#include "conversion.h"
#include "cherenkov.h"
#include "common.h"

using namespace std;
using namespace kMathConstants;
using namespace kPhysicalConstants;



TCherenkov::TCherenkov(const vector<TAtmosphere> & atmosphere, TShower * shower, double waveMin, double waveMax)
{
  fAtmosphere = atmosphere;
  fShower = shower;
  fWaveMin = waveMin; // in cm
  fWaveMax = waveMax; // in cm

  if( !fShower->GetStatus() ) {cout << "Generate shower first. EXITING." << endl; exit(0);}
}



TCherenkov::~TCherenkov()
{
  if( fShower ) delete fShower;
}



void TCherenkov::ComputeTotalNumberPhotons(vector<double> & T, vector<double> & Nc)
{
  /* Shower */
  // Longitudinal development
  unsigned int size_shower = fShower->GetStep();
  vector<double> Ne;
  fShower->GetLongitudinalProfile(T,Ne);
  // Depth at maximum development
  double Tmax = fShower->GetTmax();
  // Incoming direction
  double theta, phi;
  fShower->GetIncomingDirection(theta,phi);

  /* Slant depth to age */
  vector<double> age = depth2age(T,Tmax);

  /* Slant depth to altitude */
  vector<double> altitude(size_shower);
  for(unsigned int i = 0; i < size_shower; i++) altitude[i] = depth2altitude(T[i]*X0*cos(theta*DTOR));

  /* Linear interpolation of density and delta at altitude */
  DECLARE_VECTOR(double,altitude_table,fAtmosphere,fAltitude);
  DECLARE_VECTOR(double,density_table,fAtmosphere,fDensity);
  DECLARE_VECTOR(double,delta_table,fAtmosphere,fDelta);
  vector<double> density = Interpol(altitude_table,density_table,altitude);
  vector<double> delta = Interpol(altitude_table,delta_table,altitude);

  /* Total number of produced Cherenkov photons */
  Nc.resize(size_shower);
  unsigned int size_spectrum = 100;
  vector<double> Ee = Bins(size_spectrum,1,10000,true); // Electrons energy between 1 MeV and 10 GeV
  for(unsigned int i = 0; i < size_shower; i++)
    {
      // Normalized differential electron energy spectrum at altitude 
      vector<double> Se = ElectronEnergySpectrum(Ee,age[i]);

      // Normalized total number of Cherenkov photons produced
      vector<double> LogEe, Sc; // Electron energy above threshold and corresponding number of Cherenkov photons produced
      double yield = 0.; // Cherenkov yield
      for(unsigned int j = 0; j < size_spectrum; j++)
        {
          if( Ee[j] > EnergyThreshold(delta[i])) // Cherenkov condition
            {
              yield = Yield(Ee[j],delta[i],density[i]);
              LogEe.push_back(log(Ee[j]));
              Sc.push_back(Se[j]*yield);
            }
        }
      double NormalizedNc = Integrate_nc5(LogEe,Sc);

      // Total number of Cherenkov photons produced
      Nc[i] = Ne[i]*NormalizedNc;
    }
}


void TCherenkov::ComputeAngularDistribution(vector<double> & T, vector<double> & angle, vector<vector<double> > & distribution)
{
  /* Shower */
  // Longitudinal development
  unsigned int size_shower = fShower->GetStep();
  vector<double> Ne;
  fShower->GetLongitudinalProfile(T,Ne);
  // Depth at maximum development
  double Tmax = fShower->GetTmax();
  // Incoming direction
  double theta, phi;
  fShower->GetIncomingDirection(theta,phi);

  /* Slant depth to age */
  vector<double> age = depth2age(T,Tmax); 

  /* Slant depth to altitude */
  vector<double> altitude(size_shower);
  for(unsigned int i = 0; i < size_shower; i++) altitude[i] = depth2altitude(T[i]*X0*cos(theta*DTOR));

  /* Linear interpolation of density and delta at altitude */
  DECLARE_VECTOR(double,altitude_table,fAtmosphere,fAltitude);
  DECLARE_VECTOR(double,delta_table,fAtmosphere,fDelta);
  vector<double> delta = Interpol(altitude_table,delta_table,altitude);

  // Normalized angular distribution
  unsigned int size_angle = 180;
  angle = Bins(size_angle,0.,180.);
  distribution.resize(size_shower);
  for(unsigned int i = 0; i < size_shower; i++)
    {
      vector<double> distribution_tmp = AngularDistribution(angle,age[i],delta[i]);
      distribution[i].resize(size_angle);
      for(unsigned int j = 0; j < size_angle; j++) distribution[i][j] = distribution_tmp[j];
    }
}



double TCherenkov::EnergyThreshold(double delta)
{
  double Eth = Me/sqrt(2*delta);

  return Eth;
}



double TCherenkov::Yield(double energy, double delta, double density)
{
  double yield = 0.;

  // Below the Cherenkov energy threshold
  if( energy < EnergyThreshold(delta) ) return yield;

  // Above the Cherenkov energy threshold. See Eq. 2 in Nerling et al. (2006)
  unsigned int size = 1000;
  vector<double> wave = Bins(size, fWaveMin, fWaveMax);
  vector<double> integrand(size);
  for(unsigned int i = 0; i < size; i++) integrand[i] = (2.*delta-pow(Me,2)/pow(energy,2))/pow(wave[i],2);
  yield = (TwoPi*alpha/density)*Integrate_nc5(wave,integrand); 

  return yield;
}



vector<double> TCherenkov::AngularDistribution(vector<double> & angle, double age, double delta)
{
  unsigned int size = angle.size();
  vector<double> angle_rad(size);
  for(unsigned int i = 0; i < size; i++) angle_rad[i] = angle[i]*DTOR; 

  // Parametrization from Neirling et al. (2006)
  double a0 = 0.42489, a1 = 0.58371, a2 = -0.082373;
  double a = a0+a1*age+a2*pow(age,2);

  double b0 = 0.055108, b1 = -0.095587, b2 = 0.056952;
  double b = b0+b1*age+b2*pow(age,2);

  double Eth = EnergyThreshold(delta);
  double theta_c = 0.62694*pow(Eth,-0.60590);
  double theta_cc = (10.509-4.9644*age)*theta_c;

  vector<double> distribution(size);
  for(unsigned int i = 0; i < size; i++) distribution[i] = a*(1./theta_c)*exp(-angle_rad[i]/theta_c)+b*(1./theta_cc)*exp(-angle_rad[i]/theta_cc);

  // Normalization
  double norm = Integrate_nc5(angle_rad,distribution);
  for(unsigned int i = 0; i < size; i++) distribution[i] = distribution[i]*DTOR/norm;

  return distribution;
}



double * CherenkovEnergyThreshold(unsigned int size, const double * delta)
{
  double * Eth = new double[size];
  for(unsigned int i = 0; i < size; i++) Eth[i] = Me/sqrt(2*delta[i]);

  return Eth;
}



