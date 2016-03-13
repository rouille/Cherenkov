#include "shower.h"
#include "common.h"
#include "conversion.h"

#include <sys/time.h>
#include <cmath>
#include <iostream>



using namespace kPhysicalConstants;


TShower::TShower(double energy, double * coord, unsigned int step)
{
  fEnergy = energy;
  fTheta = coord[0];
  fPhi = coord[1];
  fStep = step;
  fStatus = false;

  Init();
}



TShower::~TShower()
{
  if( fRandom ) delete fRandom;
}



void TShower::Init()
{
  // Random generator
  unsigned int seed;
  struct timeval MyTimeVal;
  struct timezone MyTimeZone;
  gettimeofday(&MyTimeVal,&MyTimeZone);
  seed = (unsigned int) (MyTimeVal.tv_usec+(MyTimeVal.tv_sec % 1000)*1000000);
  fRandom = new TRandom3(seed);

  // Depth of the first interaction
  fT1 = -Tint*log(fRandom->Rndm());

  // Number of radiation length
  fT = Bins(fStep,0.1,40);
}


void TShower::GenerateShower()
{
  // I Y Crewther and R J Protheroe (1990)

  // Depth of the shower maximum in unit of radiation length
  double y = log(fEnergy/Ec);

  // Number of radiation length measured from first interaction
  vector<double> Tprime(fStep);
  for(unsigned int i = 0; i < fStep; i++) Tprime[i] = fT[i]-fT1;

  // Fluctuations
  double Sprime = 0, F = 0, N1 = 0, Sigma = 0, Mu = 0;
  vector<double> LogN(fStep);
  double RanNormal = fRandom->Gaus();

  for(unsigned int i = 0; i < fStep; i++)
    {
      if( Tprime[i] >= 0. )
        {
          Sprime = depth2age(Tprime[i],y);
          F = (0.88+0.146*Sprime)*(1-exp(-3.84*Sprime));
          N1 = (y/Tprime[i])*Greisen(Tprime[i],fEnergy)*F;
          Sigma = 0.157-0.0048*y+2.34*(Sprime-1)*(Sprime-1);
          Mu = log(N1)-(Sigma*Sigma)/2.;
          LogN[i] = Mu+Sigma*RanNormal;
        }
    }

  fNe.resize(fStep);
  for(unsigned int i = 0; i < fStep; i++) {if( Tprime[i] <= 0 ) fNe[i] = 0; else fNe[i] = exp(LogN[i]);}

  // Depth at shower maximum in unit of radiation length
  unsigned int index_max = 0;
  for(unsigned int i = 1; i < fStep; i++) if( fNe[i] > fNe[i-1] ) index_max = i;
  fTmax = fT[index_max];  

  // Status
  fStatus = true;
}



void TShower::GetLongitudinalProfile(vector<double> & T, vector<double> & Ne)
{
  if( fStatus == false ) {cout << "Call TShower::GenerateShower first. EXITING." << endl; exit(0);}
  T.resize(fStep);
  Ne.resize(fStep);
  for(unsigned int i = 0; i < fStep; i++) {T[i] = fT[i]; Ne[i] = fNe[i];}
}



TGraph * TShower::GetLongitudinalProfile()
{
  if( fStatus == false ) {cout << "Call TShower::GenerateShower first. EXITING." << endl; exit(0);}
  TGraph * gProfile = new TGraph(fStep);
  for(unsigned int i = 0; i < fStep; i++) gProfile->SetPoint(i,fT[i],fNe[i]);

  return gProfile;
}



vector<double> Greisen(vector<double> & T, double energy)
{
  // Mean depth of shower maximum in unit of radiation length
  double y = log(energy/Ec);

  // Shower age
  vector<double> age = depth2age(T,y);

  // Number of electrons/positrons as a function of T
  vector<double> Ne(T.size());
  for(unsigned int i = 0; i < T.size(); i++) Ne[i] = (0.31/sqrt(y))*exp(T[i]*(1-1.5*log(age[i])));

  return Ne;
}



double Greisen(double T, double energy)
{
  // Mean depth of shower maximum in unit of radiation length
  double y = log(energy/Ec);

  // Shower age
  double age = depth2age(T,y);

  // Number of electrons/positrons as a function of T
  double Ne = (0.31/sqrt(y))*exp(T*(1-1.5*log(age)));

  return Ne;
}



vector<double> ElectronEnergySpectrum(vector<double> & energy, double age)
{
  // valid for electrons with energy > 1 MeV
  unsigned int size = energy.size();

  double k0 = 0.145098;
  double k1 = 6.20114;
  double k2 = -0.596851;

  double a1 = 6.42522-1.53183*age;
  double a2 = 168.168-42.1368*age;
  double a0 = k0*exp(k1*age+k2*age*age);

  vector<double> spectrum(size);
  for(unsigned int i = 0; i < size; i++) spectrum[i] = a0*energy[i]/((energy[i]+a1)*pow(energy[i]+a2,age));

  vector<double> LogEnergy(size);
  for(unsigned int i = 0; i < size; i++) LogEnergy[i] = log(energy[i]);

  return spectrum;
}
