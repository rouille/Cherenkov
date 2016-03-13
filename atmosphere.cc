#include "atmosphere.h"

#include <iostream>
#include <fstream>

vector<TAtmosphere> GetAtmosphere(string fileName)
{
  vector<TAtmosphere> Atmosphere;
  ifstream atmosphereFile(fileName.c_str());
  double altitude, density, depth, delta;
  while( atmosphereFile >> altitude ) 
  {
    atmosphereFile >> density >> depth >> delta;
    TAtmosphere AtmosphereTmp;
    AtmosphereTmp.fAltitude = altitude;
    AtmosphereTmp.fDensity = density;
    AtmosphereTmp.fDepth = depth;
    AtmosphereTmp.fDelta = delta;
    AtmosphereTmp.fRefractiveIndex = 1+delta;

    Atmosphere.push_back(AtmosphereTmp);
  }

  atmosphereFile.close();

  return Atmosphere;
}
