#include <iostream>
#include <vector>
#include <cmath>

#include "atmosphere.h"
#include "common.h"
#include "cherenkov.h"
#include "shower.h"
#include "conversion.h"

#include <TObject.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TRint.h>
#include <TRoot.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TMath.h>

static char gObjName[1024];
static int gObjNumber = 0;
static char * GetObjName() {sprintf(gObjName,"Cherenkov %d",gObjNumber++); return gObjName;}


using namespace kPhysicalConstants;
using namespace kMathConstants;
using namespace std;



void Usage(string myName)
{
  cout << endl;
  cout << " Synopsis : " << endl;
  cout << myName << " <atmospheric file> <log(energy/[eV])> <zenith angle> <azimuth angle>" << endl << endl;

  cout << " Description :" << endl;
  cout << myName << " generates a shower produced by a cosmic ray of energy <log(energy/[eV])> and incoming"
       << " direction <zenith angle> and <azimuth angle> and plots the total number of Cherenkov photons"
       << " produced and the normalize angular distribution with respect to the shower axis." << endl;

  cout << endl;
  exit(0);
}



int main(int argc, char* argv[])
{
  // Command line
  if(argc != 5) Usage(argv[0]);
  string AtmosphereFile = argv[1];
  if( !CheckFile(AtmosphereFile) ) {cerr << "Exiting" << endl; exit(0);}
  double energy = pow(10,atof(argv[2]));
  double coord[2] = {atof(argv[3]),atof(argv[4])};

  // ROOT
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFont(30,"TITLE");

  // Atmosphere
  vector<TAtmosphere> atmosphere = GetAtmosphere(AtmosphereFile);

  // Wavelength range for Cherenkov photons produced (in cm)
  double WaveMin = 300e-7, WaveMax = 400e-7;

  /* Let's go */
  unsigned int NumberOfShower = 20;

  vector<vector<double> > Nc(NumberOfShower);
  vector<double> angle;
  vector<vector<double> > AngularDistribution;
  vector<double> T;
  double Tmax = 0;

  for(unsigned int i = 0; i < NumberOfShower; i++)
    {
      cout << "Shower #" << i+1 << "/" << NumberOfShower << endl;
      TShower * Shower = new TShower(energy,coord,50);
      Shower->GenerateShower();
      TCherenkov * Cherenkov = new TCherenkov(atmosphere,Shower,WaveMin,WaveMax);
      if( i == 0 ) {Tmax = Shower->GetTmax(); Cherenkov->ComputeAngularDistribution(T,angle,AngularDistribution);}
      vector<double> T_tmp, Nc_tmp;
      Cherenkov->ComputeTotalNumberPhotons(T_tmp,Nc_tmp);
      Nc[i].resize(T_tmp.size());
      for(unsigned j = 0; j < T_tmp.size(); j++) Nc[i][j] = Nc_tmp[j];
    }

  // Total number of Cherenkov photons produced
  TCanvas * cCherenkov = new TCanvas(GetObjName(),"Cherenkov",600,600);
  for(unsigned int i = 0; i < NumberOfShower; i++)
    {
      TGraphErrors * gCherenkov = new TGraphErrors(T.size());
      for(unsigned int j = 0; j < T.size(); j++) gCherenkov->SetPoint(j,T[j]*X0,Nc[i][j]);
      if( i == 0 )
        {
          TGaxis::SetMaxDigits(4);
          SetPlotAttributes(gCherenkov,0,1400,"","X [g . cm^{-2}]","dN_{#gamma}/dX per g . cm^{-2}");
          gCherenkov->GetYaxis()->SetRangeUser(0.,1.5*TMath::MaxElement(T.size(),gCherenkov->GetY()));
          gCherenkov->SetLineWidth(2);
          gCherenkov->Draw("AC");
        }
      else
        {
          gCherenkov->SetLineWidth(2);
          gCherenkov->Draw("same C");
        }
      cCherenkov->Update();
    }

  char LegendEnergy[128]; sprintf(LegendEnergy,"Log(E_{0} / [eV]) = %d",(unsigned int)log10(energy));
  char LegendTheta[128]; sprintf(LegendTheta,"#theta = %d#circ",(unsigned int)coord[0]);
  char LegendWave[128]; sprintf(LegendWave,"#lambda = [%g,%g] cm",WaveMin,WaveMax);

  TLegend * lCherenkov = new TLegend(0.5,0.7,0.85,0.85);
  lCherenkov->SetFillColor(0); lCherenkov->SetBorderSize(0); lCherenkov->SetTextSize(0.03); lCherenkov->SetTextFont(132);
  lCherenkov->AddEntry((TObject*)0, LegendEnergy, "");
  lCherenkov->AddEntry((TObject*)0, LegendTheta, "");
  lCherenkov->AddEntry((TObject*)0, LegendWave, "");
  lCherenkov->Draw();
  cCherenkov->Update();

  cCherenkov->SaveAs("cherenkov.pdf");

  // Angular distribution
  double age[3] = {0.8,1.0,1.2};
  unsigned int color[3] = {3,2,1};
  double factor[3] = {0.25,1.0,4.0};
  TCanvas * cAngularDistribution = new TCanvas(GetObjName(),"Angular Distribution",600,600);
  cAngularDistribution->SetLogy();
  TLegend * lAngularDistribution = new TLegend(0.7,0.7,0.85,0.85);
  lAngularDistribution->SetFillColor(0); lAngularDistribution->SetBorderSize(0);
  lAngularDistribution->SetTextSize(0.035); lAngularDistribution->SetTextFont(132);
  for(unsigned int i = 0; i < 3; i++)
    {
      char entry[128]; sprintf(entry,"s = %3.1f",age[i]);
      TGraphErrors * gAngularDistribution = new TGraphErrors(angle.size());
      for(unsigned int j = 0; j < T.size(); j++)
        {
          if( depth2age(T[j],Tmax) > age[i] )
            {
              for(unsigned int k = 0; k < angle.size(); k++) gAngularDistribution->SetPoint(k,angle[k],factor[i]*AngularDistribution[j][k]);
              break;
            }
        }
      if( i == 0 )
        {
          SetPlotAttributes(gAngularDistribution,5,60,"","Angle to shower axis [deg]","1/N_{#gamma} dN_{#gamma}/d#theta [1/deg]");
          gAngularDistribution->GetYaxis()->SetRangeUser(1.e-5,1.);
          gAngularDistribution->SetLineStyle(2);
          gAngularDistribution->SetLineWidth(3);
          gAngularDistribution->SetLineColor(color[i]);
          gAngularDistribution->Draw("AL");
          TLatex * Text = new TLatex(10,5.e-4,"#times 1/4");
          Text->SetTextSize(0.03); Text->SetTextColor(color[i]);
          Text->Draw();
        }
      else
        {
          gAngularDistribution->SetLineStyle(2);
          gAngularDistribution->SetLineWidth(3);
          gAngularDistribution->SetLineColor(color[i]);
          gAngularDistribution->Draw("same L");
          if( i == 2 )
            {
              TLatex * Text = new TLatex(20,3.e-2,"#times 4");
              Text->SetTextSize(0.03); Text->SetTextColor(color[i]);
              Text->Draw();
            }
        }
      lAngularDistribution->AddEntry(gAngularDistribution, entry, "L");
      lAngularDistribution->Draw();
      cAngularDistribution->Update();
    }

  cAngularDistribution->SaveAs("angle.pdf");

  cout << "Program Finished Normally" << endl;
}
