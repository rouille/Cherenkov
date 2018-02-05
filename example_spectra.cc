#include <iostream>
#include <vector>
#include <cmath>

#include "atmosphere.h"
#include "common.h"
#include "cherenkov.h"
#include "shower.h"
#include "conversion.h"

#include <TRint.h>
#include <TRoot.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TLegend.h>

static char gObjName[1024];
static int gObjNumber = 0;
static char * GetObjName() {sprintf(gObjName,"Electron Energy Spectra %d",gObjNumber++); return gObjName;}


using namespace kPhysicalConstants;
using namespace std;



void Usage(string myName)
{
  cout << endl;
  cout << " Synopsis : " << endl;
  cout << myName << endl << endl;

  cout << " Description :" << endl;
  cout << myName << " repoduces Fig. 7 in Nerling et al. (2006) i.e plots the electron energy "
                 << " spectrum for the three shower ages s = 0.8, s = 1.0 and s = 1.2." << endl; 

  cout << endl;
  exit(0);
}



int main(int argc, char* argv[])
{
  // ROOT
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFont(30,"TITLE");

  // Electron energy spectra between 1 MeV and 10 GeV in MeV
  unsigned int size = 100;
  vector<double> energy = Bins(size,1,10000,true);
  double age[3] = {0.8,1.0,1.2};

  TCanvas * cSpectra = new TCanvas(GetObjName(),"Electron energy spectra",600,600);
  cSpectra->SetLogx();

  vector<double> spectrum_first = ElectronEnergySpectrum(energy,age[0]);
  TGraphErrors * gSpectrum_first = new TGraphErrors(size);
  for(unsigned int i = 0; i < size; i++) gSpectrum_first->SetPoint(i,energy[i],spectrum_first[i]);

  vector<double> spectrum_second = ElectronEnergySpectrum(energy,age[1]);
  TGraphErrors * gSpectrum_second = new TGraphErrors(size);
  for(unsigned int i = 0; i < size; i++) gSpectrum_second->SetPoint(i,energy[i],spectrum_second[i]);

  vector<double> spectrum_third = ElectronEnergySpectrum(energy,age[2]);
  TGraphErrors * gSpectrum_third = new TGraphErrors(size);
  for(unsigned int i = 0; i < size; i++) gSpectrum_third->SetPoint(i,energy[i],spectrum_third[i]);

  gSpectrum_first->GetYaxis()->SetRangeUser(0.,0.25);
  SetPlotAttributes(gSpectrum_first,1,10000,"","E [MeV]","1/N_{e} dN_{e}/dLn E");
  gSpectrum_first->SetMarkerStyle(kFullCircle);
  gSpectrum_first->SetMarkerColor(kGreen);
  gSpectrum_first->Draw("AP");
  cSpectra->Update();

  gSpectrum_second->SetMarkerStyle(kFullCircle);
  gSpectrum_second->SetMarkerColor(kRed);
  gSpectrum_second->Draw("same P");
  cSpectra->Update();

  gSpectrum_third->SetMarkerStyle(kFullCircle);
  gSpectrum_third->SetMarkerColor(kBlack);
  gSpectrum_third->Draw("same P");
  cSpectra->Update();

  TLegend * legend = new TLegend(0.7,0.7,0.85,0.85);
  legend->SetFillColor(0); legend->SetBorderSize(0); legend->SetTextSize(0.035); legend->SetTextFont(132);
  legend->AddEntry(gSpectrum_first, "s = 0.8", "P");
  legend->AddEntry(gSpectrum_second, "s = 1.0", "P");
  legend->AddEntry(gSpectrum_third, "s = 1.2", "P");
  legend->Draw();
  cSpectra->Update();

  cSpectra->SaveAs("spectra.pdf");

  cout << "Program Finished Normally" << endl;
}
